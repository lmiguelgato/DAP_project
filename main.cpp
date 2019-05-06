/**
 * Final Project of Digital Audio Processing
 *
 * Objectives:
 * 1- Estimation of the angle of arrival of a sound source by generalized cross-correlation
 * 2- Sound source separation (spatial filtering) by beamforming
 *
 * Author: 		Luis Miguel Gato Diaz 		e-mail: lmiguelgato@gmail.com
 * Professor:	Caleb Rascon Estebane		e-mail: caleb.rascon@gmail.com
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// JACK: professional sound server daemon that provides real-time, 
//       low-latency connections for both audio and MIDI data between applications that use its API.
#include <jack/jack.h>

// FFTW: library for computing the discrete Fourier transform of arbitrary input size, 
//       and of both real and complex data.
#include <complex.h> 	// needs to be included before fftw3.h for compatibility
#include <fftw3.h>

// Libsndfile: library designed to allow the reading and writing of many different sampled sound file formats
//             through one standard library interface.
#include <sndfile.h>

// libsamplerate (aka Secret Rabbit Code): library for performing sample rate conversion of audio data.
#include <samplerate.h>

// Eigen: C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
#include <Eigen/Eigen>

// JACK:
jack_port_t **input_ports;
jack_port_t **output_ports;
jack_client_t *client;

// FFTW:
std::complex<double> *i_fft, *i_time, *o_fft, *o_time;
fftw_plan i_forward, o_inverse;

// Libsndfile:
SNDFILE *audio_file;
SF_INFO audio_info;
unsigned int audio_position = 0;

// libsamplerate:
#define DEFAULT_CONVERTER SRC_SINC_MEDIUM_QUALITY
float * samplerate_buff_in;
SRC_STATE * samplerate_conv;
SRC_DATA samplerate_data;

// default parameters:
double sample_rate  = 48000.0;			// default sample rate
int nframes 		= 1024;				// default number of frames per jack buffer
int window_size 	= 4096;				// default fft size (window_size must be four times nframes)

unsigned int n_out_channels = 2;		// number of output channels
unsigned int n_in_channels = 3;			// number of input channels

jack_default_audio_sample_t *hann;		// array to store the Hann window

// overlap-add registers:
jack_default_audio_sample_t **X_late;
jack_default_audio_sample_t **X_early;
jack_default_audio_sample_t **X_full;


/**
 * The process callback for this JACK application is called in a
 * special realtime thread once for each audio cycle.
 *
 * This client does nothing more than copy data from its input
 * port to its output port. It will exit when stopped by 
 * the user (e.g. using Ctrl-C on a unix-ish operating system)
 */
int jack_callback (jack_nframes_t nframes, void *arg){

	int i, j;

	jack_default_audio_sample_t **in;
	jack_default_audio_sample_t **out;

	in = (jack_default_audio_sample_t **)malloc(n_in_channels*sizeof(jack_default_audio_sample_t *));
	for(i = 0; i < n_in_channels; ++i)
		in[i] = (jack_default_audio_sample_t *)jack_port_get_buffer (input_ports[i], nframes);

	out = (jack_default_audio_sample_t **)malloc(n_out_channels*sizeof(jack_default_audio_sample_t *));
	for(i = 0; i < n_out_channels; ++i)
		out[i] = (jack_default_audio_sample_t *)jack_port_get_buffer (output_ports[i], nframes);

	// shift windows to perform overlap-add:
	for (j = 0; j < n_in_channels; ++j) {
		for (i = 0; i < nframes; ++i)
		{
			X_full[j][i] 				= X_full[j][nframes + i];
			X_full[j][nframes + i] 		= X_full[j][2*nframes + i];
			X_full[j][2*nframes + i] 	= X_full[j][3*nframes + i];
			X_full[j][3*nframes + i] 	= X_full[j][4*nframes + i];
			X_full[j][4*nframes + i] 	= X_full[j][5*nframes + i];
			X_full[j][5*nframes + i] 	= in[j][i];
		}
	}

	samplerate_data.data_in += samplerate_data.input_frames_used;
	samplerate_data.input_frames -= samplerate_data.input_frames_used;

	//Print Audio position
	audio_position += samplerate_data.output_frames_gen;
	printf("\rAudio file in position: %d (%0.2f secs)", audio_position, (double)audio_position/sample_rate);

	// ---------------------------- 1st window ------------------------------------------

	// FFT of the 1st window:
	for(i = 0; i < window_size; ++i){
		i_time[i] = X_full[i]*hann[i];
	}
	fftw_execute(i_forward);
	
	// processing of the 1st window in frequency domain:
	for(i = 0; i < window_size; ++i){
		o_fft[i] = i_fft[i];
	}
	
	// i-FFT of the 1st window:
	fftw_execute(o_inverse);
	for(i = 0; i < window_size; ++i){
		X_late[i] = real(o_time[i])/window_size; //fftw3 requires normalizing its output
	}

	// ---------------------------- 2nd window ------------------------------------------

	// FFT of the 2nd window:
	for(i = 0; i < window_size; ++i){
		i_time[i] = X_full[2*nframes+i]*hann[i];
	}
	fftw_execute(i_forward);
	
	// processing of the 2nd window in frequency domain:
	for(i = 0; i < window_size; ++i){
		o_fft[i] = i_fft[i];
	}
	
	// i-FFT of the 2nd window:
	fftw_execute(o_inverse);
	for(i = 0; i < window_size; ++i){
		X_early[i] = real(o_time[i])/window_size; //fftw3 requires normalizing its output
	}

	// --------------------------------------------------------------------------------

	// perform overlap-add:
	for (int j = 0; j < n_out_channels; ++j) {
		for (i = 0; i < nframes; ++i)
		{
			out[j][i] = X_late[i+2*nframes+nframes/2] + X_early[i+nframes/2];
		}
	}
	
	return 0;
}


/**
 * JACK calls this shutdown_callback if the server ever shuts down or
 * decides to disconnect the client.
 */
void jack_shutdown (void *arg){
	exit (1);
}


int main (int argc, char *argv[]) {

	int i;

	if(argc != 2){
		printf ("Audio File Path not provided.\n");
		exit(1);
	}
	
	char audio_file_path[100];
	strcpy(audio_file_path,argv[1]);
	printf("\nTrying to open audio File: %s\n", audio_file_path);
	
	// read audio file:
	audio_file = sf_open (audio_file_path, SFM_READ, &audio_info);
	if(audio_file == NULL){
		printf("%s\n",sf_strerror(NULL));
		exit(1);
	}else{
		printf("\nAudio file info:\n");
		printf("\tSample Rate: %d\n",audio_info.samplerate);
		printf("\tNumber of channels: %d\n",audio_info.channels);
		SF_FORMAT_INFO audio_format_info;
		sf_command(NULL,SFC_GET_FORMAT_INFO,&audio_format_info, sizeof (SF_FORMAT_INFO));
		printf("\tFormat: %d - %s\n\n",audio_info.format, audio_format_info.name);		
	}

	const char *client_name = "jack_doa_beamformer";
	jack_options_t options = JackNoStartServer;
	jack_status_t status;
	
	/* open a client connection to the JACK server */
	client = jack_client_open (client_name, options, &status);
	if (client == NULL){
		/* if connection failed, say why */
		printf ("jack_client_open() failed, status = 0x%2.0x\n", status);
		if (status & JackServerFailed) {
			printf ("Unable to connect to JACK server.\n");
		}
		exit (1);
	}
	
	/* if connection was successful, check if the name we proposed is not in use */
	if (status & JackNameNotUnique){
		client_name = jack_get_client_name(client);
		printf ("Warning: other agent with our name is running, `%s' has been assigned to us.\n", client_name);
	}
	
	/* tell the JACK server to call 'jack_callback()' whenever there is work to be done. */
	jack_set_process_callback (client, jack_callback, 0);
	
	
	/* tell the JACK server to call 'jack_shutdown()' if it ever shuts down,
	   either entirely, or if it just decides to stop calling us. */
	jack_on_shutdown (client, jack_shutdown, 0);

	// obtain here the delay from user and store it in 'delay' 
	nframes 	= (int) jack_get_buffer_size (client);
	window_size = 4*nframes;

	// initialization of internal buffers
	// - overlap-add buffers
	X_late	= (jack_default_audio_sample_t **) calloc(n_in_channels, sizeof(jack_default_audio_sample_t*));
	X_early	= (jack_default_audio_sample_t **) calloc(n_in_channels, sizeof(jack_default_audio_sample_t*));
	X_full	= (jack_default_audio_sample_t **) calloc(n_in_channels, sizeof(jack_default_audio_sample_t*));

    for (i = 0; i < n_in_channels; ++i) {
        X_late[i]	= (jack_default_audio_sample_t *) calloc(window_size, sizeof(jack_default_audio_sample_t));
		X_early[i]	= (jack_default_audio_sample_t *) calloc(window_size, sizeof(jack_default_audio_sample_t));
		X_full[i]	= (jack_default_audio_sample_t *) calloc(window_size + window_size/2, sizeof(jack_default_audio_sample_t));
    }	

	// - FFTW3 buffers
	i_fft = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size);
	i_time = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size);
	o_fft = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size);
	o_time = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size);
	
	i_forward = fftw_plan_dft_1d(window_size, reinterpret_cast<fftw_complex*>(i_time), reinterpret_cast<fftw_complex*>(i_fft), FFTW_FORWARD, FFTW_MEASURE);
	o_inverse = fftw_plan_dft_1d(window_size, reinterpret_cast<fftw_complex*>(o_fft), reinterpret_cast<fftw_complex*>(o_time), FFTW_BACKWARD, FFTW_MEASURE);
	
	sample_rate = (double) jack_get_sample_rate(client);

	// - hann window
	hann = (jack_default_audio_sample_t *) calloc(window_size, sizeof(jack_default_audio_sample_t)); 
	for(i = 0; i < window_size; ++i) {
		hann[i] = 0.5 - 0.5*cos(2.0*M_PI* ((double) i/(window_size-1)));
	}
	
	/* display the current sample rate. */
	printf ("JACK client info:\n");
	printf ("\tEngine sample rate: %d", jack_get_sample_rate (client));
	if (audio_info.samplerate != jack_get_sample_rate (client)) {
		printf ("\tWarning: sampling rate mismatch.");
	}
	printf ("\n");
	printf ("\tWindow size: %d\n\n", jack_get_buffer_size (client));	

	// creating sample rate converter:
	printf("Creating the sample rate converter...\n");
	int samplerate_error;
	samplerate_conv = src_new (DEFAULT_CONVERTER,1,&samplerate_error);
	if(samplerate_conv == NULL){
		printf("%s\n",src_strerror (samplerate_error));
		exit(1);
	}
	
	samplerate_data.src_ratio = (double)(((double)sample_rate)/((double)audio_info.samplerate));
	printf("Using samplerate ratio: %f\n", samplerate_data.src_ratio);
	if (src_is_valid_ratio (samplerate_data.src_ratio) == 0){
		printf ("Error : Sample rate change out of valid range.\n") ;
		sf_close (audio_file) ;
		exit (1) ;
	}
	samplerate_buff_in = (float *)malloc(jack_get_buffer_size(client)*sizeof(float)); //necessary to avoid overlapping buffers
	samplerate_data.data_in = samplerate_buff_in;
	samplerate_data.data_out = (float *)malloc(jack_get_buffer_size(client)*sizeof(float));
	samplerate_data.input_frames = 0;
	samplerate_data.output_frames = jack_get_buffer_size(client);
	samplerate_data.end_of_input = 0;
	
	char portname[10];
	output_ports = (jack_port_t**) malloc(n_out_channels*sizeof(jack_port_t*));
	for(i = 0; i < n_out_channels; ++i) {
		sprintf(portname, "output_%d", i+1);
		output_ports[i] = jack_port_register (client, portname, JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
		if (output_ports[i] == NULL) {
			printf("No more JACK ports available after creating port number %d\n",i);
			exit (1);
		}
	}	
	
	/* Tell the JACK server that we are ready to roll.
	   Our jack_callback() callback will start running now. */
	if (jack_activate (client)) {
		printf ("Cannot activate client.");
		exit (1);
	}
	
	printf ("Agent activated.\n");
	
	/* Connect the ports.  You can't do this before the client is
	 * activated, because we can't make connections to clients
	 * that aren't running.  Note the confusing (but necessary)
	 * orientation of the driver backend ports: playback ports are
	 * "input" to the backend, and capture ports are "output" from
	 * it.
	 */
	printf ("Connecting ports... ");

	const char **serverports_names;
	
	/* Assign our output port to a server input port*/
	// Find possible input server port names
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsInput);
	if (serverports_names == NULL) {
		printf("No available physical playback (server input) ports.\n");
		exit (1);
	}
	for(i = 0; i<n_out_channels; ++i) {
		// Connect the first available to our output port
		if (jack_connect (client, jack_port_name (output_ports[i]), serverports_names[i])) {
			printf ("Cannot connect output ports %d.\n", i);
			exit (1);
		}
	}
	// free serverports_names variable, we're not going to use it again
	free (serverports_names);
	
	
	printf ("done.\n\n");
	/* keep running until stopped by the user */
	sleep (-1);
	
	
	/* this is never reached but if the program
	   had some other way to exit besides being killed,
	   they would be important to call.
	*/
	jack_client_close (client);
	exit (0);
}