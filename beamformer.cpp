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
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

// JACK: professional sound server daemon that provides real-time, 
//       low-latency connections for both audio and MIDI data between applications that use its API.
#include <jack/jack.h>

// FFTW: library for computing the discrete Fourier transform of arbitrary input size, 
//       and of both real and complex data.
#include <complex.h> 	// needs to be included before fftw3.h for compatibility
#include <cmath>
#include <fftw3.h>

// Eigen: C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
#include <Eigen/Eigen>

// Libsndfile: library designed to allow the reading and writing of many different sampled sound file formats
//             through one standard library interface.
#include <sndfile.h>

#define RAD2DEG 57.295779513082323f		// useful to convert from radians to degrees			// verbose

double *DOA_valid;						// store valid DOAs

// JACK:
jack_port_t **input_ports;
jack_port_t **output_ports;
jack_client_t *client;

// Libsndfile:
SNDFILE ** audio_file;
SF_INFO audio_info;
float *write_buffer;

// FFTW:
std::complex<double> *i_fft_4N, *i_time_4N, *o_fft_4N, *o_time_4N;
std::complex<double> *i_fft_2N, *i_time_2N, *o_fft_2N, *o_time_2N;
fftw_plan i_forward_4N, o_inverse_4N, i_forward_2N, o_inverse_2N;
std::complex<double> ig(0.0, 1.0); 		// imaginary unit

// default parameters:
double sample_rate  = 48000.0;			// default sample rate [Hz]
int nframes 		= 1024;				// default number of frames per jack buffer
int window_size, window_size_2, nframes_2;	
double mic_separation = 0.1;			// default microphone separation [meters]
double c = 343.364;						// default sound speed [m/s]
int n_sources = 1; 						// default number of sources to be detected
double fRes; 							// frequency resolution

int state = 0;							// beamformer state

unsigned int n_out_channels = 2;		// number of output channels
unsigned int n_in_channels = 3;			// number of input channels

jack_default_audio_sample_t *hann;		// array to store the Hann window

// overlap-add registers:
jack_default_audio_sample_t **X_full;	// store the 6 last buffers of 'nframes' samples
jack_default_audio_sample_t **X_late;	// store the 4 latest buffers from X_full
jack_default_audio_sample_t **X_early;	// store the 4 earliest buffers from X_full

// GCC registers:
std::complex<double> **X_gcc;			// store GCC result (length 2 times 'nframes')
std::complex<double> *Aux_gcc;			// store axuliary GCC result (length 2 times 'nframes')

/**
 * The process callback for this JACK application is called in a
 * special realtime thread once for each audio cycle.
 *
 * This client does nothing more than copy data from its input
 * port to its output port. It will exit when stopped by 
 * the user (e.g. using Ctrl-C on a unix-ish operating system)
 */
int jack_callback (jack_nframes_t nframes, void *arg){

	int i, j, k;

	jack_default_audio_sample_t **in;
	jack_default_audio_sample_t **out;

	in = (jack_default_audio_sample_t **)malloc(n_in_channels*sizeof(jack_default_audio_sample_t *));
	for(i = 0; i < n_in_channels; ++i)
		in[i] = (jack_default_audio_sample_t *)jack_port_get_buffer (input_ports[i], nframes);

	out = (jack_default_audio_sample_t **)malloc(n_out_channels*sizeof(jack_default_audio_sample_t *));
	for(i = 0; i < n_out_channels; ++i)
		out[i] = (jack_default_audio_sample_t *)jack_port_get_buffer (output_ports[i], nframes);

	for (j = 0; j < n_in_channels; ++j) {

		// shift-register (useful for overlap-add method):
		for (i = 0; i < nframes; ++i)
		{
			X_full[j][i] 						= X_full[j][nframes + i];
			X_full[j][nframes + i] 				= X_full[j][window_size_2 + i];
			X_full[j][window_size_2 + i] 		= X_full[j][window_size-nframes + i];
			X_full[j][window_size-nframes + i] 	= X_full[j][window_size + i];
			X_full[j][window_size + i] 			= X_full[j][window_size+nframes + i];
			X_full[j][window_size+nframes + i] 	= in[j][i];	
		}
	}

	int write_count;
	int delay[2];

	// filter all
	for (int s2f = 1; s2f <= n_sources; ++s2f)
	{		
		delay[0] = (int) (mic_separation*sin(DOA_valid[s2f-1]/RAD2DEG)/c*sample_rate);
		delay[1] = (int) (mic_separation*sin((120.0-DOA_valid[s2f-1])/RAD2DEG)/c*sample_rate);

		for (k = 1; k < n_in_channels; ++k)
		{
			// ---------------------------- 1st window ------------------------------------------

			// FFT of the 1st window:
			for(i = 0; i < window_size; i++){
				i_time_4N[i] = X_full[k][i]*hann[i];
			}
			fftw_execute(i_forward_4N);
			
			// delay the 1st window in frequency domain:
			for(i = 0; i < window_size; i++){
				o_fft_4N[i] = i_fft_4N[i]*exp(ig*fRes*((double) i*delay[k-1]));
			}
			
			// i-FFT of the 1st window:
			fftw_execute(o_inverse_4N);
			for(i = 0; i < window_size; i++){
				X_late[k][i] = real(o_time_4N[i])/window_size; //fftw3 requires normalizing its output
			}

			// ---------------------------- 2nd window ------------------------------------------

			// FFT of the 2nd window:
			for(i = 0; i < window_size; i++){
				i_time_4N[i] = X_full[k][window_size_2+i]*hann[i];
			}
			fftw_execute(i_forward_4N);
			
			// delay the 2nd window in frequency domain:
			for(i = 0; i < window_size; i++){
				o_fft_4N[i] = i_fft_4N[i]*exp(ig*fRes*((double) i*delay[k-1]));
			}
			
			// i-FFT of the 2nd window:
			fftw_execute(o_inverse_4N);
			for(i = 0; i < window_size; i++){
				X_early[k][i] = real(o_time_4N[i])/window_size; //fftw3 requires normalizing its output
			}

			// perform overlap-add:		
			for (i = 0; i < nframes; ++i) {
				write_buffer[i] = X_full[0][i+window_size_2+nframes_2];;
				for (k = 1; k < n_in_channels; ++k)
				{
					write_buffer[i] += X_late[k][i+window_size_2+nframes_2] + X_early[k][i+nframes_2];
				}
				write_buffer[i] /= 3.0;
			}

			if (in[0][nframes_2] == 0.0f && in[1][nframes_2+10] == 0.0f && in[2][nframes_2-10] == 0.0f && write_buffer[nframes_2] == 0.0f)
			{
				if(state == 1) exit(1);
				/*for (i = 0; i < nframes; ++i) {
					write_buffer[i] = 0.0;
				}
				write_count = sf_write_float(audio_file[s2f-1],write_buffer,nframes);*/
			} else {
				if(state == 0) state = 1;
				write_count = sf_write_float(audio_file[s2f-1],write_buffer,nframes);
				//Check for writing error
				if(write_count != nframes){
					printf("\nEncountered I/O error. Exiting.\n");
					sf_close(audio_file[s2f-1]);
					jack_client_close (client);
					exit (1);
				}
			}

			
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
		printf ("Usage:\ndap_project s\ns: path to the beamformer initialization file.\n");
		exit(1);
	}

	string input_file_path = string(argv[1]);

	ifstream InputFile;

	InputFile.open(input_file_path);

	if (!InputFile.is_open()) {
		printf ("Error opening file.\n");
		exit(1);
	}

	string str;

	if(std::getline(InputFile, str)) {		
		mic_separation = atof (str.c_str()); 	
		printf ("\n\n");
		printf ("Microphone separation: %f meters.\n", mic_separation);
	}

	printf ("DOAs: ");
	int countSources = 0;
	DOA_valid	= (double *) calloc(10, sizeof(double));
	std::string doa; 
	std::istringstream iss(doa); 
	while (std::getline(InputFile, doa, ',')) 
    { 	
    	++countSources;
        DOA_valid[countSources-1] = atof(doa.c_str());
        
        printf ("\t%f", DOA_valid[countSources-1]);
    }
    printf ("\n\n");


    n_sources = countSources;

	system("mkdir -p output");

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
	nframes_2   = nframes/2;
	window_size = 4*nframes;
	window_size_2 = 2*nframes;

	fRes = 2.0*M_PI/window_size;

	// initialization of internal buffers
	// - overlap-add buffers
	X_late		= (jack_default_audio_sample_t **) calloc(n_in_channels, sizeof(jack_default_audio_sample_t*));
	X_early		= (jack_default_audio_sample_t **) calloc(n_in_channels, sizeof(jack_default_audio_sample_t*));
	X_full		= (jack_default_audio_sample_t **) calloc(n_in_channels, sizeof(jack_default_audio_sample_t*));
	
	write_buffer = (float *) calloc(nframes, sizeof(float));

	for (i = 0; i < n_in_channels; ++i) {
        X_late[i]	= (jack_default_audio_sample_t *) calloc(window_size, sizeof(jack_default_audio_sample_t));
		X_early[i]	= (jack_default_audio_sample_t *) calloc(window_size, sizeof(jack_default_audio_sample_t));
		X_full[i]	= (jack_default_audio_sample_t *) calloc(window_size + window_size_2, sizeof(jack_default_audio_sample_t));
    }	

	// - FFTW3 buffers
	i_fft_4N = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size);
	i_time_4N = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size);
	o_fft_4N = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size);
	o_time_4N = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size);
	
	i_forward_4N = fftw_plan_dft_1d(window_size, reinterpret_cast<fftw_complex*>(i_time_4N), reinterpret_cast<fftw_complex*>(i_fft_4N), FFTW_FORWARD, FFTW_MEASURE);
	o_inverse_4N = fftw_plan_dft_1d(window_size, reinterpret_cast<fftw_complex*>(o_fft_4N), reinterpret_cast<fftw_complex*>(o_time_4N), FFTW_BACKWARD, FFTW_MEASURE);
	
	sample_rate = (double) jack_get_sample_rate(client);

	// - hann window
	hann = (jack_default_audio_sample_t *) calloc(window_size, sizeof(jack_default_audio_sample_t)); 
	for(i = 0; i < window_size; ++i) {
		hann[i] = 0.5 - 0.5*cos(2.0*M_PI* ((double) i/(window_size-1)));
	}
	
	/* display the current sample rate. */
	printf ("JACK client info:\n");
	printf ("\tEngine sample rate: %d\n", jack_get_sample_rate (client));
	printf ("\tWindow size: %d\n\n", jack_get_buffer_size (client));
	
	char portname[13];
	output_ports = (jack_port_t**) malloc(n_out_channels*sizeof(jack_port_t*));
	for(i = 0; i < n_out_channels; ++i) {
		sprintf(portname, "speaker_%d", i+1);
		output_ports[i] = jack_port_register (client, portname, JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
		if (output_ports[i] == NULL) {
			printf("No more JACK ports available after creating output port number %d\n",i);
			exit (1);
		}
	}	

	input_ports = (jack_port_t**) malloc(n_in_channels*sizeof(jack_port_t*));
	for(i = 0; i < n_in_channels; ++i) {
		sprintf(portname, "wav_mic%d", i+1);
		input_ports[i] = jack_port_register (client, portname, JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
		if (input_ports[i] == NULL) {
			printf("No more JACK ports available after creating input port number %d\n",i);
			exit (1);
		}
	}	

	audio_file = (SNDFILE **) calloc(n_sources, sizeof(SNDFILE*));

	char audio_file_path[40];
	for (i = 0; i < n_sources; ++i)
	{
		printf("Trying to open audio File: audio_%d_(%d).wav",n_sources,i+1);
		audio_info.samplerate = sample_rate;
		audio_info.channels = 1;
		audio_info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_32;
		sprintf(audio_file_path, "output/audio_%d_(%d).wav", n_sources, i+1);
		audio_file[i] = sf_open (audio_file_path,SFM_WRITE,&audio_info);
		if(audio_file[i] == NULL){
			printf("%s\n",sf_strerror(NULL));
			exit(1);
		}else{
			printf("Audio file info:\n");
			printf("\tSample Rate: %d\n",audio_info.samplerate);
			printf("\tChannels: %d\n",audio_info.channels);
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
			printf ("Cannot connect input port number %d.\n", i);
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