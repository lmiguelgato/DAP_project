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

// JACK:
jack_port_t *input_port;
jack_port_t **output_ports;
jack_client_t *client;

// FFTW:
double complex *i_fft, *i_time, *o_fft, *o_time;
fftw_plan i_forward, o_inverse;

// Libsndfile:
SNDFILE *audio_file;
SF_INFO audio_info;
unsigned int audio_position = 0;

// default parameters:
double sample_rate  = 48000.0;			// default sample rate
int nframes 		= 1024;				// default number of frames per jack buffer
int window_size 	= 4096;				// default fft size (window_size must be four times nframes)
int num_ports 		= 2;				// number of output ports

jack_default_audio_sample_t *hann;		// array to store the Hann window

// overlap-add registers:
jack_default_audio_sample_t *X_late;
jack_default_audio_sample_t *X_early;
jack_default_audio_sample_t *X_full;


/**
 * The process callback for this JACK application is called in a
 * special realtime thread once for each audio cycle.
 *
 * This client does nothing more than copy data from its input
 * port to its output port. It will exit when stopped by 
 * the user (e.g. using Ctrl-C on a unix-ish operating system)
 */
int jack_callback (jack_nframes_t nframes, void *arg){

	int i;
	jack_default_audio_sample_t *in, *out1, *out2;
	
	in 		= (jack_default_audio_sample_t *)jack_port_get_buffer (input_port, 		nframes);
	out1 	= (jack_default_audio_sample_t *)jack_port_get_buffer (output_ports[0], nframes);
	out2 	= (jack_default_audio_sample_t *)jack_port_get_buffer (output_ports[1], nframes);

	// shift windows to perform overlap-add:
	for (i = 0; i < nframes; ++i)
	{
		X_full[i] 				= X_full[nframes + i];
		X_full[nframes + i] 	= X_full[2*nframes + i];
		X_full[2*nframes + i] 	= X_full[3*nframes + i];
		X_full[3*nframes + i] 	= X_full[4*nframes + i];
		X_full[4*nframes + i] 	= X_full[5*nframes + i];
		X_full[5*nframes + i] 	= in[i];
	}

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
		X_late[i] = creal(o_time[i])/window_size; //fftw3 requires normalizing its output
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
		X_early[i] = creal(o_time[i])/window_size; //fftw3 requires normalizing its output
	}

	// --------------------------------------------------------------------------------

	// perform overlap-add:
	for (i = 0; i < nframes; ++i)
	{
		out1[i] = X_late[i+2*nframes+nframes/2] + X_early[i+nframes/2];

		// turn stereo to mono:
		out2[i] = out1[i];
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

	const char *client_name = "jack_fft";
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
	X_late	= (jack_default_audio_sample_t *) calloc(window_size, sizeof(jack_default_audio_sample_t));
	X_early	= (jack_default_audio_sample_t *) calloc(window_size, sizeof(jack_default_audio_sample_t));
	X_full	= (jack_default_audio_sample_t *) calloc(window_size + window_size/2, sizeof(jack_default_audio_sample_t));

	// - FFTW3 buffers
	i_fft 	= (double complex *) fftw_malloc(sizeof(double complex) * window_size);
	i_time 	= (double complex *) fftw_malloc(sizeof(double complex) * window_size);
	o_fft 	= (double complex *) fftw_malloc(sizeof(double complex) * window_size);
	o_time 	= (double complex *) fftw_malloc(sizeof(double complex) * window_size);

	sample_rate = (double) jack_get_sample_rate(client);	
	
	i_forward = fftw_plan_dft_1d(window_size, i_time, i_fft , FFTW_FORWARD, FFTW_MEASURE);
	o_inverse = fftw_plan_dft_1d(window_size, o_fft , o_time, FFTW_BACKWARD, FFTW_MEASURE);

	// - hann window
	hann = (jack_default_audio_sample_t *) calloc(window_size, sizeof(jack_default_audio_sample_t)); 
	for(i = 0; i < window_size; ++i) {
		hann[i] = 0.5 - 0.5*cos(2.0*M_PI* ((double) i/(window_size-1)));
	}
	
	/* display the current sample rate. */
	printf ("Sample rate: %d\n", jack_get_sample_rate (client));
	printf ("Window size: %d\n", jack_get_buffer_size (client));	
	
	/* create the agent input port */
	input_port = jack_port_register (client, "input", JACK_DEFAULT_AUDIO_TYPE,JackPortIsInput, 0);
	
	char portname[10];
	output_ports = (jack_port_t**) malloc(num_ports*sizeof(jack_port_t*));
	for(i = 0; i<num_ports; ++i) {
		sprintf(portname, "output_%d", i+1);
		output_ports[i] = jack_port_register (client, portname, JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);
	}
	
	/* check that both ports were created succesfully */
	if (input_port == NULL) {
		printf("Could not create input agent ports. Have we reached the maximum amount of JACK agent ports?\n");
		exit (1);
	}

	for(i = 0; i<num_ports; ++i) {
		if (output_ports[i] == NULL) {
			printf("Could not create output agent port number %d. Have we reached the maximum amount of JACK agent ports?\n", i);
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
	 
	/* Assign our input port to a server output port*/
	// Find possible output server port names
	const char **serverports_names;
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsOutput);
	if (serverports_names == NULL) {
		printf("No available physical capture (server output) ports.\n");
		exit (1);
	}
	// Connect the first available to our input port
	if (jack_connect (client, serverports_names[0], jack_port_name (input_port))) {
		printf("Cannot connect input port.\n");
		exit (1);
	}
	// free serverports_names variable for reuse in next part of the code
	free (serverports_names);
	
	
	/* Assign our output port to a server input port*/
	// Find possible input server port names
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsInput);
	if (serverports_names == NULL) {
		printf("No available physical playback (server input) ports.\n");
		exit (1);
	}
	for(i = 0; i<num_ports; ++i) {
		// Connect the first available to our output port
		if (jack_connect (client, jack_port_name (output_ports[i]), serverports_names[i])) {
			printf ("Cannot connect output ports %d.\n", i);
			exit (1);
		}
	}
	// free serverports_names variable, we're not going to use it again
	free (serverports_names);
	
	
	printf ("done.\n");
	/* keep running until stopped by the user */
	sleep (-1);
	
	
	/* this is never reached but if the program
	   had some other way to exit besides being killed,
	   they would be important to call.
	*/
	jack_client_close (client);
	exit (0);
}