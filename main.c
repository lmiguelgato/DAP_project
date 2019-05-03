/**
 * A simple example of how to do FFT with FFTW3 and JACK.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <jack/jack.h>

// Include FFTW header
#include <complex.h> //needs to be included before fftw3.h for compatibility
#include <fftw3.h>

double complex *i_fft, *i_time, *o_fft, *o_time;
fftw_plan i_forward, o_inverse;

jack_port_t *input_port;
jack_port_t **output_ports;
jack_client_t *client;

double sample_rate;					// sample rate
int overlap_size = 512;				// default overlap_size in samples
int idx = 0;						// register index
int num_ports = 2;					// number of output ports
int windowSize = 1024;				// default window size
double *hann;						// Hann (Hanning) window
jack_default_audio_sample_t *aux;	// overlap register
jack_default_audio_sample_t *x_0;	
jack_default_audio_sample_t *x_1;	
jack_default_audio_sample_t *x_2;	


/**
 * The process callback for this JACK application is called in a
 * special realtime thread once for each audio cycle.
 *
 * This client does nothing more than copy data from its input
 * port to its output port. It will exit when stopped by 
 * the user (e.g. using Ctrl-C on a unix-ish operating system)
 */
int jack_callback (jack_nframes_t nframes, void *arg){
	jack_default_audio_sample_t *in, *out1, *out2;
	int i;
	
	in = (jack_default_audio_sample_t *)jack_port_get_buffer (input_port, nframes);
	out1 = (jack_default_audio_sample_t *)jack_port_get_buffer (output_ports[0], nframes);
	out2 = (jack_default_audio_sample_t *)jack_port_get_buffer (output_ports[1], nframes);

	for (int i = 0; i < nframes; ++i)
	{
		aux[i+nframes] 	= in[i];
		x_0[i] 			= aux[i]*hann[i];
		x_1[i]   		= aux[i+overlap_size]*hann[i];
		x_2[i]   		= aux[i+nframes]*hann[i];
		aux[i] 			= aux[i+nframes];
	}

	// Obteniendo la transformada de Fourier de este periodo
	for(i = 0; i < nframes; i++){
		i_time[i] = x_0[i];
	}
	fftw_execute(i_forward);
	
	// AquÃi podriamos hacer algo con i_fft
	for(i = 0; i < nframes; i++){
		o_fft[i] = i_fft[i];
	}
	
	// Regresando al dominio del tiempo
	fftw_execute(o_inverse);
	for(i = 0; i < nframes; i++){
		x_0[i] = creal(o_time[i])/nframes; //fftw3 requiere normalizar su salida real de esta manera
	}

	// Obteniendo la transformada de Fourier de este periodo
	for(i = 0; i < nframes; i++){
		i_time[i] = x_1[i];
	}
	fftw_execute(i_forward);
	
	// AquÃi podriamos hacer algo con i_fft
	for(i = 0; i < nframes; i++){
		o_fft[i] = i_fft[i];
	}
	
	// Regresando al dominio del tiempo
	fftw_execute(o_inverse);
	for(i = 0; i < nframes; i++){
		x_1[i] = creal(o_time[i])/nframes; //fftw3 requiere normalizar su salida real de esta manera
	}

	// Obteniendo la transformada de Fourier de este periodo
	for(i = 0; i < nframes; i++){
		i_time[i] = x_2[i];
	}
	fftw_execute(i_forward);
	
	// AquÃi podriamos hacer algo con i_fft
	for(i = 0; i < nframes; i++){
		o_fft[i] = i_fft[i];
	}
	
	// Regresando al dominio del tiempo
	fftw_execute(o_inverse);
	for(i = 0; i < nframes; i++){
		x_2[i] = creal(o_time[i])/nframes; //fftw3 requiere normalizar su salida real de esta manera
	}

	for (int i = 0; i < overlap_size; ++i)
	{
		out1[i] = x_0[i+overlap_size] + x_1[i];
		out2[i] = out1[i];
		out1[i+overlap_size] = x_1[i+overlap_size] + x_2[i];
		out2[i+overlap_size] = out1[i+overlap_size];
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

	/*
	printf("Introduce delay in samples: ");
	scanf("%d", &delay);
	time_delay = ((float)delay)*1000.0/((float) jack_get_sample_rate (client));
	printf("\nThe delay was set to %f ms\n\n", time_delay);
	*/

	// obtain here the delay from user and store it in 'delay' 
	windowSize = (int) jack_get_buffer_size (client);
	overlap_size = (int) windowSize / 2;
	aux = (jack_default_audio_sample_t *) calloc(windowSize*2, sizeof(jack_default_audio_sample_t)); // calloc initializes with zeros
	x_0 = (jack_default_audio_sample_t *) calloc(windowSize, sizeof(jack_default_audio_sample_t)); // calloc initializes with zeros
	x_1 = (jack_default_audio_sample_t *) calloc(windowSize, sizeof(jack_default_audio_sample_t)); // calloc initializes with zeros
	x_2 = (jack_default_audio_sample_t *) calloc(windowSize, sizeof(jack_default_audio_sample_t)); // calloc initializes with zeros
	hann = (double *) calloc(windowSize, sizeof(double)); // calloc initializes with zeros

	for(i = 0; i < windowSize; ++i) {
		hann[i] = 0.5 - 0.5*cos(2.0*M_PI* ((double) i/windowSize));
	}
	
	/* display the current sample rate. */
	printf ("Sample rate: %d\n", jack_get_sample_rate (client));
	printf ("Window size: %d\n", jack_get_buffer_size (client));
	sample_rate = (double)jack_get_sample_rate(client);
	int nframes = jack_get_buffer_size (client);
	
	//preparing FFTW3 buffers
	i_fft = (double complex *) fftw_malloc(sizeof(double complex) * nframes);
	i_time = (double complex *) fftw_malloc(sizeof(double complex) * nframes);
	o_fft = (double complex *) fftw_malloc(sizeof(double complex) * nframes);
	o_time = (double complex *) fftw_malloc(sizeof(double complex) * nframes);
	
	i_forward = fftw_plan_dft_1d(nframes, i_time, i_fft , FFTW_FORWARD, FFTW_MEASURE);
	o_inverse = fftw_plan_dft_1d(nframes, o_fft , o_time, FFTW_BACKWARD, FFTW_MEASURE);
	
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