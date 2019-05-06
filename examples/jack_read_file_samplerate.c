/**
 * A simple wav reader, using libsndfile and libsamplerate.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <time.h>

#include <jack/jack.h>

//To read audio files
#include <sndfile.h>

//To convert between sample rates
#include <samplerate.h>

jack_port_t **output_port;
jack_client_t *client;

//sndfile stuff
SNDFILE * audio_file;
SF_INFO audio_info;
unsigned int audio_position = 0;

//samplerate stuff
#define DEFAULT_CONVERTER SRC_SINC_MEDIUM_QUALITY
float * samplerate_buff_in;
SRC_STATE * samplerate_conv;
SRC_DATA samplerate_data;

unsigned int channels = 2;

double sample_rate;

/**
 * The process callback for this JACK application is called in a
 * special realtime thread once for each audio cycle.
 *
 * This client does nothing more than copy data from its input
 * port to its output port. It will exit when stopped by 
 * the user (e.g. using Ctrl-C on a unix-ish operating system)
 */
int jack_callback (jack_nframes_t nframes, void *arg){
	jack_default_audio_sample_t **out;
	int i,j,error;
	
	//Preparing buffers
	out = (jack_default_audio_sample_t **)malloc(channels*sizeof(jack_default_audio_sample_t *));
	for(i = 0; i < channels; i++)
		out[i] = (jack_default_audio_sample_t *)jack_port_get_buffer (output_port[i], nframes);
	
	/* If the input buffer is empty, refill it. */
	if (samplerate_data.input_frames == 0){
		samplerate_data.input_frames = sf_readf_float (audio_file, samplerate_buff_in, nframes);
		samplerate_data.data_in = samplerate_buff_in; //necessary to avoid overlapping buffers
		
		/* The last read will not be a full buffer, so snd_of_input. */
		if (samplerate_data.input_frames < nframes)
			samplerate_data.end_of_input = SF_TRUE;
	}
	
	if ((error = src_process (samplerate_conv, &samplerate_data))){
		printf ("\nError : %s\n", src_strerror (error)) ;
		exit (1) ;
	}
	
	/* Terminate if done. */
	if (samplerate_data.end_of_input && samplerate_data.output_frames_gen == 0){
		printf("\nFinished reading file.\n");
		sf_close(audio_file);
		src_delete(samplerate_conv);
		jack_client_close (client);
		exit (1);
	}
	
	//Output to channels
	for (i=0;i<nframes;i++){
		for (j=0;j<channels;j++){
			if (samplerate_data.input_frames != nframes && i >= samplerate_data.input_frames){
				// Finished reading file
				// Completing the buffer with silence
				out[j][i] = 0.0;
			}else{
				out[j][i] = (jack_default_audio_sample_t)samplerate_data.data_out[i];
			}
		}
	}
	
	samplerate_data.data_in += samplerate_data.input_frames_used;
	samplerate_data.input_frames -= samplerate_data.input_frames_used;
	
	//Print Audio position
	audio_position += samplerate_data.output_frames_gen;
	printf("\rAudio file in position: %d (%0.2f secs)", audio_position, (double)audio_position/sample_rate);
	
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
	
	if(argc != 2){
		printf ("Audio File Path not provided.\n");
		exit(1);
	}
	
	//printf("%s %s\n",argv[1],argv[2]);
	char audio_file_path[100];
	strcpy(audio_file_path,argv[1]);
	printf("Trying to open audio File: %s\n",audio_file_path);
	
	audio_file = sf_open (audio_file_path,SFM_READ,&audio_info);
	if(audio_file == NULL){
		printf("%s\n",sf_strerror(NULL));
		exit(1);
	}else{
		printf("Audio file info:\n");
		printf("\tSample Rate: %d\n",audio_info.samplerate);
		printf("\tChannels: %d\n",audio_info.channels);
		SF_FORMAT_INFO audio_format_info;
		sf_command(NULL,SFC_GET_FORMAT_INFO,&audio_format_info, sizeof (SF_FORMAT_INFO));
		printf("\tFormat: %s\n",audio_format_info.name);
		
	}
	
	const char *client_name = "read_audio_file_samplerate";
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
	
	
	/* display the current sample rate. */
	sample_rate = (double)jack_get_sample_rate(client);
	printf ("Engine sample rate: %0.0f\n", sample_rate);
	
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
	
	/* create the agent output port */
	output_port = malloc(channels*sizeof(jack_port_t *));
	char port_name[50];
	for(int i = 0;i <channels;i++){
		sprintf(port_name,"output_%d",i);
		output_port[i] = jack_port_register (client, port_name, JACK_DEFAULT_AUDIO_TYPE,JackPortIsOutput, 0);
		if (output_port[i] == NULL) {
			printf("no more JACK ports available after %d\n",i);
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
		printf("no available physical playback (server input) ports.\n");
		exit (1);
	}
	// Connect the first available to our output ports
	for(int i = 0;i <channels;i++){
		if (jack_connect (client, jack_port_name (output_port[i]), serverports_names[i])) {
			printf ("cannot connect output ports\n");
			exit(1);
		}
	}
	// free serverports_names variable, we're not going to use it again
	free (serverports_names);
	
	printf ("done.\n");
	
	/* keep running until finished reading file */
	sleep(-1);
	
	jack_client_close (client);
	exit (0);
}
