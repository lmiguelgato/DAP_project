/*
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
#include <math.h>
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

// Useful functions developed by myself:
#include "max.h"
#include "unwrap.h"
#include "phat.h"
#include "angleManipulation.h"
#include "kmeans.h"
#include "kalman.h"

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

#define RAD2DEG 57.295779513082323f		// useful to convert from radians to degrees
#define GCC_STYLE 4						// 1: GCC, 2:GCC (frequency restrained), 3:GCC-PHAT, 4:GCC-PHAT (frequency restrained)
#define GCC_TH 100.0f					// correlation threshold (to avoid false alarms)
#define REDUNDANCY_TH 20.0f				// redundancy threshold (for DOA estimation)
#define DYNAMIC_GCC_TH 1				// enable a dynamic GCC threshold (0: disabled, 1: mean peak values, 2: max peak values)
#define MOVING_AVERAGE 1				// enable a moving average on kmeans centroids (0: disabled, 1: finite memory, 2: infinite memory)
#define MOVING_FACTOR 1					// allow variations in DOA if the sources are moving (how many times the standard deviation)
#define MEMORY_FACTOR 5				// memory of the k-means algorithm
#define VERBOSE false					// display additional information

// JACK:
jack_port_t **input_ports;
jack_port_t **output_ports;
jack_client_t *client;

// Libsndfile:
SNDFILE * saudio_file;
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
unsigned int window_size, window_size_2, nframes_2;	
double mic_separation = 0.18;			// default microphone separation [meters]
double c = 343.364;						// default sound speed [m/s]
double dt_max, N_max;					// maximum delay between microphones [s, samples]
double doa;								// direction of arrival
double mean_doa = 0.0;					// average direction of arrival
double std_doa = 90.0;					// standard deviation of the direction of arrival
double std_cum = 90.0;					// standard deviation of the direction of arrival (cummulative)
unsigned int n_sources = 1; 			// default number of sources to be detected
int source2filter = 0;					// default source to filter
double gcc_th = GCC_TH;					// default GCC threshold
double fRes; 							// frequency resolution

double f_min = 1000.0;					// minimum frequency of the desired signal [Hz]
double f_max = 4000.0;					// maximum frequency of the desired signal [Hz]
int kmin, kmax;							// discrete minimum and maximum frequencies of the desired signal

unsigned int n_out_channels = 2;		// number of output channels
unsigned int n_in_channels = 3;			// number of input channels

int state = 0;							// beamformer state

jack_default_audio_sample_t *hann;		// array to store the Hann window

// overlap-add registers:
jack_default_audio_sample_t **X_full;	// store the 6 last buffers of 'nframes' samples
jack_default_audio_sample_t **X_late;	// store the 4 latest buffers from X_full
jack_default_audio_sample_t **X_early;	// store the 4 earliest buffers from X_full

// GCC registers:
std::complex<double> **X_gcc;			// store GCC result (length 2 times 'nframes')
std::complex<double> *Aux_gcc;			// store axuliary GCC result (length 2 times 'nframes')

// DOA registers:
double *DOA_hist;						// store DOA history
unsigned int *DOA_class;				// store classification of DOA history
double *DOA_kmean;						// store DOA mean for each source (kmeans algorithm)
double *DOA_mean;						// store DOA mean for each source
double *DOA_stdev;						// store DOA standard deviation for each source
double *DOA_valid;						// store valid DOAs
unsigned int *counter;					// store number of detections for each source
int    *dcounter;						// number of valid DOAs
int    *ecounter;						// number of smoothed valid DOAs
int    icounter = 0;					// number of detections
int    ccounter = 0;					// number of cycles
int    hist_length;						// number of cycles
bool   firstDetection = true;

double **kalmanState;
double **covMatrix;

ofstream outputFile;					// save results for data analysis
ofstream outputKalman;					// save results for data analysis

/**
 * The process callback for this JACK application is called in a
 * special realtime thread once for each audio cycle.
 *
 * This client does nothing more than copy data from its input
 * port to its output port. It will exit when stopped by 
 * the user (e.g. using Ctrl-C on a unix-ish operating system)
 */
int jack_callback (jack_nframes_t nframes, void *arg){

	unsigned int i, j, k;

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

		// cross-correlation in four steps:
		// 1- zero padding:
		for (i = 0; i < nframes; ++i)
			i_time_2N[i] = in[j][i];

		for (i = nframes; i < window_size_2; ++i)
			i_time_2N[i] = 0.0;

		// 2- apply FFT:
		fftw_execute(i_forward_2N);

		for (i = 0; i < window_size_2; ++i)
			X_gcc[j][i] = i_fft_2N[i];
	}

	// 3- multiply pairs of FFTs (time reversing one of them), and
	// 4- apply iFFT
	phat (o_fft_2N, X_gcc[0], X_gcc[1], window_size_2, kmin, kmax, GCC_STYLE);
	fftw_execute(o_inverse_2N);

	for (i = 0; i < window_size_2; ++i)
		Aux_gcc[i] = o_time_2N[i];

	phat (o_fft_2N, X_gcc[1], X_gcc[2], window_size_2, kmin, kmax, GCC_STYLE);
	fftw_execute(o_inverse_2N);

	for (i = 0; i < window_size_2; ++i) {
		X_gcc[1][i] = Aux_gcc[i];
		Aux_gcc[i] = o_time_2N[i];
	}

	phat (o_fft_2N, X_gcc[2], X_gcc[0], window_size_2, kmin, kmax, GCC_STYLE);
	fftw_execute(o_inverse_2N);

	for (i = 0; i < window_size_2; ++i) {
		X_gcc[2][i] = Aux_gcc[i];
		X_gcc[0][i] = o_time_2N[i];
	}

	double max_val12, max_val23, max_val31;

	// find maximum of the cross-correlations, and estimate DOA:
	double theta[6] = {asin(  unwrap(  max(X_gcc[1], window_size_2, N_max, &max_val12), nframes, N_max  )/N_max  )*RAD2DEG,
					   asin(  unwrap(  max(X_gcc[2], window_size_2, N_max, &max_val23), nframes, N_max  )/N_max  )*RAD2DEG,
					   asin(  unwrap(  max(X_gcc[0], window_size_2, N_max, &max_val31), nframes, N_max  )/N_max  )*RAD2DEG,
					   0.0,
					   0.0,
					   0.0};

	angleTranslation(theta);	// use a coherent reference to measure DOA

	double thetaRedundant[3] = {0.0, 0.0, 0.0};

	doa = angleRedundancy (theta, thetaRedundant, REDUNDANCY_TH);

	double max_mean, max_max;

	switch (DYNAMIC_GCC_TH) { //enable a dynamic GCC threshold
		case 1:	// mean peak values
			max_mean = (max_val12 + max_val23 + max_val31)/3;

			if (max_mean > GCC_TH) {
				++ccounter;
				//gcc_th = (gcc_th*ccounter + max_mean)/(ccounter+1);
				gcc_th = max_mean;
			}
		break;

		case 2:	// max peak values
			max_max = 0.0;

			if (max_val12 > max_max)
				max_max = max_val12;
			if (max_val23 > max_max)
				max_max = max_val23;
			if (max_val31 > max_max)
				max_max = max_val31;

			if (max_max > 1.0) {
				++ccounter;
				gcc_th = (gcc_th*(ccounter-1) + max_max)/ccounter;
			}
		break;

		default:
			break;
	}

	if (VERBOSE)
	{
		printf("theta1 = [%1.5f, %1.5f];\ttheta2 = [%1.5f, %1.5f];\ttheta3 = [%1.5f, %1.5f]\n", theta[0], theta[3], theta[1], theta[4], theta[2], theta[5]);
		printf("val1 = %1.5f;\tval2 = %1.5f;\tval3 = %1.5f\n", max_val12, max_val23, max_val31);
		printf("*gcc_th = %1.5f\n", gcc_th);
		printf("thetaR = [%1.5f, %1.5f, %1.5f]\n", thetaRedundant[0], thetaRedundant[1], thetaRedundant[2]);
	}

	if (doa != 181.0 && max_val12 > 0.9*gcc_th && max_val23 > 0.9*gcc_th && max_val31 > 0.9*gcc_th) { 	// are these DOAs valid?
		outputFile << icounter << ':' << ' ';
		DOA_hist[icounter%hist_length] = doa; 	// save into shift register
		++icounter;

		if (icounter >= hist_length) {	// is the shift register full?

			if (icounter == hist_length && firstDetection) {
				firstDetection = false;

				DOA_kmean[0] = mean(DOA_hist, hist_length);
				
				for (i = 1; i < n_sources; ++i)
				{
					DOA_kmean[i] = doa + 360.0/n_sources*i;		

					if (DOA_kmean[i] > 180.0) {
						DOA_kmean[i] -= 360.0;
					}
				}

				// initialize Kalman state:
				double initialThisState[2];
				for (i = 0; i < n_sources; ++i) {
					angle2state(DOA_kmean[i], initialThisState);

					kalmanState[0][i] = initialThisState[0];
					kalmanState[1][i] = initialThisState[1];

					printf("kalman initialization[%d]: %1.1f\n", i, DOA_kmean[i]);
				}
			}

			// group similar DOAs into clusters using the k-means algorithm:
			kmeans (DOA_hist, DOA_class, DOA_kmean, counter, n_sources, hist_length);

			double measurement[2];
			double state[4];
			double cov[4][4];

			for (i = 0; i < n_sources; ++i)
			{
				//if (counter[i] > 0) {	// any DOA in this cluster?

					if (DOA_class[icounter%hist_length] == i) {
						angle2state (DOA_hist[icounter%hist_length], measurement);

						state[0] = kalmanState[0][i];	state[1] = kalmanState[1][i];	state[2] = kalmanState[2][i];	state[3] = kalmanState[3][i];

						for (j = 0; j < 4; ++j) {
							cov[j][0] = covMatrix[4*i+j][0];	cov[j][1] = covMatrix[4*i+j][1];	cov[j][2] = covMatrix[4*i+j][2];	cov[j][3] = covMatrix[4*i+j][3];
						}

						kalman (measurement, state, cov);

						outputKalman << DOA_class[icounter%hist_length] << ' ';
						outputKalman << setprecision(2) << DOA_hist[icounter%hist_length];			// save results into text file
						outputKalman << ' ';
						outputKalman << setprecision(2) << state2angle (state);			// save results into text file
						outputKalman << endl;

						kalmanState[0][i] = state[0];	kalmanState[1][i] = state[1];	kalmanState[2][i] = state[2];	kalmanState[3][i] = state[3];
						for (j = 0; j < 4; ++j) {
							covMatrix[4*i+j][0] = cov[j][0];	covMatrix[4*i+j][1] = cov[j][1];	covMatrix[4*i+j][2] = cov[j][2];	covMatrix[4*i+j][3] = cov[j][3];
						}
					}
				//}	
			}

			for (i = 0; i < n_sources; ++i)
			{
				if (counter[i] > 0) {	// any DOA in this cluster?
					++dcounter[i];

					if (MOVING_AVERAGE == 2) { //enable moving average
						DOA_mean[i] = (DOA_mean[i]*(dcounter[i]-1) + DOA_kmean[i])/dcounter[i];		// moving average
						DOA_stdev[i] += pow(DOA_kmean[i]-DOA_mean[i], 2);							// standard deviation
					} else {
						DOA_mean[i] = (DOA_mean[i] + DOA_kmean[i])/2.0;		// moving average
						DOA_stdev[i] += pow(DOA_kmean[i]-DOA_mean[i], 2);							// standard deviation
					}

					if (abs(DOA_kmean[i]-DOA_mean[i]) < MOVING_FACTOR*sqrt(DOA_stdev[i]/dcounter[i])) {		// avoid outsiders
						++ecounter[i];
						if (MOVING_AVERAGE == 0) {
							DOA_valid[i] = DOA_kmean[i];
						} else {
							DOA_valid[i] = DOA_mean[i];
						}

						printf("DOA[%d] = %1.1f\n", i, DOA_valid[i]);
						outputFile << setprecision(2) << DOA_valid[i];			// save results into text file
					}

				} else {
						outputFile << 181;			// save results into text file
				}
				outputFile << ", ";
			}
		} 
		else {
			for (i = 0; i < n_sources; ++i)
			{
				outputFile << 181 << ',' << ' ';
			}
		}
		outputFile << endl;
	}

	unsigned int write_count;
	int delay[2];

	if (source2filter != 0) {
		if (source2filter != -1) {
			if (ecounter[source2filter-1] > 0) 	// if at least one valid DOA was found of the target source, apply beamforming
			{
				delay[0] = (int) (mic_separation*sin(DOA_valid[source2filter-1]/RAD2DEG)/c*sample_rate);
				delay[1] = (int) (mic_separation*sin((120.0-DOA_valid[source2filter-1])/RAD2DEG)/c*sample_rate);

				//printf("doa = %1.5f\n", DOA_valid[source2filter-1]);

				//----  BEAMFORMING:

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
				}
				
				max_max = 0.0f;

					if (max_val12 > max_max)
						max_max = max_val12;
					if (max_val23 > max_max)
						max_max = max_val23;
					if (max_val31 > max_max)
						max_max = max_val31;

					if (max_max == 0.0f) {
						if(state == 1) exit(1);
					}
					else {
						if(state == 0) state = 1;

						// perform overlap-add:		
						for (i = 0; i < nframes; ++i) {
							out[0][i] = X_full[0][i+window_size_2+nframes_2];
							for (k = 1; k < n_in_channels; ++k)
							{
								out[0][i] += X_late[k][i+window_size_2+nframes_2] + X_early[k][i+nframes_2];
							}
							out[0][i] /= 3.0;
							out[1][i] = 0.0;

							write_buffer[i] = out[0][i];
						}

						write_count = sf_write_float(saudio_file,write_buffer,nframes);

						//Check for writing error
						if(write_count != nframes){
							printf("\nEncountered I/O error. Exiting.\n");
							sf_close(saudio_file);
							jack_client_close (client);
							exit (1);
					}
				}
				
			} else {	// don't start writing if no audio is being played
				max_max = 0.0;

				if (max_val12 > max_max)
					max_max = max_val12;
				if (max_val23 > max_max)
					max_max = max_val23;
				if (max_val31 > max_max)
					max_max = max_val31;

				if (max_max == 0.0f) {
					if(state == 1) exit(1);
				} else {
					if(state == 0) state = 1;

					for (i = 0; i < nframes; ++i) {
						write_buffer[i] = X_full[0][i+window_size_2+nframes_2];
					}
					write_count = sf_write_float(saudio_file,write_buffer,nframes);

					//Check for writing error
					if(write_count != nframes){
						printf("\nEncountered I/O error. Exiting.\n");
						sf_close(saudio_file);
						jack_client_close (client);
						exit (1);
					}
				}
			}
		} else {
			max_max = 0.0f;

			if (max_val12 > max_max)
				max_max = max_val12;
			if (max_val23 > max_max)
				max_max = max_val23;
			if (max_val31 > max_max)
				max_max = max_val31;

			if (max_max == 0.0f) {
				if(state == 1) exit(1);
			}
			else {
				if(state == 0) state = 1;

				// filter none
				for (i = 0; i < nframes; ++i) {
					out[0][i] = X_full[0][i+window_size_2+nframes_2];
					out[1][i] = 0.0;	
				}
			}
		}

	} else {
		// filter all
		for (unsigned int s2f = 1; s2f <= n_sources; ++s2f)
		{
			if (ecounter[s2f-1] > 0) 	// if at least one valid DOA was found of the target source, apply beamforming
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

					max_max = 0.0f;

					if (max_val12 > max_max)
						max_max = max_val12;
					if (max_val23 > max_max)
						max_max = max_val23;
					if (max_val31 > max_max)
						max_max = max_val31;

					if (max_max == 0.0f) {
						if(state == 1) exit(1);
					}
					else {
						if(state == 0) state = 1;

						// perform overlap-add:		
						for (i = 0; i < nframes; ++i) {
							write_buffer[i] = X_full[0][i+window_size_2+nframes_2];;
							for (k = 1; k < n_in_channels; ++k)
							{
								write_buffer[i] += X_late[k][i+window_size_2+nframes_2] + X_early[k][i+nframes_2];
							}
							write_buffer[i] /= 3.0;
						}

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
			} else {	// don't start writing if no audio is being played
				max_max = 0.0f;

				if (max_val12 > max_max)
					max_max = max_val12;
				if (max_val23 > max_max)
					max_max = max_val23;
				if (max_val31 > max_max)
					max_max = max_val31;

				if (max_max == 0.0f) {
					if(state == 1) exit(1);
				}
				else {
					if(state == 0) state = 1;

					for (i = 0; i < nframes; ++i) {
						write_buffer[i] = X_full[0][i+window_size_2+nframes_2];
					}
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

	unsigned int i, j;

	if(argc != 4 && argc != 3){		
		printf ("Usage:\ngcc_beamformer d N k\nd: Microphone separation (in meters).\nN: Maximum number of sources.\nk: Which source to filter (optional parameter, the default option is filtering all sources).\n");
		exit(1);
	}

	mic_separation = atof(argv[1]);
	printf("\nMicrophone separation: %1.4f meters.", mic_separation);	
	n_sources = atoi(argv[2]);
	printf("\nMaximum number of sources: %d.\n\n", n_sources);

	if (argc == 4)
	{
		source2filter = atoi(argv[3]);

		if (source2filter > (int) n_sources || source2filter < -1) {
			printf("The source to be filtered is not between 1 and the maximum number of sources provided. Quitting ...\n\n");
			exit(1);
		}

		if (source2filter == 0)
			printf("\nAll sources are going to be filtered.\n\n");
		else
			printf("\nSource number %d is going to be filtered.\n\n", source2filter);
	}

	system("mkdir -p output");
	char audio_file_path[100];
	sprintf(audio_file_path, "./output/bf_audio_%d_of_%d.wav", source2filter, n_sources);

	char text_file_path[100];
	if (n_sources == 1)
		sprintf(text_file_path, "./output/track_%d_source.txt", n_sources);
	else		
		sprintf(text_file_path, "./output/track_%d_sources.txt", n_sources);
	outputFile.open(text_file_path);

	char text_kalman_path[100];
	if (n_sources == 1)
		sprintf(text_kalman_path, "./output/track_%d_sourceKalman.txt", n_sources);
	else		
		sprintf(text_kalman_path, "./output/track_%d_sourcesKalman.txt", n_sources);
	outputKalman.open(text_kalman_path);

	dt_max = mic_separation/c;
	N_max = dt_max*sample_rate;

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
	kmin = (int) (f_min/sample_rate*window_size_2);
	kmax = (int) (f_max/sample_rate*window_size_2);

	fRes = 2.0*M_PI/window_size;

	hist_length = MEMORY_FACTOR*n_sources;

	// initialization of internal buffers
	// - overlap-add buffers
	X_late		= (jack_default_audio_sample_t **) calloc(n_in_channels, sizeof(jack_default_audio_sample_t*));
	X_early		= (jack_default_audio_sample_t **) calloc(n_in_channels, sizeof(jack_default_audio_sample_t*));
	X_full		= (jack_default_audio_sample_t **) calloc(n_in_channels, sizeof(jack_default_audio_sample_t*));
	X_gcc		= (std::complex<double> **) calloc(n_in_channels, sizeof(std::complex<double>*));
	Aux_gcc		= (std::complex<double> *) calloc(window_size_2, sizeof(std::complex<double>));
	DOA_hist	= (double *) calloc(hist_length, sizeof(double));
	DOA_class	= (unsigned int *) calloc(hist_length, sizeof(unsigned int));
	DOA_kmean	= (double *) calloc(n_sources, sizeof(double));
	DOA_mean	= (double *) calloc(n_sources, sizeof(double));
	DOA_stdev	= (double *) calloc(n_sources, sizeof(double));
	DOA_valid	= (double *) calloc(n_sources, sizeof(double));
	kalmanState = (double **) calloc(4, sizeof(double*));
	covMatrix   = (double **) calloc(4*n_sources, sizeof(double*));

	counter		= (unsigned int *) calloc(n_sources, sizeof(unsigned int));
	dcounter	= (int *) calloc(n_sources, sizeof(int));
	ecounter	= (int *) calloc(n_sources, sizeof(int));
	write_buffer = (float *) calloc(nframes, sizeof(float));

	std::default_random_engine generator;
  	std::uniform_real_distribution<double> distribution(-180.0,180.0);

	for (i = 0; i < n_sources; ++i)
	{
		DOA_kmean[i] = 360.0/n_sources*i; //+ distribution(generator); 360.0/2.0/n_sources;	// "optimal" initialization for the k-means algorithm
		//DOA_kmean[i] = distribution(generator);	// random initialization for the k-means algorithm
		if (DOA_kmean[i] > 180.0) {
			DOA_kmean[i] -= 360.0;
		}
		//DOA_kmean[i] = distribution(generator);	// random initialization for the k-means algorithm
		DOA_stdev[i] = 360.0/2.0/n_sources;
	}

	for (i = 0; i < 4; ++i) {
		kalmanState[i] = (double *) calloc(n_sources, sizeof(double));
		for (j = 0; j < n_sources; ++j) {
			covMatrix[j*4+i] = (double *) calloc(4, sizeof(double));
		}
	}

	// initialize Kalman state:
	double initialState[2];
	for (i = 0; i < n_sources; ++i) {
		angle2state(DOA_kmean[i], initialState);

		kalmanState[0][i] = initialState[0];
		kalmanState[1][i] = initialState[1];
		
		/*
		covMatrix[i*4+0][0] = 0.000002493289818;	covMatrix[i*4+0][2] = 0.000116781724500;
		covMatrix[i*4+1][1] = 0.000002493289818;	covMatrix[i*4+1][3] = 0.000116781724500;
		covMatrix[i*4+2][0] = 0.000116781724500;	covMatrix[i*4+2][2] = 0.005469870000000;
		covMatrix[i*4+3][1] = 0.000116781724500;	covMatrix[i*4+3][3] = 0.005469870000000;
		*/
	}


    for (i = 0; i < n_in_channels; ++i) {
        X_late[i]	= (jack_default_audio_sample_t *) calloc(window_size, sizeof(jack_default_audio_sample_t));
		X_early[i]	= (jack_default_audio_sample_t *) calloc(window_size, sizeof(jack_default_audio_sample_t));
		X_full[i]	= (jack_default_audio_sample_t *) calloc(window_size + window_size_2, sizeof(jack_default_audio_sample_t));
		X_gcc[i]	= (std::complex<double> *) calloc(window_size_2, sizeof(std::complex<double>));
    }	

	// - FFTW3 buffers
	i_fft_4N = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size);
	i_time_4N = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size);
	o_fft_4N = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size);
	o_time_4N = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size);
	
	i_forward_4N = fftw_plan_dft_1d(window_size, reinterpret_cast<fftw_complex*>(i_time_4N), reinterpret_cast<fftw_complex*>(i_fft_4N), FFTW_FORWARD, FFTW_MEASURE);
	o_inverse_4N = fftw_plan_dft_1d(window_size, reinterpret_cast<fftw_complex*>(o_fft_4N), reinterpret_cast<fftw_complex*>(o_time_4N), FFTW_BACKWARD, FFTW_MEASURE);

	i_fft_2N = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size_2);
	i_time_2N = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size_2);
	o_fft_2N = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size_2);
	o_time_2N = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * window_size_2);

	i_forward_2N = fftw_plan_dft_1d(window_size_2, reinterpret_cast<fftw_complex*>(i_time_2N), reinterpret_cast<fftw_complex*>(i_fft_2N), FFTW_FORWARD, FFTW_MEASURE);
	o_inverse_2N = fftw_plan_dft_1d(window_size_2, reinterpret_cast<fftw_complex*>(o_fft_2N), reinterpret_cast<fftw_complex*>(o_time_2N), FFTW_BACKWARD, FFTW_MEASURE);
	
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

	if (source2filter == 0) {
		for (i = 0; i < n_sources; ++i)
		{
			printf("Trying to open audio File: ./output/bf_audio_%d_of_%d.wav\n", i+1, n_sources);
			audio_info.samplerate = sample_rate;
			audio_info.channels = 1;
			audio_info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_32;

			sprintf(audio_file_path, "./output/bf_audio_%d_of_%d.wav", i+1, n_sources);
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
	} else {
		if (source2filter != -1)
		{
			printf("Trying to open audio File: %s\n",audio_file_path);
			audio_info.samplerate = sample_rate;
			audio_info.channels = 1;
			audio_info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_32;

			saudio_file = sf_open (audio_file_path,SFM_WRITE,&audio_info);
			if(saudio_file == NULL){
				printf("%s\n",sf_strerror(NULL));
				exit(1);
			}else{
				printf("Audio file info:\n");
				printf("\tSample Rate: %d\n",audio_info.samplerate);
				printf("\tChannels: %d\n",audio_info.channels);
			}
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
	//sf_close(audio_file);
	outputFile.close();
	outputKalman.close();
	exit (0);
}