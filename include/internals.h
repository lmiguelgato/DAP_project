#ifndef __SYNC_H__
#define __SYNC_H__

// Useful functions developed by myself:
#include "max.h"
#include "unwrap.h"
#include "phat.h"
#include "angleManipulation.h"
#include "kmeans.h"
#include "kalman.h"

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
#define GCC_TH 140.0f					// correlation threshold (to avoid false alarms)
#define REDUNDANCY_TH 20.0f				// redundancy threshold (for DOA estimation)
#define DYNAMIC_GCC_TH 0				// enable a dynamic GCC threshold (0: disabled, 1: mean peak values, 2: max peak values)
#define MOVING_AVERAGE 1				// enable a moving average on kmeans centroids (0: disabled, 1: finite memory, 2: infinite memory)
#define MOVING_FACTOR 1					// allow variations in DOA if the sources are moving (how many times the standard deviation)
#define MEMORY_FACTOR 5				    // memory of the k-means algorithm
#define n_sources 4

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
double sample_rate;			            // sample rate [Hz]
int nframes = 1024;				        // default number of frames per jack buffer
unsigned int window_size, window_size_2, nframes_2;
double c = 343.364;						// default sound speed [m/s]
int N_max[3];				        // maximum delay between microphones [samples]
double doa;								// direction of arrival
double mean_doa = 0.0;					// average direction of arrival
double std_doa = 90.0;					// standard deviation of the direction of arrival
double std_cum = 90.0;					// standard deviation of the direction of arrival (cummulative)
double gcc_th = GCC_TH;					// default GCC threshold
double fRes; 							// frequency resolution

double f_min = 1000.0;					// minimum frequency of the desired signal [Hz]
double f_max = 4000.0;					// maximum frequency of the desired signal [Hz]
int kmin, kmax;							// discrete minimum and maximum frequencies of the desired signal

unsigned int n_in_channels;			    // number of input channels

int state = 0;							// beamformer state

float *hann;		// array to store the Hann window

// overlap-add registers:
float **X_full;	// store the 6 last buffers of 'nframes' samples
float **X_late;	// store the 4 latest buffers from X_full
float **X_early;	// store the 4 earliest buffers from X_full

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
double initialState[2];

float microphone_positions[6];
float sides[3];
float angles[3];

char audio_file_path[1000];
char output_file_path[1000];
char text_file_path[20];
char text_kalman_path[20];

void init(void);    // Initialize all internal registers, variables and memory allocation
void process_audio (float*);
void measure_array_geometry (float*, float*, float*);

#endif