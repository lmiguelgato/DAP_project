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
#include <string>
#include <iomanip>
using namespace std;

ofstream outputFile;					// save results for data analysis
ofstream outputKalman;					// save results for data analysis

#include "internals.h"

#define VERBOSE false					// display additional information

int main (int argc, char *argv[]) {

	if(argc != 2 && argc != 3){		
		printf ("Usage:\ngcc_beamformer_offline <audio file path> <output path>\n");
		exit(1);
	}

	string line;
	ifstream settings_file ("array_settings.txt");
	if (settings_file.is_open())
	{
		getline (settings_file,line);
		getline (settings_file,line);
		microphone_positions[0] = stof(line);
		getline (settings_file,line);
		microphone_positions[1] = stof(line);

		getline (settings_file,line);
		getline (settings_file,line);
		microphone_positions[2] = stof(line);
		getline (settings_file,line);
		microphone_positions[3] = stof(line);

		getline (settings_file,line);
		getline (settings_file,line);
		microphone_positions[4] = stof(line);
		getline (settings_file,line);
		microphone_positions[5] = stof(line);

		settings_file.close();

		measure_array_geometry(microphone_positions, sides, angles);
	}
	else {
		cout << "Unable to open settings file.\n";
		exit(1);
	}

	char* full_path;
	const char* str1 = "gcc_ssl.txt";
	const char* str2 = "gcc_sst.txt";

	if(argc == 2){
		printf ("No output path specified. The current directory will be used ...\n");
		sprintf(output_file_path, "%s", "./");

		const char* str3 = "./";
		full_path = (char*) malloc(strlen(str3)+strlen(str1)+1); 
		strcpy(full_path, str3);
		strcat(full_path, str1);
		outputFile.open(full_path);

		full_path = (char*) malloc(strlen(str3)+strlen(str2)+1); 
		strcpy(full_path, str3);
		strcat(full_path, str2);
		outputKalman.open(full_path);
	} else {		
		sprintf(output_file_path, "%s", argv[2]);
		
		full_path = (char*) malloc(strlen(argv[2])+strlen(str1)+1); 
		strcpy(full_path, argv[2]);
		strcat(full_path, str1);
		outputFile.open(full_path);

		full_path = (char*) malloc(strlen(argv[2])+strlen(str2)+1); 
		strcpy(full_path, argv[2]);
		strcat(full_path, str2);
		outputKalman.open(full_path);
	}

	sprintf(audio_file_path, "%s", argv[1]);
	FILE* audiodata = fopen(audio_file_path, "r");
	if (audiodata ==  NULL) {
		printf("\nERROR: Could not open audio data file.\n\n");
		exit(2);
	}
    fseek(audiodata, 22, SEEK_SET);

	unsigned char buffer2[2];

	if(fread(&buffer2, sizeof(buffer2), 1, audiodata)) {
		n_in_channels = buffer2[0] | (buffer2[1] << 8);
		printf("%d audio channels detected ...\n", n_in_channels);
		if (n_in_channels != 3) {
			printf("This demo only supports 3 input channels. Aborting ...\n");
			return 0;
		}
	} else {
		printf("Error reading audio file. Aborting ...\n");
		exit(2);
	}

	unsigned char buffer4[4];

	if (fread(&buffer4, sizeof(buffer4), 1, audiodata)) {
		sample_rate = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) | (buffer4[3] << 24);
		printf("%d Hz sampling rate ...\n", (unsigned int) sample_rate);
	} else {
		printf("Error reading audio file. Aborting ...\n");
		exit(2);
	}
	
    int16_t recvsample_s;
	nframes = 1024;
	init();		// Initialize all internal registers, variables and memory allocation

	float * in = (float *)malloc(nframes*n_in_channels*sizeof(float));

	fseek(audiodata, 44, SEEK_SET);

	while (!feof(audiodata)) {
		
		for (unsigned int i = 0; i < (unsigned int) nframes; i++) {
			for (unsigned int j = 0; j < n_in_channels; j++) {
				if(fread(&recvsample_s, sizeof(int16_t), 1, audiodata)) {
					in[i*n_in_channels + j] = recvsample_s / 32767.0f;
				}
			}
		}
		process_audio(in);
	}
	
	outputFile.close();
	outputKalman.close();
	return 0;
}

void process_audio (float *in){

	unsigned int i, j;

	for (j = 0; j < n_in_channels; ++j) {

		// shift-register (useful for overlap-add method):
		for (i = 0; i < (unsigned int) nframes; ++i)
		{
			X_full[j][i] 						= X_full[j][nframes + i];
			X_full[j][nframes + i] 				= X_full[j][window_size_2 + i];
			X_full[j][window_size_2 + i] 		= X_full[j][window_size-nframes + i];
			X_full[j][window_size-nframes + i] 	= X_full[j][window_size + i];
			X_full[j][window_size + i] 			= X_full[j][window_size+nframes + i];
			X_full[j][window_size+nframes + i] 	= in[i*n_in_channels + j];	
		}

		// cross-correlation in four steps:
		// 1- zero padding:
		for (i = 0; i < (unsigned int) nframes; ++i) {
			i_time_2N[i] = in[i*n_in_channels + j];
			//printf("%6.12f, ", real(i_time_2N[i]));
		}

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
	double theta[6] = {asin(  unwrap(  max(X_gcc[1], window_size_2, N_max[0], &max_val12), nframes, N_max[0]  )/N_max[0]  )*RAD2DEG,
					   asin(  unwrap(  max(X_gcc[2], window_size_2, N_max[1], &max_val23), nframes, N_max[1]  )/N_max[1]  )*RAD2DEG,
					   asin(  unwrap(  max(X_gcc[0], window_size_2, N_max[2], &max_val31), nframes, N_max[2]  )/N_max[2]  )*RAD2DEG,
					   0.0,
					   0.0,
					   0.0};

	angleTranslation(theta, angles);	// use a coherent reference to measure DOA

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

if (VERBOSE)
{
					printf("kalman initialization[%d]: %1.1f\n", i, DOA_kmean[i]);
}
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
if (VERBOSE)
{
						printf("DOA[%d] = %1.1f\n", i, DOA_valid[i]);
}
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
}

void measure_array_geometry (float* mic_pos, float* side, float* angle) {
	float dx, dy;
	float d12, d23, d31;

	dx = mic_pos[0] - mic_pos[2];
	dy = mic_pos[1] - mic_pos[3];
	side[0] = sqrt( dx*dx + dy*dy );
	d12 = side[0];

	dx = mic_pos[2] - mic_pos[4];
	dy = mic_pos[3] - mic_pos[5];
	side[1] = sqrt( dx*dx + dy*dy );
	d23 = side[1];

	dx = mic_pos[4] - mic_pos[0];
	dy = mic_pos[5] - mic_pos[1];
	side[2] = sqrt( dx*dx + dy*dy );
	d31 = side[2];

	angle[0] = acos((d12*d12 + d31*d31 - d23*d23)/2/d12/d31);
	angle[1] = acos((d12*d12 + d23*d23 - d31*d31)/2/d12/d23);
	angle[2] = M_PI - angle[0] - angle[1];
}

void init(void) {

	N_max[0] = (int) (sides[0]/c*sample_rate);
	N_max[1] = (int) (sides[1]/c*sample_rate);
	N_max[2] = (int) (sides[2]/c*sample_rate);

	// obtain here the delay from user and store it in 'delay'
	nframes_2   = nframes/2;
	window_size = 4*nframes;
	window_size_2 = 2*nframes;
	kmin = (int) (f_min/sample_rate*window_size_2);
	kmax = (int) (f_max/sample_rate*window_size_2);

	fRes = 2.0*M_PI/window_size;

	hist_length = MEMORY_FACTOR*n_sources;

	// initialization of internal buffers
	// - overlap-add buffers
	X_late		= (float **) calloc(n_in_channels, sizeof(float*));
	X_early		= (float **) calloc(n_in_channels, sizeof(float*));
	X_full		= (float **) calloc(n_in_channels, sizeof(float*));
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

	for (unsigned int i = 0; i < n_sources; ++i)
	{
		DOA_kmean[i] = 360.0/n_sources*i; //+ distribution(generator); 360.0/2.0/n_sources;	// "optimal" initialization for the k-means algorithm
		//DOA_kmean[i] = distribution(generator);	// random initialization for the k-means algorithm
		if (DOA_kmean[i] > 180.0) {
			DOA_kmean[i] -= 360.0;
		}
		//DOA_kmean[i] = distribution(generator);	// random initialization for the k-means algorithm
		DOA_stdev[i] = 360.0/2.0/n_sources;
	}

	for (unsigned int i = 0; i < 4; ++i) {
		kalmanState[i] = (double *) calloc(n_sources, sizeof(double));
		for (unsigned int j = 0; j < n_sources; ++j) {
			covMatrix[j*4+i] = (double *) calloc(4, sizeof(double));
		}
	}

	// initialize Kalman state:
	for (unsigned int i = 0; i < n_sources; ++i) {
		angle2state(DOA_kmean[i], initialState);

		kalmanState[0][i] = initialState[0];
		kalmanState[1][i] = initialState[1];
	}


    for (unsigned int i = 0; i < n_in_channels; ++i) {
        X_late[i]	= (float *) calloc(window_size, sizeof(float));
		X_early[i]	= (float *) calloc(window_size, sizeof(float));
		X_full[i]	= (float *) calloc(window_size + window_size_2, sizeof(float));
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

	// - hann window
	hann = (float *) calloc(window_size, sizeof(float)); 
	for(unsigned int i = 0; i < window_size; ++i) {
		hann[i] = 0.5 - 0.5*cos(2.0*M_PI* ((double) i/(window_size-1)));
	}
}