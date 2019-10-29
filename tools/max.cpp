#include "max.h"

int max (std::complex<double> *signal, int signal_length, int Nmax, double *max_value) {
	int max_index = -1;
	double tmp2 = -1.0;

	max_value[0] = -1.0;
	for (int i = 0; i < Nmax; ++i) {
		tmp2 = abs(signal[i]);
		if(tmp2 > max_value[0]) {
			max_index = i;
			max_value[0] = tmp2;
		}
	}
	for (int i = signal_length-Nmax; i < signal_length; ++i) {
		tmp2 = abs(signal[i]);
		if(tmp2 > max_value[0]) {
			max_index = i;
			max_value[0] = tmp2;
		}
	}

	return max_index;
}