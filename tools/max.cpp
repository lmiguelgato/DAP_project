#include "max.h"
//#define DELTA 0.0000000001f

int max (std::complex<double> *signal, int signal_length, double *max_value) {
	int max_index = -1;
	double tmp2 = -1.0;

	max_value[0] = -1.0;
	for (int i = 0; i < signal_length; ++i) {
		tmp2 = abs(signal[i]);
		if(tmp2 > max_value[0]) {
			max_index = i;
			max_value[0] = tmp2;
		}
	}

	return max_index;
}