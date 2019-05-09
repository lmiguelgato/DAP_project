#include "max.h"
#define DELTA 0.0000000001f

int max (std::complex<double> *signal, int signal_length) {
	int max_index = -1;
	double tmp1 = -1.0;
	double tmp2 = -1.0;

	tmp1 = -1.0;
	for (int i = 0; i < signal_length; ++i) {
		tmp2 = abs(signal[i]);
		if(tmp2 > tmp1+DELTA) {
			max_index = i;
			tmp1 = tmp2;
		}
	}

	return max_index;
}