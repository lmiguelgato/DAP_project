#include "xcorr.h"

void xcorr (std::complex<double> ** X, int X_length, fftw_plan ifft, std::complex<double> *time, std::complex<double> *freq) {
	int i;

	for (i = 0; i < X_length; ++i)
	{
		freq[i] = X[0][i] * conj(X[1][i]);
	}
	fftw_execute(ifft);
	for (i = 0; i < 2*nframes; ++i)
	{
		Aux_gcc[i] = o_time_2N[i];
	}
}