#ifndef XCORR
#define XCORR
// FFTW: library for computing the discrete Fourier transform of arbitrary input size, 
//       and of both real and complex data.
#include <complex.h> 	// needs to be included before fftw3.h for compatibility
#include <fftw3.h>

void xcorr (std::complex<double> **, int, fftw_plan, std::complex<double> *, std::complex<double> *);

#endif
