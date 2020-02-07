#include "phat.h"

void phat (std::complex<double>* Z, std::complex<double>* X, std::complex<double>* Y, int X_length, int kmin, int kmax, int type) {

	int i;
	std::complex<double> temp;

	switch ( type ) {
		case 2 :	// GCC (frequency restrained)
			for (i = 0; i < kmin; ++i)							Z[i] = 0.0;
			for (i = kmax; i < X_length-kmax; ++i)				Z[i] = 0.0;
			for (i = X_length-kmin; i < X_length; ++i)			Z[i] = 0.0;	
			for (i = kmin; i < kmax; ++i)						Z[i] = X[i] * conj(Y[i]);
			for (i = X_length-kmax; i < X_length-kmin; ++i)		Z[i] = X[i] * conj(Y[i]);
			break;
		case 3 :	// GCC-PHAT
			for (i = 0; i < X_length; ++i) {
				temp = X[i] * conj(Y[i]);
				(abs(temp) == 0.0) ? (Z[i] = 0.0) : (Z[i] = temp/abs(temp));
			}
			break;
		case 4 :	// GCC-PHAT (frequency restrained)
			for (i = 0; i < kmin; ++i)							Z[i] = 0.0;
			for (i = kmax; i < X_length-kmax; ++i)				Z[i] = 0.0;
			for (i = X_length-kmin; i < X_length; ++i) 			Z[i] = 0.0;		
			for (i = kmin; i < kmax; ++i) {
				temp = X[i] * conj(Y[i]);
				(abs(temp) == 0.0) ? (Z[i] = 0.0) : (Z[i] = temp/abs(temp));
			}
			for (i = X_length-kmax; i < X_length-kmin; ++i)	{
				temp = X[i] * conj(Y[i]);
				(abs(temp) == 0.0) ? (Z[i] = 0.0) : (Z[i] = temp/abs(temp));
			}
			break;
		default : 	// GCC
      		for (i = 0; i < X_length; ++i)	Z[i] = X[i] * conj(Y[i]);
			break;
	}
}