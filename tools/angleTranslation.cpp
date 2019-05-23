#include "angleTranslation.h"

void angleTranslation (double *theta) {

	int i;

	for (i = 0; i < 3; ++i) {
		if (theta[i] > 0.0) {
			theta[i+3] = 180.0 - theta[i];
		} else {
			theta[i+3] = -180.0 - theta[i];
		}
	}

	theta[1] -= 120.0;
	theta[2] += 120.0;
	theta[4] -= 120.0;
	theta[5] += 120.0;

	for (i = 0; i < 6; ++i) {
		if (theta[i] > 180.0) {
			theta[i] -= 360.0;
		}
		if (theta[i] < -180.0) {
			theta[i] += 360.0;
		}
	}

	return;
}