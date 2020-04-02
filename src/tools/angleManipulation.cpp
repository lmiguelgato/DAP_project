#include "angleManipulation.h"

double state2angle (double* state) {

	/*
		Convert a state (Cartesian coordinates) to the corresponding angle

		input: 	state (Cartesian coordinates)
		output: (the corresponding angle)
	*/

	if (abs(state[1]) <= 0.001) {
		if (state[0] >= 0.0)
			return 90.0;
		else
			return -90.0;
	} else {
		if (state[1] >= 0) {
			return atan(state[0]/state[1])*RAD2DEG;
		} else {	
			if (state[0] > 0) {
				return atan(state[0]/state[1])*RAD2DEG + 180.0;
			} else {
				if (state[0] < 0)
					return atan(state[0]/state[1])*RAD2DEG - 180.0;
				else
					return 180.0;
			}
		}
	}
	
}

void angle2state (double angle, double* state) {

	/*
		Convert an angle to its corresponding state (Cartesian coordinates)

		input: 	an angle (in radians)
		output: state (the corresponding Cartesian coordinates)
	*/

	state[0] = sin(angle*DEG2RAD);
	state[1] = cos(angle*DEG2RAD);
	
}

double angleRedundancy (double* theta, double* thetaRedundant, double Ethresh) {
	int i, j, k;
	double theta12, theta23, theta31;
	double Epqr, minAngleDiff;
	double minEpqr = 100000;
	double doa;
	double xy_1[2];
	double xy_2[2];
	double xy_diff[2];
	bool found = false;

	for (i = 0; i < 2; ++i)	{
		theta12 = theta[i*3];
		for (j = 0; j < 2; ++j)	{
			theta23 = theta[1+j*3];
			for (k = 0; k < 2; ++k)	{
				theta31 = theta[2+k*3];

				angle2state (theta12, xy_1);
				angle2state (theta23, xy_2);

				xy_diff[0] = xy_1[0]-xy_2[0];
				xy_diff[1] = xy_1[1]-xy_2[1];

				Epqr = xy_diff[0]*xy_diff[0] + xy_diff[1]*xy_diff[1];
				if (Epqr < minEpqr)
				{
					minEpqr = Epqr;
					minAngleDiff = abs(theta12-theta23);

					if (minAngleDiff <= Ethresh) {
						found = true;
						doa = (theta12+theta23)/2;
					}
					if (minAngleDiff >= (360.0-Ethresh)) {
						found = true;
						doa = (abs(theta12)+abs(theta23))/2;
					}
					
				}

				angle2state (theta31, xy_1);

				xy_diff[0] = xy_1[0]-xy_2[0];
				xy_diff[1] = xy_1[1]-xy_2[1];

				Epqr = xy_diff[0]*xy_diff[0] + xy_diff[1]*xy_diff[1];
				if (Epqr < minEpqr)
				{
					minEpqr = Epqr;
					minAngleDiff = abs(theta31-theta23);
					if (minAngleDiff <= Ethresh) {
						found = true;
						doa = (theta31+theta23)/2;
					}
					if (minAngleDiff >= (360.0-Ethresh)) {
						found = true;
						doa = (abs(theta31)+abs(theta23))/2;
					}
				}

				angle2state (theta12, xy_2);

				xy_diff[0] = xy_1[0]-xy_2[0];
				xy_diff[1] = xy_1[1]-xy_2[1];

				Epqr = xy_diff[0]*xy_diff[0] + xy_diff[1]*xy_diff[1];
				if (Epqr < minEpqr)
				{
					minEpqr = Epqr;
					minAngleDiff = abs(theta31-theta12);
					if (minAngleDiff <= Ethresh) {
						found = true;
						doa = (theta12+theta31)/2;
					}
					if (minAngleDiff >= (360.0-Ethresh)) {
						found = true;
						doa = (abs(theta12)+abs(theta31))/2;
					}
				}
			}
		}
	}

	if (!found) 	doa = 181.0;

	return doa;

}

void angleTranslation (double *theta, float *alpha) {

	int i;

	for (i = 0; i < 3; ++i) {
		if (theta[i] >= 0.0) {
			theta[i+3] = 180.0 - theta[i];
		} else {
			theta[i+3] = -180.0 - theta[i];
		}
	}

	theta[1] -= 180.0 - alpha[1]*RAD2DEG;
	theta[2] += 180.0 - alpha[0]*RAD2DEG;
	theta[4] -= 180.0 - alpha[1]*RAD2DEG;
	theta[5] += 180.0 - alpha[0]*RAD2DEG;

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

void angleTranslation (double *theta) {

	int i;

	for (i = 0; i < 3; ++i) {
		if (theta[i] >= 0.0) {
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