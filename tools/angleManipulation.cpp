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
	double Epqr;
	double minEpqr = 1000;
	double doa;
	bool found = false;

	for (i = 0; i < 2; ++i)	{
		theta12 = theta[i];
		for (j = 0; j < 2; ++j)	{
			theta23 = theta[2+i];
			for (k = 0; k < 2; ++k)	{
				theta31 = theta[4+i];
				Epqr = (abs(theta12-theta23) + abs(theta23-theta31) + abs(theta31-theta12))/3.0;
				if (Epqr < minEpqr)
				{
					minEpqr = Epqr;
					thetaRedundant[0] = theta12;
					thetaRedundant[1] = theta23;
					thetaRedundant[2] = theta31;
				}
			}
		}
	}



	minEpqr = Ethresh;
	Epqr = abs(thetaRedundant[0]-thetaRedundant[1]);
	if (Epqr < minEpqr) {
		doa = (thetaRedundant[0]+thetaRedundant[1])/2.0;
		minEpqr = Epqr;
		found = true;
	}

	Epqr = abs(thetaRedundant[1]-thetaRedundant[2]);
	if (Epqr < minEpqr) {
		doa = (thetaRedundant[1]+thetaRedundant[2])/2.0;
		minEpqr = Epqr;
		found = true;
	}

	Epqr = abs(thetaRedundant[2]-thetaRedundant[0]);
	if (Epqr < minEpqr) {
		doa = (thetaRedundant[2]+thetaRedundant[0])/2.0;
		found = true;
	}

	if (!found) 	doa = 181.0;

	return doa;

}

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