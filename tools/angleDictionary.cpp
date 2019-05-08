#include "angleDictionary.h"

void angleDictionary (double *theta) {

	double temp = theta[1];

	theta[1] = temp - 120.0;
	theta[4] = 60.0 - temp;

	if (theta[1] > 180.0)
		theta[1] -= 360.0;
	if (theta[1] < -180.0)
		theta[1] += 360.0;
	if (theta[4] > 180.0)
		theta[4] -= 360.0;
	if (theta[4] < -180.0)
		theta[4] += 360.0;

	temp = theta[2];

	theta[2] = temp + 120.0;
	theta[5] = temp - 60.0;

	if (theta[2] > 180.0)
		theta[2] -= 360.0;
	if (theta[2] < -180.0)
		theta[2] += 360.0;
	if (theta[5] > 180.0)
		theta[5] -= 360.0;
	if (theta[5] < -180.0)
		theta[5] += 360.0;

	theta[3] = 180.0 - theta[0];

	if (theta[2] > 180.0)
		theta[2] -= 360.0;
	if (theta[2] < -180.0)
		theta[2] += 360.0;
	if (theta[5] > 180.0)
		theta[5] -= 360.0;
	if (theta[5] < -180.0)
		theta[5] += 360.0;

/*
	if (theta[1] > -60.0 && theta[1] <= 90.0)
	{
		theta[1] = temp - 120.0;
		theta[4] = 60.0 - temp;
	} else {
		theta[1] = 120.0 - temp;
		theta[4] = 60.0 - temp; 
	}

	temp = theta[2];

	if (theta[2] > -90.0 && theta[2] <= 60.0)
	{
		theta[2] = 120.0 + temp;
		theta[5] = -60.0 - temp;
	} else {
		theta[2] = -240 + temp;
		theta[5] = -60.0 - temp;
	}
	*/

	return;
}