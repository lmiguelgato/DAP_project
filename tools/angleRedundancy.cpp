#include "angleRedundancy.h"

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