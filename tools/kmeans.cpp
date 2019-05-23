#include "kmeans.h"
#include <stdio.h>
#include <math.h>
#define DEG2RAD 0.017453292519943		// useful to convert from degrees to radians
#define RAD2DEG 57.295779513082323f		// useful to convert from radians to degrees

void kmeans (double* DOA_hist, double* DOA_mean, int* class_count, int n_sources, int hist_length) {
	int i, j;
	double min_dist, dist;
	double x, y;

	int *DOA_class;
	DOA_class = (int *) calloc(hist_length, sizeof(int));

	for (i = 0; i < n_sources; ++i)
	{
		class_count[i] = 0;
	}

	for (j = 0; j < hist_length; ++j)
	{
		min_dist = 10; 	// a number grater than 4
		for (i = 0; i < n_sources; ++i)
		{			
			//dist = abs(DOA_hist[j] - DOA_mean[i]);
			dist = pow(cos(DOA_hist[j]*DEG2RAD) - cos(DOA_mean[i]*DEG2RAD), 2) + pow(sin(DOA_hist[j]*DEG2RAD) - sin(DOA_mean[i]*DEG2RAD), 2);
			if (dist < min_dist) {
				min_dist = dist;
				DOA_class[j] = i;
			}
		}
		++class_count[ DOA_class[j] ];
	}

	for (i = 0; i < n_sources; ++i)
	{			
		//DOA_mean[i] = 0.0;
		//printf("i = %d\n", i);
		for (j = 0; j < hist_length; ++j)
		{
			if 	(DOA_class[j] == i) {
				DOA_mean[i] += DOA_hist[j];
				//printf("DOA_hist[%d] = %1.5f\n", j, DOA_hist[j]);
			}
		}
		if (class_count[i] != 0)
			DOA_mean[i] /= (class_count[i]+1);

		/*if (class_count[i] != 0) {
			DOA_mean[i] /= (class_count[i]);
		} else {
			DOA_mean[i] = 360.0/n_sources*i + 360.0/2.0/n_sources;
			if (DOA_mean[i] > 180.0) {
				DOA_mean[i] -= 360.0;
			}
		}*/
	}

	/*for (i = 0; i < n_sources; ++i)
	{	
		x = 0.0;
		y = 0.0;
		for (j = 0; j < hist_length; ++j)
		{
			if 	(DOA_class[j] == i) {
				x += cos(DOA_hist[j]*DEG2RAD);
				y += sin(DOA_hist[j]*DEG2RAD);
			}
		}
		if (class_count[i] != 0) {
			if (x < 0) {
				DOA_mean[i] = atan(y/x)*RAD2DEG - 180.0;
			} else {
				DOA_mean[i] = atan(y/x)*RAD2DEG;
			}
		}
	}*/
}