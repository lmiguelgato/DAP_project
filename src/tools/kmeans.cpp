#include "kmeans.h"
#include <math.h>
#define DEG2RAD 0.017453292519943		// useful to convert from degrees to radians
#define RAD2DEG 57.295779513082323f		// useful to convert from radians to degrees

void kmeans (double* DOA_hist, double* DOA_mean, unsigned int* class_count, unsigned int n_sources, unsigned int hist_length) {
	unsigned int i, j;
	double min_dist, dist;

	unsigned int *DOA_class;
	DOA_class = (unsigned int *) calloc(hist_length, sizeof(unsigned int));

	for (i = 0; i < n_sources; ++i)
	{
		class_count[i] = 0;
	}

	for (j = 0; j < hist_length; ++j)
	{
		min_dist = 10; 	// a number grater than 4
		for (i = 0; i < n_sources; ++i)
		{
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
		for (j = 0; j < hist_length; ++j)
		{
			if 	(DOA_class[j] == i) {
				DOA_mean[i] += DOA_hist[j];
			}
		}
		if (class_count[i] != 0)
			DOA_mean[i] /= (class_count[i]+1);
	}
}

void kmeans (double* DOA_hist, unsigned int *DOA_class, double* DOA_mean, unsigned int* class_count, unsigned int n_sources, unsigned int hist_length) {
	unsigned int i, j;
	double min_dist, dist;

	for (i = 0; i < n_sources; ++i)
	{
		class_count[i] = 0;
	}

	for (j = 0; j < hist_length; ++j)
	{
		min_dist = 10; 	// a number grater than 4
		for (i = 0; i < n_sources; ++i)
		{
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
		for (j = 0; j < hist_length; ++j)
		{
			if 	(DOA_class[j] == i) {
				DOA_mean[i] += DOA_hist[j];
			}
		}
		if (class_count[i] != 0)
			DOA_mean[i] /= (class_count[i]+1);
	}
}

void kmeans (double* DOA_hist, unsigned int *DOA_class, double* DOA_mean, unsigned int n_sources, unsigned int hist_length) {
	unsigned int i, j;
	double min_dist, dist;

	for (j = 0; j < hist_length; ++j)
	{
		min_dist = 10; 	// a number grater than 4
		for (i = 0; i < n_sources; ++i)
		{
			dist = pow(cos(DOA_hist[j]*DEG2RAD) - cos(DOA_mean[i]*DEG2RAD), 2) + pow(sin(DOA_hist[j]*DEG2RAD) - sin(DOA_mean[i]*DEG2RAD), 2);
			if (dist < min_dist) {
				min_dist = dist;
				DOA_class[j] = i;
			}
		}
	}
}