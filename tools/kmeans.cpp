#include "kmeans.h"
#include <math.h>

void kmeans (double* DOA_hist, double* DOA_mean, int* class_count, int n_sources, int hist_length) {
	int i, j;
	double min_dist, dist;

	int *DOA_class;
	DOA_class = (int *) calloc(hist_length, sizeof(int));

	for (i = 0; i < n_sources; ++i)
	{
		class_count[i] = 0;
	}

	for (j = 0; j < hist_length; ++j)
	{
		min_dist = 360.0;
		for (i = 0; i < n_sources; ++i)
		{			
			dist = abs(DOA_hist[j] - DOA_mean[i]);
			if (dist < min_dist) {
				min_dist = dist;
				DOA_class[j] = i;
			}
		}
		++class_count[ DOA_class[j] ];
	}

	for (i = 0; i < n_sources; ++i)
	{			
		DOA_mean[i] = 0.0;

		for (j = 0; j < hist_length; ++j)
		{
			if 	(DOA_class[j] == i) {
				DOA_mean[i] += DOA_hist[j];
			}
		}
		if (class_count[i] != 0)
			DOA_mean[i] /= class_count[i];		
	}
}