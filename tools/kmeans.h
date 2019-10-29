#ifndef KMEANS
#define KMEANS

void kmeans (double* DOA_hist, double* DOA_mean, int* class_count, int n_sources, int hist_length);
void kmeans (double* DOA_hist, int *DOA_class, double* DOA_mean, int* class_count, int n_sources, int hist_length);
void kmeans (double* DOA_hist, int *DOA_class, double* DOA_mean, int n_sources, int hist_length);

#endif