#ifndef KMEANS
#define KMEANS

void kmeans (double* DOA_hist, double* DOA_mean, unsigned int* class_count, unsigned int n_sources, unsigned int hist_length);
void kmeans (double* DOA_hist, unsigned int *DOA_class, double* DOA_mean, unsigned int* class_count, unsigned int n_sources, unsigned int hist_length);
void kmeans (double* DOA_hist, unsigned int *DOA_class, double* DOA_mean, unsigned int n_sources, unsigned int hist_length);

#endif