#include "unwrap.h"

double unwrap (int index, int N, double max_index) {
	double temp;

	if (index < N)
		temp = -((double) index);
	else
		temp = ((double) 2*N - index);

	if (temp > max_index)
		return max_index;
	if (temp < -max_index)
		return -max_index;

	return temp;
}