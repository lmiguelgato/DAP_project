#ifndef KALMAN
#define KALMAN

// Eigen: C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
#include <Eigen/Dense>
using namespace Eigen;

void kalman (double y_array[2], double *x_array, double (*P_array)[4]);

#endif