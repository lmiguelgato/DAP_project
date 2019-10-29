#ifndef KALMAN
#define KALMAN

// Eigen: C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
#include <Eigen/Dense>
using namespace Eigen;

template <typename Derived, typename OtherDerived>
void kalman (const MatrixBase<Derived>& y, MatrixBase<OtherDerived> const & x, MatrixBase<OtherDerived> const & P);

#endif