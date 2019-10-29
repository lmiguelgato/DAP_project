#include "kalman.h"

#define TAU 0.0427	// 2048/48000: observation time
#define ACC 1.0		// modeled acceleration
#define NVAR 0.001	// measurement noise variance

// state transition model:
const MatrixXd F = (MatrixXd(4,4) << 1.0, 0.0, TAU, 0.0, 0.0, 1.0, 0.0, TAU, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0).finished();
const MatrixXd Q = (MatrixXd(4,4) << ACC*TAU*TAU*TAU*TAU/4, 0.0, ACC*TAU*TAU*TAU/2, 0.0, 0.0, ACC*TAU*TAU*TAU*TAU/4, 0.0, ACC*TAU*TAU*TAU/2, ACC*TAU*TAU*TAU/2, 0.0, ACC*TAU*TAU, 0.0, 0.0, ACC*TAU*TAU*TAU/2, 0.0, ACC*TAU*TAU).finished();
// measurement model:
const MatrixXd H = (MatrixXd(2,4) << 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0).finished();
const MatrixXd R = (MatrixXd(2,2) << NVAR, 0.0, 0.0, NVAR).finished();

// 4x4 identity matrix:
const MatrixXd I4x4 = (MatrixXd(4,4) << 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0).finished();

template <typename Derived, typename OtherDerived>
void kalman (const MatrixBase<Derived>& y, MatrixBase<OtherDerived> const & x, MatrixBase<OtherDerived> const & P) {

	/*
		Kalman filter

	 	input:	x   (previous state)
	    		y   (current measurement)
				P   (previous-state error posterior covariance matrix)
	*/

	MatrixXd P_(4,4);		// previous-state error prior covariance matrix
	MatrixXd K(4,2);		// Kalman gain
	MatrixXd temp2x2(2,2);	// temporary 2x2 matrix
	MatrixXd temp4x1(4,1);	// temporary 4x1 matrix

	// prediction:
	temp4x1	= F*x;
	P_ = F*P*F.transpose() + Q;

	// update:
	temp2x2 = H*P_*H.transpose() + R;
	K = (P_ * H.transpose()) * temp2x2.inverse();

	const_cast< MatrixBase<OtherDerived>& >(x) = temp4x1 + K*(y - H*temp4x1);

	const_cast< MatrixBase<OtherDerived>& >(P) = (I4x4 - K*H)*P_;
}