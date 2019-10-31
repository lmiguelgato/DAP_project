#include "kalman.h"

#define TAU 0.0213333	// 1024/48000: observation time
#define ACC 1.0			// modeled acceleration
#define NVAR 0.005	// measurement noise variance
//#define ACC 0.05		// modeled acceleration
//#define NVAR 0.000001	// measurement noise variance

// state transition model:
const MatrixXd F = (MatrixXd(4,4) << 1.0, 0.0, TAU, 0.0, 0.0, 1.0, 0.0, TAU, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0).finished();
const MatrixXd Q = (MatrixXd(4,4) << ACC*TAU*TAU*TAU*TAU/4, 0.0, ACC*TAU*TAU*TAU/2, 0.0, 0.0, ACC*TAU*TAU*TAU*TAU/4, 0.0, ACC*TAU*TAU*TAU/2, ACC*TAU*TAU*TAU/2, 0.0, ACC*TAU*TAU, 0.0, 0.0, ACC*TAU*TAU*TAU/2, 0.0, ACC*TAU*TAU).finished();
// measurement model:
const MatrixXd H = (MatrixXd(2,4) << 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0).finished();
const MatrixXd R = (MatrixXd(2,2) << NVAR, 0.0, 0.0, NVAR).finished();

// 4x4 identity matrix:
const MatrixXd I4x4 = (MatrixXd(4,4) << 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0).finished();

void kalman (double y_array[2], double *x_array, double (*P_array)[4]) {

	/*
		Kalman filter

	 	input:	x   (previous state)
	    		y   (current measurement)
				P   (previous-state error posterior covariance matrix)
	*/

	// initialization:
	MatrixXd y(2,1);
	MatrixXd x(4,1);
	MatrixXd P(4,4);

	y(0,0) = y_array[0];	y(1,0) = y_array[1];

	for (int i = 0; i < 4; ++i)
	{
		x(i,0) = x_array[i];

		for (int j = 0; j < 4; ++j)
		{
			P(i,j) = P_array[i][j];
		}
	}

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

	x = temp4x1 + K*(y - H*temp4x1);
	P = (I4x4 - K*H)*P_;

	for (int i = 0; i < 4; ++i)
	{
		x_array[i] = x(i,0);

		for (int j = 0; j < 4; ++j)
		{
			P_array[i][j] = P(i,j);
		}
	}
}