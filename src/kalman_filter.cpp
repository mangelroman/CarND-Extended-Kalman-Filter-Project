#include "kalman_filter.h"

#include <math.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(const VectorXd &x_in, const MatrixXd &P_in) {
  x_ = x_in;
  P_ = P_in;
}

void KalmanFilter::Predict(const MatrixXd &F_in, const MatrixXd &Q_in) {
  x_ = F_in * x_;
	P_ = F_in * P_ * F_in.transpose() + Q_in;
}

void KalmanFilter::Update(const VectorXd &z, const MatrixXd &H_in, const MatrixXd &R_in) {
	VectorXd y = z - H_in * x_;

  _update(y, H_in, R_in);
}

void KalmanFilter::UpdateEKF(const VectorXd &z, const MatrixXd &H_in, const MatrixXd &R_in) {
  double px = x_(0);
	double py = x_(1);
	double vx = x_(2);
	double vy = x_(3);

	//check division by zero
	if(px == 0 && py == 0) {
	    cout << "Error: Division by zero in UpdateEKF()" << endl;
	    return;
	}

  VectorXd z_pred(3);
  double rho = sqrt(px*px + py*py);
  double phi = atan2(py, px);

  z_pred << rho,
            phi,
            (rho != 0) ? (px*vx + py*vy) / rho : 0;

  VectorXd y = z - z_pred;

  y(1) = tools.SubtractAngles(z(1), phi);

  _update(y, H_in, R_in);
}

void KalmanFilter::_update(const VectorXd &y, const MatrixXd &H, const MatrixXd &R) {
	MatrixXd Ht = H.transpose();
	MatrixXd S = H * P_ * Ht + R;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H) * P_;
}
