#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
    cout << "Error with input vector sizes in CalculateRMSE()" << endl;
    return rmse;
  }
  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    VectorXd aux = (estimations[i]-ground_truth[i]);
    aux = aux.array() * aux.array();
    rmse += aux;
  }

  //calculate the mean
  rmse /= estimations.size();
  //calculate the squared root
  rmse = rmse.array().sqrt();
  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

	//recover state parameters
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	//check division by zero
	if(px == 0 && py == 0) {
    cout << "Error: Division by zero in CalculateJacobian()" << endl;
    return MatrixXd::Zero(3, 4);
	}

	//compute the Jacobian matrix
  double a = px*px + py*py;
  double b = sqrt(a);
  double c = a * b;
  double d = (vx*py - vy*px) / c;

  MatrixXd Hj(3,4);
  Hj <<  px / b,  py / b, 0,      0,
        -py / a,  px / a, 0,      0,
         py * d, -px * d, px / b, py / b;

	return Hj;
}

double Tools::SubtractAngles(double a, double b) {
    double diff = a - b;
    if (diff < -M_PI) {
        diff += 2 * M_PI;
    }
    else if (diff > M_PI) {
      diff -= 2 * M_PI;
    }
    return diff;
}
