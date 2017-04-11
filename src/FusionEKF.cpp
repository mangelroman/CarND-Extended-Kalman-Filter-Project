#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {

  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // H_radar_ will be calculated on every measurement update using the
  // CalculateJacobian() helper function
  H_radar_ = MatrixXd(3, 4);

  // process noise covariance matrix
  double noise_ax = 9, noise_ay = 9;
  Qv_ = MatrixXd(2, 2);
  Qv_ << noise_ax, 0,
         0,        noise_ay;

  // measurement noise covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  // measurement noise covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    VectorXd x_init(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double rho   = measurement_pack.raw_measurements_[0];
      double phi   = measurement_pack.raw_measurements_[1];
      double d_rho = measurement_pack.raw_measurements_[2];

      x_init << rho * cos(phi), rho * sin(phi), d_rho * cos(phi), d_rho * sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      double px = measurement_pack.raw_measurements_[0];
      double py = measurement_pack.raw_measurements_[1];

      x_init << px, py, 0, 0;
    }

    MatrixXd P_init(4, 4);
    P_init << 1, 0, 0,    0,
              0, 1, 0,    0,
              0, 0, 1000, 0,
              0, 0, 0,    1000;

    ekf_.Init(x_init, P_init);
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
   //Modify the F matrix so that the time is integrated
   double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0f;
   previous_timestamp_ = measurement_pack.timestamp_;

   MatrixXd F(4, 4);
   F << 1,  0, dt,  0,
        0,  1,  0, dt,
        0,  0,  1,  0,
        0,  0,  0,  1;

  double dt2 = (dt * dt) / 2;

  MatrixXd G(4, 2);
  G << dt2, 0,
       0,   dt2,
       dt,  0,
       0,   dt;

  MatrixXd Q = G * Qv_ * G.transpose();

  ekf_.Predict(F, Q);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  // Update the state and covariance matrices.
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    //printf("P:%8.4f %8.4f %8.4f %8.4f | ", ekf_.x_(0), ekf_.x_(1), ekf_.x_(2), ekf_.x_(3));
    H_radar_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_, H_radar_, R_radar_);
    //printf("U:%8.4f %8.4f %8.4f %8.4f\n", ekf_.x_(0), ekf_.x_(1), ekf_.x_(2), ekf_.x_(3));
    // Radar updates
  } else {
    ekf_.Update(measurement_pack.raw_measurements_, H_laser_, R_laser_);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
