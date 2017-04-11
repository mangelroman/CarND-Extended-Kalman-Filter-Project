#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#include "tools.h"
#include "Eigen/Dense"

class KalmanFilter {
public:

  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // state transistion matrix
  Eigen::MatrixXd F_;

  // process covariance matrix
  Eigen::MatrixXd Q_;

  // measurement matrix
  Eigen::MatrixXd H_;

  // measurement covariance matrix
  Eigen::MatrixXd R_;

  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param x_in Initial state
   * @param P_in Initial state covariance
   * @param F_in Transition matrix
   * @param H_in Measurement matrix
   * @param R_in Measurement covariance matrix
   * @param Q_in Process covariance matrix
   */
  void Init(const Eigen::VectorXd &x_in, const Eigen::MatrixXd &P_in);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict(const Eigen::MatrixXd &F_in, const Eigen::MatrixXd &Q_in);

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(
    const Eigen::VectorXd &z,
    const Eigen::MatrixXd &H_in,
    const Eigen::MatrixXd &R_in);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(
    const Eigen::VectorXd &z,
    const Eigen::MatrixXd &H_in,
    const Eigen::MatrixXd &R_in);

private:
  void _update(
    const Eigen::VectorXd &y,
    const Eigen::MatrixXd &H,
    const Eigen::MatrixXd &R);

  Tools tools;
};


#endif /* KALMAN_FILTER_H_ */
