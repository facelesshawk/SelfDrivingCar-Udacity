#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  MatrixXd FT_ = F_.transpose();
  x_ = F_*x_;
  P_ = (F_*P_*FT_) + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  MatrixXd Y_ = z-(H_*x_);
  MatrixXd HT_ = H_.transpose();
  MatrixXd S_ = (H_*P_*HT_)+R_;
  MatrixXd SI_ = S_.inverse();
  MatrixXd K_ = P_*HT_*SI_;
  x_ = x_+(K_*Y_);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I-(K_*H_))*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  VectorXd hj = VectorXd(3);
  hj(0) = sqrt( x_(0) *  x_(0) + x_(1) * x_(1) );
  hj(1) = atan2( x_(1) , x_(0) );
  hj(2) = (x_(0) * x_(2) + x_(1) * x_(3)) / hj(0);
  VectorXd Y_ = z - hj;
  Y_(1) = atan2(sin(Y_(1)),cos(Y_(1)));
  MatrixXd HT_ = H_.transpose();
  MatrixXd S_ = (H_*P_*HT_)+R_;
  MatrixXd SI_ = S_.inverse();
  MatrixXd K_ = P_*HT_*SI_;
  x_ = x_+(K_*Y_);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I-(K_*H_))*P_;
}
