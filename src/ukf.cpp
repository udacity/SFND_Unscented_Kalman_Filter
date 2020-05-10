#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // Set UKF parameters

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // Augmented sigma point spreading parameter
  lambda_aug_ = 3 - n_aug_;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial augmented state vector
  x_aug_ = VectorXd(n_aug_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // initial augmented covariance matrix
  P_aug_ = MatrixXd(n_aug_, n_aug_);

  // initial sigma point matrix
  Xsigma_ = MatrixXd(n_x_, 2 * n_x_ + 1);

  // initial augmented sigma point matrix
  Xsigma_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // initial predicted sigma point matrix
  Xsigma_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values
   */
}

UKF::~UKF() {}

void UKF::GenerateSigmaPoints()
{
  MatrixXd A = P_.llt().matrixL();
  Xsigma_.col(0) = x_;
  for (int i = 0; i < n_x_; ++i)
  {
    Xsigma_.col(i + 1) = x_ + sqrt(lambda_ + n_x_) * A.col(i);
    Xsigma_.col(i + 1 + n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
  }
}

void UKF::GenerateAugmentedSigmaPoints()
{
  // add augmented part to mean state
  x_aug_.head(n_x_) = x_;
  x_aug_(n_x_) = std_a_;
  x_aug_(n_x_ + 1) = std_yawdd_;

  // add augmented part to covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(n_x_, n_x_) = std_a_ * std_a_;
  P_aug_(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  MatrixXd A_aug = P_aug_.llt().matrixL();
  Xsigma_aug_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; ++i)
  {
    Xsigma_aug_.col(i + 1) = x_aug_ + sqrt(lambda_aug_ + n_aug_) * A_aug.col(i);
    Xsigma_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_aug_ + n_aug_) * A_aug.col(i);
  }
}

void UKF::PredictSigmaPoints(const double dt)
{
  double epsilon = 1e-5; // equal to zero threshold for yaw_rate
  for (int i = 0; i < Xsigma_pred_.cols(); i++)
  {
    // Read data from augmented sigma point
    double pos_x = Xsigma_aug_(0, i);
    double pos_y = Xsigma_aug_(1, i);
    double vel_abs = Xsigma_aug_(2, i);
    double yaw_angle = Xsigma_aug_(3, i);
    double yaw_rate = Xsigma_aug_(4, i);
    double std_a = Xsigma_aug_(5, i);
    double std_yawdd = Xsigma_aug_(6, i);

    // Sigma point that will be predicted
    double pos_x_pred;
    double pos_y_pred;
    double vel_abs_pred;
    double yaw_angle_pred;
    double yaw_rate_pred;

    // Apply CTRV motion model
    if (std::abs(yaw_rate) < epsilon) // if yaw_rate is 0
    {
      pos_x_pred = pos_x + vel_abs * dt * std::cos(yaw_angle);
      pos_y_pred = pos_y + vel_abs * dt * std::sin(yaw_angle);
    }
    else
    {
      pos_x_pred = pos_x + (vel_abs / yaw_rate) * (std::sin(yaw_angle + yaw_rate * dt) - std::sin(yaw_angle));
      pos_y_pred = pos_y + (vel_abs / yaw_rate) * (std::cos(yaw_angle) - std::cos(yaw_angle + yaw_rate * dt));
    }
    vel_abs_pred = vel_abs;
    yaw_angle_pred = yaw_angle + yaw_rate * dt;
    yaw_rate_pred = yaw_rate;

    // Add noise
    pos_x_pred = pos_x_pred + 0.5 * dt * dt * std::cos(yaw_angle) * std_a;
    pos_y_pred = pos_y_pred + 0.5 * dt * dt * sin(yaw_angle) * std_a;
    vel_abs_pred = vel_abs_pred + std_a * dt;
    yaw_angle_pred = yaw_angle_pred + 0.5 * dt * dt * std_yawdd;
    yaw_rate_pred = yaw_rate_pred + std_yawdd * dt;

    // Write the predicted sigma point
    Xsigma_pred_(0, i) = pos_x_pred;
    Xsigma_pred_(1, i) = pos_y_pred;
    Xsigma_pred_(2, i) = vel_abs_pred;
    Xsigma_pred_(3, i) = yaw_angle_pred;
    Xsigma_pred_(4, i) = yaw_rate_pred;
  }
}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
}

void UKF::Prediction(double dt)
{
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */

  GenerateSigmaPoints();
  GenerateAugmentedSigmaPoints();
}

void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}