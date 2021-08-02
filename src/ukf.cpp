#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3.0;
  
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
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  n_x_ = 5;
  n_aug_ = 7;
  n_sigma_ = 2*n_x_ - 1;
  lambda_ = 3.0 - n_aug_;
  weights_[0] = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i<weights_.size(); ++i){
    weights_[i] = 1.0 / (2.0*(lambda_ + n_aug_));
  }
  Xsig_pred_ = MatrixXd((int) n_x_, (int) n_sigma_);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (is_initialized_) {
    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    Prediction(dt);
    switch (meas_package.sensor_type_)
    {
      case MeasurementPackage::LASER:
        if (use_laser_) { UpdateLidar(meas_package); }
        break;
      case MeasurementPackage::RADAR:
        if (use_radar_) { UpdateRadar(meas_package); }
        break;
      default:
        break;
    }
  }
  else {
    double x, y;
    switch (meas_package.sensor_type_)
    {
      case MeasurementPackage::LASER:
        x = meas_package.raw_measurements_[0];
        y = meas_package.raw_measurements_[1];
        P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
              0, std_laspy_ * std_laspy_, 0, 0, 0,
              0, 0, 1, 0, 0,
              0, 0, 0, 1, 0,
              0, 0, 0, 0, 1;
        x_ << x, y, 0, 0, 0;
        break;
      case MeasurementPackage::RADAR:
      
        auto Rho = meas_package.raw_measurements_(0);
        auto Phi = meas_package.raw_measurements_(1);
        float Rho_dot = meas_package.raw_measurements_(2);
        x = Rho * cos(Phi);
        y = Rho * sin(Phi);
        P_ << std_laspy_ * std_radr_, 0, 0, 0, 0,
              0, std_radr_ * std_radr_,  0, 0, 0,
              0, 0, std_radrd_ * std_radrd_, 0, 0,
              0, 0, 0, std_radphi_ * std_radphi_, 0,
              0, 0, 0, 0, 1;
        //x_ << x, y, 0, 0, 0;
      
        break;
      //default:
      //  exit(-1);
      //  break;
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
  }  
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  /* Step 1: Augment the state vector and covariance matricies with noise parameters */
  VectorXd x_a = VectorXd(n_aug_);
  x_a << x_, 
         0, 
         0;
  MatrixXd P_a = MatrixXd(n_aug_, n_aug_);
  P_a.fill(0.0);
  P_a.topLeftCorner(n_x_, n_x_) = P_;
  P_a(5,5) = std_a_;
  P_a(6,6) = std_yawdd_;
  
  /* Step 2: Generate sigma points */
  MatrixXd X_a = GenerateSigmaPoints(x_a, P_a);

  /* Step 3: Propagate the sigma points */
  for (int i=0; i<X_a.cols(); ++i){
    VectorXd x_col = X_a.col(i);
    Xsig_pred_.col(i) = Propagate(x_col, delta_t);
  }

  /* Step 4: Predict mean and covariance */
  x_.fill(0.0); // Reset
  for (int i=0; i<n_sigma_; ++i){
    x_ += weights_[i] * Xsig_pred_.col(i); // weighted average of sigma points
  }
  P_.fill(0.0); // Reset
  for (int i=0; i<n_sigma_; ++i){
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // Angle normalizations
    while (x_diff[3] > M_PI) { x_diff[3] -= 2.*M_PI; }
    while (x_diff[3] <-M_PI) { x_diff[3] += 2.*M_PI; }
    P_ += weights_[i] * x_diff * x_diff.transpose();
  }
}

Eigen::VectorXd UKF::Propagate(Eigen::VectorXd& x, double& delta_t){
  // UKF convention is to have the state estimate and noise parameters combined into one 'augmented' state
  Eigen::VectorXd x_k, nu_k;
  nu_k << 0.5*x[5]*delta_t*delta_t*cos(x[3]),
          0.5*x[6]*delta_t*delta_t*sin(x[3]),
          x[5]*delta_t,
          0.5*x[6]*delta_t*delta_t,
          x[6]*delta_t;

  if ( fabs(x(4)) > 0.001){
    x_k << x[0] + x[2]/x[4] * (sin(x[3] + x[4]*delta_t) - sin(x[3])),
           x[1] + x[2]/x[4] * (cos(x[3]) - cos(x[3] + x[4]*delta_t)),
           x[2],
           x[3] + x[4]*delta_t,
           x[4];
  }
  else {
    x_k << x[0] + x[2]*delta_t*cos(x[3]);
           x[1] + x[2]*delta_t*sin(x[3]);
           x[2],
           x[3] + x[4]*delta_t,
           x[4];
  }
  return x_k + nu_k;
}

Eigen::MatrixXd UKF::GenerateSigmaPoints(Eigen::VectorXd& x_augmented, Eigen::MatrixXd& P_augmented){
  MatrixXd sigma_points = MatrixXd(x_augmented.rows()-2, x_augmented.cols());
  sigma_points.col(0) = x_augmented.col(0);
  MatrixXd L = P_augmented.llt().matrixL();
  for (int i = 0; i < n_aug_; ++i){
    sigma_points.col(i+1) = x_augmented + sqrt(lambda_ + n_aug_) * L.col(i);
    sigma_points.col(i+1+n_aug_) = x_augmented + sqrt(lambda_ + n_aug_) * L.col(i);
  }
  return sigma_points;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}