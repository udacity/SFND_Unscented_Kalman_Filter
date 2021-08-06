#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.30;
  
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
  lambda_ = 3.0 - n_x_;
  n_sigma_ = 2*n_x_ + 1;
  weights_.push_back(lambda_ / (n_x_ + lambda_));
  for (int i=1; i<n_sigma_; ++i){
    weights_.push_back(0.5 / (n_x_ + lambda_));
  }
  time_us_ = 0;
  Xsig_pred_ = MatrixXd( n_x_, n_sigma_);
  Xsig_pred_.fill(0.0);
  R_Laser = MatrixXd(2, 2);
  R_Laser << std_laspx_*std_laspx_, 0,
             0, std_laspy_*std_laspy_;
  R_Radar = MatrixXd(3, 3);
  R_Radar << std_radr_*std_radr_, 0,                       0,
             0,                   std_radphi_*std_radphi_, 0,
             0,                   0,                       std_radrd_*std_radrd_;
  
  std::cout << "UKF constructed.\n";
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (is_initialized_){
    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    Prediction(dt);
    switch (meas_package.sensor_type_){
    case MeasurementPackage::LASER:{
      if (use_laser_){ UpdateLidar(meas_package); }
      break;
    }
    case MeasurementPackage::RADAR:
    { 
      if (use_radar_) { UpdateRadar(meas_package); }
      break;
    }
    default:
      std::cout << "Measurement not recognized.\n";
      break;
    }
  }
  else{
    //double x, y;
    switch (meas_package.sensor_type_)
    {
      case MeasurementPackage::LASER:
      {
        std::cout << "Initializing off laser measurement!\n";
        x_ << meas_package.raw_measurements_[0],
              meas_package.raw_measurements_[1],
              0,
              0,
              0;
        P_ = Eigen::MatrixXd::Identity(5,5);
        P_.topLeftCorner(2,2) = R_Laser;
        is_initialized_ = true;
        time_us_ = meas_package.timestamp_;
      }
      case MeasurementPackage::RADAR:
      {
        std::cout << "Initializing off radar measurement!\n";
        double rho = meas_package.raw_measurements_(0);
        double phi = meas_package.raw_measurements_(1);
        x_ << rho * cos(phi),
              rho * sin(phi),
              0,
              0,
              0;
        P_ = Eigen::MatrixXd::Identity(5,5);
        P_.topLeftCorner(3,3) = R_Radar;
        is_initialized_ = true;
        time_us_ = meas_package.timestamp_;
      }
      default:
        std::cout << "Measurement not recognized. Initialization incomplete.\n";
        break;
    }
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  std::cout << "Predict\n";
  /* 
    Unfortunately, we can't do straight matrix summations with Eigne in C++ source code
    like what can be done in the command prompt on Matlab (or in the Python bindings, 
    I think). So we need to decompose what a matrix and vector sum actually is: a 
    column-wise sum of the vector plus the column of the matrix.
  */
  // Generate augmented sigma points
  VectorXd x_aug = VectorXd(7);                   // Augmented state estimate
  x_aug << x_, 0, 0;
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);      // Augmented covariance
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;
  MatrixXd L = P_aug.llt().matrixL();             // Approximate matrix square root
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_); // Augmented sigma points
  Xsig_aug.fill(0.0);
  Xsig_aug.col(0) = x_aug;
  std::cout << "Generating sigma points...";
  for (int i = 0; i < n_x_; ++i)
  {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_x_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  std::cout << "Done!\n";  
  // Propagate each sigma point
  std::cout << "Propagating sigma points...";
  for (int i = 0; i < n_sigma_; ++i)
  {
    // extract values for better readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001)
    {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else
    {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }
    // Propagate and add noise
    Xsig_pred_(0, i) = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    Xsig_pred_(1, i) = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);;
    Xsig_pred_(2, i) = v + nu_a * delta_t;
    Xsig_pred_(3, i) = yaw + yawd * delta_t + 0.5 * nu_yawdd * delta_t * delta_t;
    Xsig_pred_(4, i) = yawd + nu_yawdd*delta_t;
  }
  std::cout << "Done!\n";
  // Determine predicted state mean and covariance
  std::cout << "Predicting state mean...";
  x_.fill(0.0);
  for (int i=0; i<n_sigma_; ++i) { x_ += weights_[i]*Xsig_pred_.col(i); }
  std::cout << "Done!\n";
  std::cout << "Predicting covariance...";
  P_.fill(0.0);
  for (int i=0; i<n_sigma_; ++i){
    VectorXd diff = Xsig_pred_.col(i) - x_;
    // Normalize the yaw angle
    while (diff[3] > M_PI)  { diff[3] -= 2.0 * M_PI; }
    while (diff[3] < -M_PI) { diff[3] += 2.0 * M_PI; }
    P_ += weights_[i] * diff*diff.transpose();
  }  
  std::cout << "Done!\n";
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  std::cout << "Lidar update\n";
  int n_z = 2;
  // Create measurement sigma points and predict the measurement
  Eigen::MatrixXd Zsig(n_z, n_sigma_);
  Zsig.fill(0.0);
  Eigen::VectorXd z_pred(n_z);
  for (int i=0; i<n_sigma_; ++i){
    Zsig(0,i) = Xsig_pred_(0,i);
    Zsig(1,i) = Xsig_pred_(1,i);
    z_pred += weights_[i]*Zsig.col(i);
  }
  // Create innovation covariance matrix S
  Eigen::MatrixXd S(n_z, n_z);
  S.fill(0.0);
  for (int i=0; i<n_sigma_; ++i){
    Eigen::VectorXd diff = Zsig.col(i) - z_pred;
    // Angle normalization
    //while (diff[1] > M_PI)  { diff[1] -= 2.0 * M_PI; }
    //while (diff[1] < -M_PI) { diff[1] += 2.0 * M_PI; }
    S += weights_[i]*diff*diff.transpose();
  }
  // Add noise
  S += R_Laser;
  // Calculate cross correlation matrix T
  Eigen::MatrixXd T(n_x_, n_z);
  T.fill(0.0);
  for (int i=0; i<n_sigma_; ++i){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //while (z_diff[1] > M_PI)  { z_diff[1] -= 2.0*M_PI; }
    //while (z_diff[1] < -M_PI) { z_diff[1] += 2.0*M_PI; }
    //std::cout << "Zdiff:"<<z_diff.transpose()<<"\n------\n";
    Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //while (x_diff[3] > M_PI)  { x_diff[3] -= 2.0*M_PI; }
    //while (x_diff[3] < -M_PI) { x_diff[3] += 2.0*M_PI; }
    //std::cout << "Xdiff: "<<x_diff.transpose()<<"\n-----\n";
    //std::cout << x_diff * z_diff.transpose();
    T = T + weights_[i] * x_diff * z_diff.transpose();
    //std::cout << "T:\n"<<T<<"\n";
  }
  //std::cout<<"MARKER\n";
  // Calculate Kalman gain
  Eigen::MatrixXd K = T * S.inverse();
  // Update the mean and covariance
  Eigen::VectorXd z = meas_package.raw_measurements_;
  x_ += K*(z - z_pred);
  P_ -= K*S*K.transpose();

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  std::cout << "Radar update\n";
  int n_z_ = 3; // Measurement vector size
  // Create measurement sigma points in measurement coordinates and predict the measurement
  Eigen::MatrixXd Zsig(n_z_, n_sigma_); // Measurement sigma points
  Zsig.fill(0.0);
  Eigen::VectorXd z_pred(n_z_);
  z_pred.fill(0.0);
  //std::cout << "marker\n";
  for (int i=0; i<n_sigma_; ++i){
    // Extract values
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    //std::cout << "Sigma point: ["<<p_x<<", "<<p_y<<", "<<v<<", "<<yaw<<"]\n";
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);    // r
    Zsig(1,i) = atan2(p_y, p_x);            // phi
    Zsig(2,i) = (p_x*cos(yaw)*v + p_y*sin(yaw)*v) / sqrt(p_x*p_x + p_y*p_y);
    z_pred += weights_[i]*Zsig.col(i);
  }
  
  // Create innovation covariance matrix S
  Eigen::MatrixXd S(n_z_, n_z_);
  S.fill(0.0);
  for (int i=0; i<n_sigma_; ++i){
    VectorXd diff = Zsig.col(i) - z_pred;
    // Angle normalization
    while (diff[1] > M_PI)  { diff[1] -= 2.0 * M_PI; }
    while (diff[1] < -M_PI) { diff[1] += 2.0 * M_PI; }
    S += weights_[i]*diff*diff.transpose();
  }
  // Add noise
  S += R_Radar;
  // Calculate cross correlation matrix T
  Eigen::MatrixXd T(n_x_, n_z_);
  T.fill(0.0);
  for (int i=0; i<n_sigma_; ++i){
    Eigen::VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff[1] > M_PI)  { z_diff[1] -= 2.0*M_PI; }
    while (z_diff[1] < -M_PI) { z_diff[1] += 2.0*M_PI; }
    //std::cout << "Zdiff:"<<z_diff.transpose()<<"\n------\n";
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff[3] > M_PI)  { x_diff[3] -= 2.0*M_PI; }
    while (x_diff[3] < -M_PI) { x_diff[3] += 2.0*M_PI; }
    //std::cout << "Xdiff: "<<x_diff.transpose()<<"\n-----\n";
    //std::cout << x_diff * z_diff.transpose();
    T = T + weights_[i] * x_diff * z_diff.transpose();
    //std::cout << "T:\n"<<T<<"\n";
  }
  
  //std::cout<<"MARKER\n";
  // Calculate Kalman gain
  Eigen::MatrixXd K = T * S.inverse();
  // Update the mean and covariance
  Eigen::VectorXd z = meas_package.raw_measurements_;
  x_ += K*(z - z_pred);
  P_ -= K*S*K.transpose();
}