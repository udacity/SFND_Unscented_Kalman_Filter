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
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  is_initialized_ = false;
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
        break;
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
        break;
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
  //Eigen::MatrixXd X_sig = GenerateSigmaPoints();
  GenerateSigmaPoints();
  
}

void UKF::GenerateSigmaPoints(){
  Eigen::MatrixXd P_sqrt = P_.llt().matrixL();  
  /* 
    Unfortunately, we can't do straight matrix summations with Eigne in C++ source code
    like what can be done in the command prompt on Matlab (or in the Python bindings, 
    I think). So we need to decompose what a matrix and vector sum actually is: a 
    column-wise sum of the vector plus the column of the matrix.
  */
  Xsig_pred_.col(0) << x_;
  for (int i=1; i<n_x_; ++i){
    Xsig_pred_.col(i)          = x_ + (sqrt(lambda_ + n_x_)) * P_sqrt.col(i);
    Xsig_pred_.col(i+n_x_) = x_ - (sqrt(lambda_ + n_x_)) * P_sqrt.col(i);
  }  
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  std::cout << "Lidar update\n";
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  std::cout << "Radar update\n";
}