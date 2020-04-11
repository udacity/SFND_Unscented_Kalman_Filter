#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <exception>
#include <math.h>       /* isnan, sqrt */


class unknownsensortype: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Unknown sensor type";
  }
} sensor_exception;

using Eigen::MatrixXd;
using Eigen::VectorXd;

void normalizeAngle(double& angle){
  while (angle> M_PI) angle -= 2.*M_PI;
  while (angle<-M_PI) angle += 2.*M_PI;  
}
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
  std_a_ = 2.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;
  
  std::cout << "[std_a_, std_yawdd_] = [" << std_a_ << ", " << std_yawdd_<<"]\n";
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
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0,std_radrd_*std_radrd_;
  
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_*std_laspx_,0,
            0,std_laspy_*std_laspy_;  
  

  is_initialized_ = false;  
  n_x_    = 5;
  n_aug_  = 7;
  lambda_ = 3 - n_aug_;
  weights_ = VectorXd(2*n_aug_+1);
  const double weight_0 = lambda_/(lambda_ + n_aug_);
  const double weight   = 0.5/(lambda_ + n_aug_);
  weights_(0) = weight_0;
  for(int i=1; i < 2*n_aug_+1; ++i){
    weights_(i) = weight;
  }

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1); 
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  if(is_initialized_){
    const double delta_t = (meas_package.timestamp_ - time_us_)/1e6;    
    Prediction(delta_t);    
  }
  time_us_ = meas_package.timestamp_;
  switch (meas_package.sensor_type_)
  {
  case MeasurementPackage::SensorType::LASER:    
    UpdateLidar(meas_package);
    break;

  case MeasurementPackage::SensorType::RADAR:    
    UpdateRadar(meas_package);
    break;
    
  default:
    throw sensor_exception;
    break;
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  //if(fabs(delta_t)<1e-4) return;
  //std::cout << "Generate sigma points" << delta_t << std::endl;
  
  MatrixXd Xsig_aug;
  Xsig_aug.fill(0.0);
  
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_aug_-2,n_aug_-2) = std_a_ * std_a_;
  P_aug(n_aug_-1,n_aug_-1) = std_yawdd_ * std_yawdd_; 
  VectorXd x = VectorXd(n_aug_);
  x.fill(0.);
  x.topRows(n_x_) = x_;  
  UKF::GenerateSigmaPoints(&Xsig_aug, x, P_aug, lambda_);
  

  // Predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    const double& p_x      = Xsig_aug(0,i);
    const double& p_y      = Xsig_aug(1,i);
    const double& v        = Xsig_aug(2,i);
    const double& yaw      = Xsig_aug(3,i);
    const double& yawd     = Xsig_aug(4,i);
    const double& nu_a     = Xsig_aug(5,i);
    const double& nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw + yawd*delta_t/2);
        py_p = p_y + v*delta_t*sin(yaw + yawd*delta_t/2);
    }

    double v_p    = v;
    double yaw_p  = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p  = v_p  + nu_a*delta_t;

    yaw_p  = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // Restore predicted mean and cov from sigma points  
  x_.fill(0.0);
  double cs = 0;
  double ss = 0;
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    VectorXd x = Xsig_pred_.col(i);
    ss += weights_(i)*sin(x(3));
    cs += weights_(i)*cos(x(3));
    x_ = x_ + weights_(i) * x;
  }   
  x_(3) = atan2(ss,cs);

  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }  
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */  
    
  const int n_z = 2;
  VectorXd z = meas_package.raw_measurements_;
  if(!is_initialized_){    
    // Init state and covariance matrix
    x_.fill(0.);
    P_.fill(0.);
    P_(0,0) = P_(1,1) = 1*1; // Init to a std of 1 m
    P_(2,2) = 10 * 10; // Init relative speed std to 10 m/s
    P_(3,3) = 1; // Angular std to 180 degree
    P_(4,4) = 1; // Angular speed to 180 degree/sec

    x_.fill(0);
    x_.topRows(2) = z;
    Prediction(0);
    P_.topLeftCorner(2,2) = R_lidar_;    
    is_initialized_ = true;
    return;
  }

  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  

    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R_lidar_;
  

  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    VectorXd z_diff = Zsig.col(i) - z_pred;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd K = Tc * S.inverse();

  VectorXd y = z - z_pred;

  x_ = x_ + K * y;
  P_ = P_ - K*S*K.transpose();  

  // Compute consistency
  std::cout << "0, " << y.transpose()*S.inverse()*y << std::endl;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
    
  static const int n_z = 3;
  VectorXd z = meas_package.raw_measurements_;

  if(!is_initialized_){    
    // Init state and covariance matrix
    x_.fill(0.);
    P_.fill(0.);
    P_(0,0) = P_(1,1) = 1*1; // Init to a std of 1 m
    P_(2,2) = 1; // Init relative speed std to 10 m/s
    P_(3,3) = 1; // Angular std to 180 degree
    P_(4,4) = 1; // Angular speed to 180 degree/sec

    const double cos_y = std::cos(z(1));
    const double sin_y = std::sin(z(1));
    const double p_x = z(0) * cos_y;
    const double p_y = z(0) * sin_y;   
    MatrixXd J = MatrixXd(2,3);
    J.fill(0.);
    J(0,0) = cos_y;
    J(1,0) = sin_y;
    J(0,1) = -z(0)*sin_y;
    J(1,1) =  z(0)*cos_y;
    x_(0)  = p_x;
    x_(1)  = p_y;
    P_.topLeftCorner(2,2) = J*R_radar_*J.transpose();
    is_initialized_ = true;
    return;
  }
  
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);  

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // extract values for better readability
    const double& p_x = Xsig_pred_(0,i);
    const double& p_y = Xsig_pred_(1,i);
    const double& v   = Xsig_pred_(2,i);
    const double& yaw = Xsig_pred_(3,i);

    const double v1 = std::cos(yaw)*v;
    const double v2 = std::sin(yaw)*v;

    // measurement model
    double rho = sqrt(p_x*p_x + p_y*p_y);                      
    //if(rho < 1e-3) rho = 1e-3;
    Zsig(0,i) = rho;                       // r 
    Zsig(1,i) = atan2(p_y,p_x);            // phi
    Zsig(2,i) = (p_x*v1 + p_y*v2) / rho;   // r_dot
  }
  
  // calculate mean predicted measurement
  double r_mean, x_mean, y_mean, r_dot_mean;
  r_mean = x_mean = y_mean = r_dot_mean = 0;
  for (int i=0; i < 2*n_aug_+1; ++i) {
    const VectorXd& col = Zsig.col(i);
    r_mean += weights_(i) * col(0);
    x_mean += weights_(i) * std::cos(col(1));
    y_mean += weights_(i) * std::sin(col(1));
    r_dot_mean += weights_(i) * col(2);
  }
  z_pred(0) = r_mean;
  z_pred(1) = atan2(y_mean, x_mean);
  z_pred(2) = r_dot_mean;
  
  // calculate innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 sigma points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;    
    // angle normalization
    normalizeAngle(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  S = S + R_radar_;

  // calculate cross correlation matrix
  Tc.fill(0.);
  for (int i=1; i<2*n_aug_+1; ++i) {  
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    normalizeAngle(z_diff(1));
    normalizeAngle(x_diff(3));
    Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * z_diff.transpose();
  }
  
  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  // update state mean and covariance matrix  
  VectorXd y = z - z_pred;
  normalizeAngle(y(1));
  x_ = x_ + K*y;
  normalizeAngle(x_(3));  
  P_ = P_ - K*S*K.transpose();


  // Compute consistency  
  std::cout << "1, " << y.transpose()*S.inverse()*y << std::endl;
}


void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out, const VectorXd x, const MatrixXd P, const double& lambda){

  //set state dimension
  const int n_x = x.rows();
  
  //create sigma point matrix  
  MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);

  //calculate square root of P
  MatrixXd A = P.llt().matrixL();

  Xsig.col(0) = x;
  const double scale = sqrt(lambda+n_x);
  for(int i=0;i<n_x;i++){
      Xsig.col(i+1)     = x + scale*A.col(i);
      Xsig.col(i+n_x+1) = x - scale*A.col(i);
  }
  
  //write result
  *Xsig_out = Xsig;

}