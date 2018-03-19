#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Angle normalization
 * @param {double} a The angle in radians
 * 
 * Returns an angle in the range [-pi, pi] 
 */
double anorm(double a) {
  while(a > M_PI)
    a -= 2.*M_PI;
  while(a <-M_PI)
    a += 2.*M_PI;
  return a;
}


/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // If this is false, laser measurements will be ignored (except during init).
  use_laser_ = true;

  // If this is false, radar measurements will be ignored (except during init).
  use_radar_ = true;

  // initial state vector
  // State is x, y, tangential velocity, yaw, yaw velocity: 
  // px, py, v, psi, psidot
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.fill(0.0);
  int i = 0;
  P_(i, i) = 1; i++; // Var x
  P_(i, i) = 1; i++; // Var y
  P_(i, i) = 2; i++; // Var v
  P_(i, i) = 1; i++; // Var yaw
  P_(i, i) = 2; i++; // Var vyaw

  // process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .1 * 2.5;

  // process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .1 * 2.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  // state dimension
  n_x_ = 5;

  // augmented state dimension
  n_aug_ = 7;

  // number of sigma points
  n_sigma_ = 2 * n_aug_ + 1;

  // (unaugmented) sigma points after motion model
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);

  // sigma point weights
  weights_ = VectorXd(n_sigma_);
  lambda_ = 3 - n_x_; // Sigma point spreading factor; somewhat arbitrarily chosen
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<n_sigma_; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
}


UKF::~UKF() {}


/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /////////////// Initialize from first measurment. //////////////////
  if(! is_initialized_) {
    previous_timestamp_ = 1e-6 * meas_package.timestamp_;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      print("Initializing with RADAR.");
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double ro, theta, x, y, v, psidot;
      ro = meas_package.raw_measurements_[0];
      theta = meas_package.raw_measurements_[1];

      x = ro * cos(theta);
      y = ro * sin(theta);

      // "Although radar gives velocity data in the form of the range rate,
      //  a radar measurement does not contain enough information 
      //  to determine the state variable velocities."
      v = 0;
      psidot = 0;
      x_ << x, y, v, anorm(theta), psidot;

      // Assume RADAR gives us a worse initial position (x, y) guess.
      P_(0, 0) *= 4; // Vx
      P_(1, 1) *= 4; // Vy

    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      print("Initializing with LIDAR.");
      double x = meas_package.raw_measurements_[0];
      double y = meas_package.raw_measurements_[1];
      x_ << x, y, 0, anorm(tan(y / x)), 0;

      // Assume LIDAR gives us a better initial position (x, y) guess.
      P_(0, 0) *= .1; // Vx
      P_(1, 1) *= .1; // Vy
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  /////////////// Do an update with the provided sensor data. //////////
  else {

    // First, do a prediction to the current time
    // (but only if the sensor in question is marked as "use_").
    double t = 1e-6 * meas_package.timestamp_;
    double dt = t - previous_timestamp_;

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
      previous_timestamp_ = t;
      Prediction(dt);
      UpdateRadar(meas_package);
    } else 
    if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
      previous_timestamp_ = t;
      Prediction(dt);
      UpdateLidar(meas_package);
    }
  }
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Generate initial augmented sigma points.
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);
  GenerateAugmentedSigmaPoints(&Xsig_aug);

  // Predict sigma points.
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    // Extract values for better readability.
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p, py_p;

    // Avoid division by zero.
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    // Define noise components.
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // Add noise.
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // Write predicted sigma point into the current column.
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  //////// Evaluate predicted mean and covariance. ////////

  // predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    // angle normalization
    // TODO: Use a state vector object that handles this automatically with a setter method.
    x_diff(3) = anorm(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}


/**
 * Generate augmented sigma points from the current state (includes two noise values).
 * @param {MatrixXd*} Pointer to array in which to write the result.
 */
void UKF::GenerateAugmentedSigmaPoints(MatrixXd* Xsig_out) {

  // augmented mean vector
  VectorXd x_aug = VectorXd(7);

  // augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  // sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
  // Write result.
  *Xsig_out = Xsig_aug;
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage &meas_package) {
  /**
  Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // TODO: Obviously there's a lot of code duplication between this and UpdateRadar; refactor to reuse.

  ///////////////// PREDICT MEASUREMENT ///////////////// 

  // Predict lidar measurment.
  // 1. Take predicted (unaugmented) sigma points, and
  //    transform to measurement space with 0 noise.
  // 2. Evaluate predicted measurement
  //  a) mean and
  //  b) covariance (with added measurement noise covariance).

  // 1
  int n_z = 2;
  // Create matrix for sigma points in measurement space.
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  // Transform sigma points into measurement space.
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    Zsig(0, i) = Xsig_pred_(0,i);
    Zsig(1, i) = Xsig_pred_(1,i);
  }

  // 2a
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // 2b
  // innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    z_diff(1) = anorm(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise covariance matrix.
  MatrixXd R = MatrixXd(n_z, n_z);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;


  /////////////////  APPLY BAYES RULE ///////////////// 

  // Extract the measurment vector.
  VectorXd z = meas_package.raw_measurements_;

  // Create matrix for cross correlation Tc.
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // Calculate cross correlation matrix.
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    z_diff(1) = anorm(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    x_diff(3) = anorm(x_diff(3));

    // weighted sum
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  z_diff(1) = anorm(z_diff(1));

  // Update state mean and covariance matrix.
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}


/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage &meas_package) {
  /**
  Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  ///////////////// PREDICT MEASUREMENT ///////////////// 

  // Predict radar measurment.
  // 1. Take predicted (unaugmented) sigma points, and
  //    transform to measurement space with 0 noise.
  // 2. Evaluate predicted measurement
  //  a) mean and
  //  b) covariance (with added measurement noise covariance).

  // 1
  int n_z = 3;
  // Create matrix for sigma points in measurement space.
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  // Transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1, i) = atan2(p_y,p_x);                                 //phi
    Zsig(2, i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }


  // 2a
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // 2b
  // innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    z_diff(1) = anorm(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise covariance matrix.
  MatrixXd R = MatrixXd(n_z, n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;


  /////////////////  APPLY BAYES RULE ///////////////// 

  // Extract the measurment vector.
  VectorXd z = meas_package.raw_measurements_;

  // Create matrix for cross correlation Tc.
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // Calculate cross correlation matrix.
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    z_diff(1) = anorm(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    x_diff(3) = anorm(x_diff(3));

    // weighted sum
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  z_diff(1) = anorm(z_diff(1));

  // Update state mean and covariance matrix.
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}
