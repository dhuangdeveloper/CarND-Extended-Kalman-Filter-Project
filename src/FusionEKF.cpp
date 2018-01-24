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
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, 1, 0,
		0, 1, 0, 1,
		0, 0, 1, 0,
		0, 0, 0, 1;		

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
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 0, 0;
    ekf_.P_ = MatrixXd(4, 4);	  
    ekf_.P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;	
	ekf_.H_ = MatrixXd(2,4);
	ekf_.H_ << 1, 0, 0, 0,
	  0, 1, 0, 0;	  
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */	  	  
	  float ro =  measurement_pack.raw_measurements_[0];
	  float theta =  measurement_pack.raw_measurements_[1];
	  float ro_dot =  measurement_pack.raw_measurements_[2];	  
	  
	  ekf_.x_ << ro * cos(theta), ro * sin(theta), 0, 0;
	  
	  // initialize P based on R_radar_ info
	  MatrixXd Hj_pxpy_inv(2,3);
	  Hj_pxpy_inv << cos(theta), -ro * sin(theta), 0,
		sin(theta), ro * cos(theta),0;
	  MatrixXd R_pxy_radar;
	  R_pxy_radar = Hj_pxpy_inv * R_radar_ * Hj_pxpy_inv.transpose();
      ekf_.P_(0, 0) = R_pxy_radar(0, 0);
	  ekf_.P_(0, 1) = R_pxy_radar(0, 1);
	  ekf_.P_(1, 0) = R_pxy_radar(1, 0);
	  ekf_.P_(1, 1) = R_pxy_radar(1, 1);	
	  
	  previous_timestamp_ = measurement_pack.timestamp_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */	  
	  
	  ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;	  
	  
      ekf_.P_(0, 0) = R_laser_(0, 0);
	  ekf_.P_(0, 1) = R_laser_(0, 1);
	  ekf_.P_(1, 0) = R_laser_(1, 0);
	  ekf_.P_(1, 1) = R_laser_(1, 1);	
	  previous_timestamp_ = measurement_pack.timestamp_;	  
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }  

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  // Code is from what I submitted in quiz, with some portion coming from prefilled-code in the quiz
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //set the process covariance matrix Q
  float noise_ax = 9;
  float noise_ay = 9;  
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0, 
    0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
    dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
    0, dt_3/2*noise_ay, 0, dt_2*noise_ay;  
  ekf_.Predict();


  /*****************************************************************************
   *  Update
   ****************************************************************************/
 
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates	
	ekf_.R_ = MatrixXd(3,3);
	ekf_.R_ = R_radar_;
	VectorXd z = VectorXd(3);
	z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[2];
	ekf_.UpdateEKF(z);
  } else {	
	ekf_.R_ = MatrixXd(2,2);
	ekf_.R_ = R_laser_;	
    // Laser updates
	VectorXd z = VectorXd(2);
	z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];		
	ekf_.Update(z);
  }

  // print the output  
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
