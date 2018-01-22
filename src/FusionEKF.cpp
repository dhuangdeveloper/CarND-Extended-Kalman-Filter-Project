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

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, 1, 0,
		0, 1, 0, 1,
		0, 0, 1, 0,
		0, 0, 0, 1;
		
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << 1, 0, 1, 0,
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
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */	  
	  cout << "RADAR" << endl;
	  float ro =  measurement_pack.raw_measurements_[0];
	  float theta =  measurement_pack.raw_measurements_[1];
	  float ro_dot =  measurement_pack.raw_measurements_[2];
	  ekf_.x_ << ro * cos(theta), ro * sin(theta), ro_dot * cos(theta), ro_dot * sin(theta);
	  // TO DO P value only a place holder
	  ekf_.P_ = MatrixXd(4, 4);
	  ekf_.P_ << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1000, 0,
		0, 0, 0, 1000;
	  previous_timestamp_ = measurement_pack.timestamp_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	  cout << "Laser" << endl;
	  ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;	  
	  previous_timestamp_ = measurement_pack.timestamp_;	  
	  // TO DO P value only a place holder	  
	  ekf_.P_ = MatrixXd(4, 4);	  
	  ekf_.P_ << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1000, 0,
		0, 0, 0, 1000;	  
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  
  cout << "Initialized" << endl;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
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
  cout << "Predict" << endl;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  cout << "Update" << endl;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	cout << "Radar Update" << endl;
	ekf_.R_ = MatrixXd(3,3);
	ekf_.R_ = R_radar_;
	VectorXd z = VectorXd(3);
	z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[2];
	//ekf_.UpdateEKF(z);
  } else {
	cout << "Laser Update" << endl;
	ekf_.R_ = MatrixXd(2,2);
	ekf_.R_ = R_laser_;	
	ekf_.H_ = MatrixXd(2,4);
	ekf_.H_ << 1, 0, 0, 0,
	  0, 1, 0, 0;
    // Laser updates
	VectorXd z = VectorXd(2);
	z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];		
	ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
