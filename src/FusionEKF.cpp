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
  Hj_ = MatrixXd(3, 4); //related to delta_t, initialized later
  ekf_.Q_ = MatrixXd(4, 4);//related to delta_t, initialized later
  ekf_.F_ = MatrixXd(4, 4);//related to delta_t, initialized later

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
  H_laser_<< 1, 0, 0, 0,
			  0, 1, 0, 0;
			  
	// initialize the initial measurement
 /* 
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
	
	//intialize the initial state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;
  */
 

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
	
	//intialize the  state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	  double x_co=measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]);
	  double y_co=measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]);

        ekf_.x_ << x_co, y_co, 0, 0;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	  	ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

		
    }
	//record the previous_timestamp_
	previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

 

	
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
    //caltulate delta_t
	double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds


	if(fabs(dt) < 0.0001){
		cout << " FusionEKF - dt equals to Zero" << endl;
		dt=0.0001;
	}

	previous_timestamp_ = measurement_pack.timestamp_;
	
	double dt_2 = dt * dt;
	double dt_3 = dt_2 * dt;
	double dt_4 = dt_3 * dt;

	//Modify the F matrix so that the time is integrated
		//the initial transition matrix F_
	
	ekf_.F_ << 1, 0, dt, 0,
			  0, 1, 0, dt,
			  0, 0, 1, 0,
			  0, 0, 0, 1;
	//ekf_.F_(0, 2) = dt;
	//ekf_.F_(1, 3) = dt;
	
	//set the process acceleration noise components
	double noise_ax = 9.0;
	double noise_ay = 9.0;

	//set the process covariance matrix Q
	
	ekf_.Q_ <<  dt_4/4*noise_ax, 0.0, dt_3/2*noise_ax, 0.0,
			   0.0, dt_4/4*noise_ay, 0.0, dt_3/2*noise_ay,
			   dt_3/2*noise_ax, 0.0, dt_2*noise_ax, 0.0,
			   0.0, dt_3/2*noise_ay, 0.0, dt_2*noise_ay;
/*	
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	ekf_.Init(ekf_.x_ , ekf_.P_ , ekf_.F_ , Hj_, R_radar_, ekf_.Q_ );
	} 
	else {
    ekf_.Init(ekf_.x_ , ekf_.P_ , ekf_.F_ , H_laser_, R_laser_ , ekf_.Q_ );
	}
*/
    //cout << "before predict() x_ = " << ekf_.x_ << endl;
	ekf_.Predict();
    //cout << "after predict() x_ = " << ekf_.x_ << endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      double x_co=measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]);
      double y_co=measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]);
      if((fabs(x_co) < 0.0001)&& (fabs(y_co) < 0.0001)){
          cout << " Discard the radar measurement update when the radar measurement shows x=y=0" << endl;
          return;
      }
	Hj_=tools.CalculateJacobian(ekf_.x_);
	ekf_.Init(ekf_.x_ , ekf_.P_ , ekf_.F_ , Hj_, R_radar_, ekf_.Q_ );
	//ekf_.R_=R_radar_;
	ekf_.UpdateEKF(measurement_pack.raw_measurements_,Hj_);
  } else {
    // Laser updates
	  //cout << "x_ = " << ekf_.x_ << endl;
	  //cout << "H_ = " << ekf_.H_ << endl;
    ekf_.Init(ekf_.x_ , ekf_.P_ , ekf_.F_ , H_laser_, R_laser_ , ekf_.Q_ );
    //ekf_.R_=R_laser_;
	ekf_.Update(measurement_pack.raw_measurements_); 
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
