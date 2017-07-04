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

  H_laser_ << 1,0,0,0,
              0,1,0,0;


  //acceleration noise component
  ekf_.F_= MatrixXd(4, 4);

  ekf_.P_ = MatrixXd(4, 4);

  ekf_.F_ << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;
 ekf_.P_ << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1000, 0,
         0, 0, 0, 1000;


  ekf_.Q_ = MatrixXd(4, 4);

  //set the acceleration noise components
  noise_ax = 9 ;
  noise_ay = 9;


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
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;

    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 0.2, 0.5;  //play around for RMSE



    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      cout << "RADAR Data: " << endl;

      double rho = measurement_pack.raw_measurements_[0]; // range
      double phi = measurement_pack.raw_measurements_[1]; // bearing
      double rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho

      // Coordinates convertion from polar to cartesian
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      ekf_.x_ << x, y, vx , vy;
      cout << "END RADAR Data: " << endl;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      cout << "LASER Data: " << endl;

        ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;





      cout << "END LASER Data: " << endl;


      //

    }

    previous_timestamp_ = measurement_pack.timestamp_;

//    cout << "Time: " + previous_timestamp_ << endl;

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

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

  cout << "DT" << endl;
  previous_timestamp_ = measurement_pack.timestamp_;
  cout << "Previous timestamp" << endl;

  float dt_2 = dt*dt;
  float dt_3 = dt_2*dt;
  float dt_4 = dt_3*dt;

  //modify the F matrix so that the time is integrated
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  //set the process covariance matrix Q
  ekf_.Q_ <<  dt_4/4.0 * noise_ax, 0, dt_3 / 2.0 * noise_ax, 0,
          0, dt_4 / 4.0  * noise_ay, 0, dt_3 / 2.0 * noise_ay,
          dt_3 / 2.0 * noise_ax, 0, dt_2 * noise_ax, 0,
          0, dt_3 / 2.0 * noise_ay, 0, dt_2 * noise_ay;
  if(dt > 0.001) {
    ekf_.Predict();
  }
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
    cout << "RADAR update: " << endl;

    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {

    // Laser updates
    cout << "LASER update: " << endl;

    ekf_.H_= H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
    cout << "END LASER update: " << endl;

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
