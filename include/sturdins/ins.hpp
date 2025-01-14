/**
 * *ins.hpp*
 *
 * =======  ========================================================================================
 * @file    sturdins/ins.hpp
 * @brief   Inertial navigation Kalman Filter equations.
 * @date    January 2025
 * @author  Daniel Sturdivant <sturdivant20@gmail.com>
 * @ref     1. "Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems", 2nd
 *              Edition, 2013 - Groves
 *          2. "Fundamentals of Inertial Navigation, Satellite-based Positioning and their
 *              Integration" - Noureldin, Karamat, & Georgy
 * =======  ========================================================================================
 */

// TODO: add dynamic IMU biases?

#ifndef STURDINS_INS_HPP
#define STURDINS_INS_HPP

#include "sturdins/strapdown.hpp"

namespace sturdins {

class Ins : public Strapdown {
 public:
  /**
   * *=== Ins ===*
   * @brief constructor
   * @param lat   Initial latitude [rad]
   * @param lon   Initial longitude [rad]
   * @param alt   Initial altitude [m]
   * @param veln  Initial north velocity [m/s]
   * @param vele  Initial east velocity [m/s]
   * @param veld  Initial down velocity [m/s]
   * @param roll  Initial roll angle [rad]
   * @param pitch Initial pitch angle [rad]
   * @param yaw   Initial yaw angle [rad]
   * @param cb    Initial clock bias [m]
   * @param cd    Initial clock drift [m/s]
   */
  Ins();
  Ins(const double lat,
      const double lon,
      const double alt,
      const double veln,
      const double vele,
      const double veld,
      const double roll,
      const double pitch,
      const double yaw,
      const double cb,
      const double cd);

  /**
   * *=== ~Ins ===*
   * @brief Destructor
   */
  ~Ins();

  /**
   * * === SetImuSpec ===
   * @brief Set the noise parameters of the IMU
   * @param Ba     Accelerometer bias instability [mg]
   * @param Ka     Accelerometer rate random walk [(m/s)/(hr*sqrt(hr))]
   * @param Na     Accelerometer random walk [m/s/sqrt(hr)]
   * @param Bg     Gyroscope bias instability [deg/hr]
   * @param Kg     Gyroscope rate random walk [deg/(hr*sqrt(hr))]
   * @param Ng     Gyroscope random walk [deg/sqrt(hr)]
   */
  void SetImuSpec(const double &Ba, const double &Na, const double &Bg, const double &Ng);

  /**
   * * === SetClockSpec ===
   * @brief Set the noise parameters of the Clock
   * @param h0
   * @param h1
   * @param h2
   */
  void SetClockSpec(const double &h0, const double &h1, const double &h2);

  /**
   * *=== SetClock ===*
   * @brief Set the clock states of the INS system
   * @param cb    Clock bias [m]
   * @param cd    Clock drift [m/s]
   */
  void SetClock(const double &cb, const double &cd);

  /**
   * *=== Propagate ===*
   * @brief Propagate the error state matrices
   * @param wb  Measured angular rates (delta thetas) in the body frame [rad/s]
   * @param fb  Measured specific forces (delta velocities) in the body frame [m/s^2]
   * @param dt  Integration time [s]
   */
  void Propagate(const Eigen::Vector3d &wb, const Eigen::Vector3d &fb, const double &dt);

  /**
   * *=== GnssUpdate ===*
   * @brief Correct state with GPS measurements
   * @param sv_pos      Satellite ECEF positions [m]
   * @param sv_vel      Satellite ECEF velocities [m/s]
   * @param psr         Pseudorange measurements [m]
   * @param psrdot      Pseudorange-rate measurements [m/s]
   * @param psr_var     Pseudorange measurement variance [m^2]
   * @param psrdot_var  Pseudorange-rate measurement variance [(m/s)^2]
   */
  void GnssUpdate(
      const Eigen::MatrixXd &sv_pos,
      const Eigen::MatrixXd &sv_vel,
      const Eigen::VectorXd &psr,
      const Eigen::VectorXd &psrdot,
      const Eigen::VectorXd &psr_var,
      const Eigen::VectorXd &psrdot_var);

  /**
   * @brief States not included from Strapdown
   */
  Eigen::Vector3d bg_;  // gyroscope bias estimate
  Eigen::Vector3d ba_;  // accelerometer bias estimate
  double cb_;           // clock bias estimate
  double cd_;           // clock drift estimate

 private:
  /**
   * @brief Kalman Filter Matrices (these have constant size)
   */
  Eigen::VectorXd x_;  // error state vector
  Eigen::MatrixXd P_;  // error state covariance
  Eigen::MatrixXd F_;  // state transition matrix
  Eigen::MatrixXd Q_;  // process covariance matrix
  Eigen::MatrixXd I17_;

  /**
   * @brief IMU allan variance parameters
   */
  double Sra_;   // accelerometer measurement noise PSD
  double Sbad_;  // accelerometer bias variation PSD
  double Srg_;   // gyroscope measurement noise PSD
  double Sbgd_;  // gyroscope bias variation PSD

  /**
   * @brief Clock allan variance parameters
   */
  double h0_;
  double h1_;
  double h2_;
};

}  // namespace sturdins

#endif