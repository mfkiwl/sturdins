/**
 * *kns.hpp*
 *
 * =======  ========================================================================================
 * @file    sturdins/kns.hpp
 * @brief   Kinematic navigation Kalman Filter equations.
 * @date    January 2025
 * @author  Daniel Sturdivant <sturdivant20@gmail.com>
 * @ref     1. "Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems", 2nd
 *              Edition, 2013 - Groves
 * =======  ========================================================================================
 */

// TODO: increase model order to include accelerations (angular rates for versions with attitude)?

#ifndef STURDINS_KNS_HPP
#define STURDINS_KNS_HPP

#include <Eigen/Dense>

namespace sturdins {

class Kns {
 public:
  /**
   * *=== Kns ===*
   * @brief constructor
   * @param lat   Initial latitude [rad]
   * @param lon   Initial longitude [rad]
   * @param alt   Initial altitude [m]
   * @param veln  Initial north velocity [m/s]
   * @param vele  Initial east velocity [m/s]
   * @param veld  Initial down velocity [m/s]
   */
  Kns();
  Kns(const double lat,
      const double lon,
      const double alt,
      const double veln,
      const double vele,
      const double veld,
      const double cb,
      const double cd);

  /**
   * *=== SetPosition ===*
   * @brief Set the position of the INS system
   * @param lat Latitude [rad]
   * @param lon Longitude [rad]
   * @param alt Altitude [m]
   */
  void SetPosition(const double &lat, const double &lon, const double &alt);

  /**
   * *=== SetVelocity ===*
   * @brief Set the velocity of the INS system
   * @param veln   North Velocity [m/s]
   * @param vele   East Velocity [m/s]
   * @param veld   Down Velocity [m/s]
   */
  void SetVelocity(const double &veln, const double &vele, const double &veld);

  /**
   * *=== SetClock ===*
   * @brief Set the clock states of the INS system
   * @param cb    Clock bias [m]
   * @param cd    Clock drift [m/s]
   */
  void SetClock(const double &cb, const double &cd);

  /**
   * * === SetClockSpec ===
   * @brief Set the noise parameters of the Clock
   * @param h0
   * @param h1
   * @param h2
   */
  void SetClockSpec(const double &h0, const double &h1, const double &h2);

  /**
   * * === SetProcessNoise ===
   * @brief Set the noise parameter of the Kinematic (constant velocity) model
   * @param Sa
   */
  void SetProcessNoise(const double &Sa);

  /**
   * *=== Propagate ===*
   * @brief Propagate the error state matrices
   * @param dt  Integration time [s]
   */
  void Propagate(const double &dt);

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

 private:
  /**
   * @brief states
   */
  double phi_;  // Latitude [rad]
  double lam_;  // Longitude [rad]
  double h_;    // Altitude [m]
  double vn_;   // North velocity [m/s]
  double ve_;   // East velocity [m/s]
  double vd_;   // Down velocity [m/s]
  double cb_;   // clock bias [m]
  double cd_;   // clock drift [m/s]

  /**
   * @brief Process noise allan variance parameters
   */
  double h0_;
  double h1_;
  double h2_;
  double Sa_;
  double halfSa_;
  double thirdSa_;

  /**
   * @brief Kalman Filter Matrices (these have constant size)
   */
  Eigen::VectorXd x_;  // error state vector
  Eigen::MatrixXd P_;  // error state covariance
  Eigen::MatrixXd F_;  // state transition matrix
  Eigen::MatrixXd Q_;  // process covariance matrix
  Eigen::MatrixXd I8_;

  /**
   * @brief Functions of latitude
   */
  double sL_;    // sin(phi)
  double cL_;    // cos(phi)
  double sLsq_;  // sin(phi)^2

  /**
   * @brief Radii of curvature
   */
  double X1ME2_;  // 1 - e^2
  double Re_;     // Meridian radius
  double Rn_;     // Transverse radius
  double Hn_;     // Rn + h
  double He_;     // Re + h
};

}  // namespace sturdins

#endif