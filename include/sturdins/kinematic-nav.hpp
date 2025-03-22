/**
 * *kinematic-nav.hpp*
 *
 * =======  ========================================================================================
 * @file    sturdins/kinematic-nav.hpp
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

class KinematicNav {
 public:
  /**
   * *=== KinematicNav ===*
   * @brief constructor
   * @param lat   Initial latitude [rad]
   * @param lon   Initial longitude [rad]
   * @param alt   Initial altitude [m]
   * @param veln  Initial north velocity [m/s]
   * @param vele  Initial east velocity [m/s]
   * @param veld  Initial down velocity [m/s]
   */
  KinematicNav();
  KinematicNav(
      const double lat,
      const double lon,
      const double alt,
      const double veln,
      const double vele,
      const double veld,
      const double cb,
      const double cd);
  KinematicNav(
      const double lat,
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
   * *=== SetAttitude ===*
   * @brief Set the attitude of the INS system
   * @param roll   Roll [rad]
   * @param pitch  Pitch [rad]
   * @param yaw    Yaw (heading) [rad]
   */
  void SetAttitude(const double &roll, const double &pitch, const double &yaw);
  void SetAttitude(const Eigen::Ref<const Eigen::Matrix3d> &C);

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
   * @param h0  white frequency modulation
   * @param h1  flicker frequency modulation
   * @param h2  random walk frequency modulation
   */
  void SetClockSpec(const double &h0, const double &h1, const double &h2);

  /**
   * * === SetProcessNoise ===
   * @brief Set the noise parameter of the Kinematic (constant velocity) model
   * @param Svel  PSD of expected acceleration white noise [(m/s^2)^2]
   * @param Satt  PSD of expected angular rate white noise [(rad/s)^2]
   */
  void SetProcessNoise(const double &Svel, const double &Satt);

  /**
   * *=== Propagate ===*
   * @brief Propagate the error state matrices
   * @param dt  Integration time [s]
   */
  void Propagate(const double &dt);
  void FalsePropagateState(
      Eigen::Ref<Eigen::Vector3d> ecef_p,
      Eigen::Ref<Eigen::Vector3d> ecef_v,
      double &cb,
      double &cd,
      const double &dt);

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
      const Eigen::Ref<const Eigen::Matrix3Xd> &sv_pos,
      const Eigen::Ref<const Eigen::Matrix3Xd> &sv_vel,
      const Eigen::Ref<const Eigen::VectorXd> &psr,
      const Eigen::Ref<const Eigen::VectorXd> &psrdot,
      const Eigen::Ref<const Eigen::VectorXd> &psr_var,
      const Eigen::Ref<const Eigen::VectorXd> &psrdot_var);

  /**
   * *=== PhasedArrayUpdate ===*
   * @brief Correct state with GPS measurements
   * @param sv_pos      Satellite ECEF positions [m]
   * @param sv_vel      Satellite ECEF velocities [m/s]
   * @param psr         Pseudorange measurements [m]
   * @param psrdot      Pseudorange-rate measurements [m/s]
   * @param psr_var     Pseudorange measurement variance [m^2]
   * @param psrdot_var  Pseudorange-rate measurement variance [(m/s)^2]
   */
  void PhasedArrayUpdate(
      const Eigen::Ref<const Eigen::Matrix3Xd> &sv_pos,
      const Eigen::Ref<const Eigen::Matrix3Xd> &sv_vel,
      const Eigen::Ref<const Eigen::VectorXd> &psr,
      const Eigen::Ref<const Eigen::VectorXd> &psrdot,
      const Eigen::Ref<const Eigen::MatrixXd> &phase,
      const Eigen::Ref<const Eigen::VectorXd> &psr_var,
      const Eigen::Ref<const Eigen::VectorXd> &psrdot_var,
      const Eigen::Ref<const Eigen::MatrixXd> &phase_var,
      const Eigen::Ref<const Eigen::Matrix3Xd> &ant_xyz,
      const int &n_ant,
      const double &lamb);

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
  Eigen::Vector3d ecef_p_;
  Eigen::Vector3d ecef_v_;
  Eigen::Vector4d q_b_l_;
  Eigen::Matrix3d C_b_l_;
  Eigen::MatrixXd P_;  // error state covariance

 private:
  /**
   * @brief Process noise allan variance parameters
   */
  double h0_;
  double h1_;
  double h2_;
  double Sv_;
  double Sa_;
  double halfSv_;
  double thirdSv_;

  /**
   * @brief Kalman Filter Matrices (these have constant size)
   */
  Eigen::VectorXd x_;  // error state vector
  Eigen::MatrixXd F_;  // state transition matrix
  Eigen::MatrixXd Q_;  // process covariance matrix
  Eigen::MatrixXd I11_;
  bool is_init_;

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

  void KalmanUpdate(
      const Eigen::Ref<const Eigen::MatrixXd> &R,
      const Eigen::Ref<const Eigen::MatrixXd> &H,
      const Eigen::Ref<const Eigen::VectorXd> &dy);
};

}  // namespace sturdins

#endif