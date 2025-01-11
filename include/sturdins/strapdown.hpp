/**
 * *strapdown.hpp*
 *
 * =======  ========================================================================================
 * @file    sturdins/strapdown.hpp
 * @brief   Inertial navigation strapdown integration equations.
 * @date    January 2025
 * @author  Daniel Sturdivant <sturdivant20@gmail.com>
 * @ref     1. "Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems", 2nd
 *              Edition, 2013 - Groves
 *          2. "Attitude and Position Integration" - VectorNav
 *              https://www.vectornav.com/resources/inertial-navigation-primer/
 * =======  ========================================================================================
 */

#ifndef STURDINS_STRAPDOWN_HPP
#define STURDINS_STRAPDOWN_HPP

#include <Eigen/Dense>
#include <navtools/constants.hpp>

namespace sturdins {

class Strapdown {
 public:
  /**
   * *=== Strapdown ===*
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
   */
  Strapdown();
  Strapdown(
      const double lat,
      const double lon,
      const double alt,
      const double veln,
      const double vele,
      const double veld,
      const double roll,
      const double pitch,
      const double yaw);

  /**
   * *=== ~Strapdown ===*
   * @brief Destructor
   */
  ~Strapdown();

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
   * @brief Set the orientation of the INS system
   * @param roll   Roll angle [rad]
   * @param pitch  Pitch angle [rad]
   * @param yaw    Yaw angle [rad]
   */
  void SetAttitude(const double &roll, const double &pitch, const double &yaw);
  void SetAttitude(const Eigen::Matrix3d &C);

  /**
   * @brief Integrate measured angular rates (delta thetas) and specific forces (delta velocities)
   * @param wb  Measured angular rates (delta thetas) in the body frame [rad/s]
   * @param fb  Measured specific forces (delta velocities) in the body frame [m/s^2]
   * @param dt  Integration time [s]
   */
  void Mechanize(const Eigen::Vector3d &wb, const Eigen::Vector3d &fb, const double &dt);

 protected:
  /**
   * @brief states
   */
  double phi_;             // Latitude [rad]
  double lam_;             // Longitude [rad]
  double h_;               // Altitude [m]
  double vn_;              // North velocity [m/s]
  double ve_;              // East velocity [m/s]
  double vd_;              // Down velocity [m/s]
  Eigen::Vector4d q_b_l_;  // Body -> Local/Nav frame quaternion
  Eigen::Matrix3d C_b_l_;  // Body -> Local/Nav frame DCM

  /**
   * @brief Functions of latitude
   */
  double sL_;    // sin(phi)
  double cL_;    // cos(phi)
  double tL_;    // tan(phi)
  double sLsq_;  // sin(phi)^2
  double cLsq_;  // cos(phi)^2

  /**
   * @brief Radii of curvature
   */
  double X1ME2_;    // 1 - e^2
  double X1ME2sq_;  // (1 - e^2)^2
  double Re_;       // Meridian radius
  double Rn_;       // Transverse radius
  double Rg_;       // Geocentric radius
  double Hn_;       // Rn + h
  double He_;       // Re + h
  double Hnsq_;     // (Rn + h)^2
  double Hesq_;     // (Re + h)^2

  /**
   * @brief Coriolis and gravity
   */
  double g0_;               // Somigliana model gravity
  Eigen::Vector3d g_;       // Gravity vector in NED frame
  Eigen::Vector3d w_en_n_;  // Rotation rate of the ECEF frame in the NED frame
  Eigen::Vector3d w_ie_n_;  // Rotation rate of earth in the NED frame
  double wR0sqRpMu_;        // (w_ie * R0)^2 * Rp / Mu

  /**
   * @brief Calculate gravity vector based on current LLA
   */
  void GravityVector();

  /**
   * @brief Calculate earth rotation rate based on current LLA
   */
  void EarthRateVector();

  /**
   * @brief Calculate transport rate of ECEF frame based on current LLA
   */
  void TransportRateVector();

  /**
   * @brief Calculate radii of curvature
   */
  template <bool CalcGeocentricRadius = false>
  void RadiiOfCurvature() {
    double t = 1.0 - navtools::WGS84_E2<> * sLsq_;
    double sqt = std::sqrt(t);
    Re_ = navtools::WGS84_R0<> / sqt;
    Rn_ = navtools::WGS84_R0<> * X1ME2_ / (t * t / sqt);
    if constexpr (CalcGeocentricRadius) {
      Rg_ = Re_ * std::sqrt(cLsq_ + sLsq_ * X1ME2sq_);
    }
  }
};

}  // namespace sturdins

#endif