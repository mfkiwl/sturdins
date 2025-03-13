/**
 * *strapdown.cpp*
 *
 * =======  ========================================================================================
 * @file    sturdins/strapdown.cpp
 * @brief   Inertial navigation strapdown integration equations.
 * @date    January 2025
 * @author  Daniel Sturdivant <sturdivant20@gmail.com>
 * @ref     1. "Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems", 2nd
 *              Edition, 2013 - Groves
 * =======  ========================================================================================
 */

#include "sturdins/strapdown.hpp"

#include <navtools/attitude.hpp>

namespace sturdins {

// *=== Strapdown ===*
Strapdown::Strapdown()
    : X1ME2_{1.0 - navtools::WGS84_E2<>},
      g_{Eigen::Vector3d::Zero()},
      w_en_n_{Eigen::Vector3d::Zero()},
      w_ie_n_{Eigen::Vector3d::Zero()},
      wR0sqRpMu_{
          navtools::WGS84_OMEGA<> * navtools::WGS84_OMEGA<> * navtools::WGS84_R0<> *
          navtools::WGS84_R0<> * navtools::WGS84_RP<> / navtools::WGS84_MU<>} {
}
Strapdown::Strapdown(
    const double lat,
    const double lon,
    const double alt,
    const double veln,
    const double vele,
    const double veld,
    const double roll,
    const double pitch,
    const double yaw)
    : phi_{lat},
      lam_{lon},
      h_{alt},
      vn_{veln},
      ve_{vele},
      vd_{veld},
      X1ME2_{1.0 - navtools::WGS84_E2<>},
      g_{Eigen::Vector3d::Zero()},
      w_en_n_{Eigen::Vector3d::Zero()},
      w_ie_n_{Eigen::Vector3d::Zero()},
      wR0sqRpMu_{
          navtools::WGS84_OMEGA<> * navtools::WGS84_OMEGA<> * navtools::WGS84_R0<> *
          navtools::WGS84_R0<> * navtools::WGS84_RP<> / navtools::WGS84_MU<>} {
  // initialize attitude
  Eigen::Vector3d rpy{roll, pitch, yaw};
  navtools::euler2quat<double>(q_b_l_, rpy, true);
  navtools::quat2dcm<double>(C_b_l_, q_b_l_);
}

// *=== ~Strapdown ===*
Strapdown::~Strapdown() {
}

// *=== SetPosition ===*
void Strapdown::SetPosition(const double &lat, const double &lon, const double &alt) {
  phi_ = lat;
  lam_ = lon;
  h_ = alt;
}

// *=== SetVelocity ===*
void Strapdown::SetVelocity(const double &veln, const double &vele, const double &veld) {
  vn_ = veln;
  ve_ = vele;
  vd_ = veld;
}

// *=== SetAttitude ===*
void Strapdown::SetAttitude(const double &roll, const double &pitch, const double &yaw) {
  Eigen::Vector3d rpy{roll, pitch, yaw};
  navtools::euler2quat<double>(q_b_l_, rpy, true);
  navtools::quat2dcm<double>(C_b_l_, q_b_l_);
}
void Strapdown::SetAttitude(const Eigen::Matrix3d &C) {
  C_b_l_ = C;
  navtools::dcm2quat<double>(q_b_l_, C_b_l_);
}

// *=== Mechanize ===*
void Strapdown::Mechanize(const Eigen::Vector3d &wb, const Eigen::Vector3d &fb, const double &dt) {
  // Sine functions of latitude
  sL_ = std::sin(phi_);
  cL_ = std::cos(phi_);
  tL_ = sL_ / cL_;
  sLsq_ = sL_ * sL_;

  // Radii of curvature
  double t = 1.0 - navtools::WGS84_E2<> * sLsq_;
  double sqt = std::sqrt(t);
  Re_ = navtools::WGS84_R0<> / sqt;
  Rn_ = navtools::WGS84_R0<> * X1ME2_ / (t * t / sqt);
  He_ = Re_ + h_;
  Hn_ = Rn_ + h_;

  // Gravity and coriolis
  g0_ = 9.7803253359 * ((1.0 + 0.001931853 * sLsq_) / sqt);
  GravityVector();
  EarthRateVector();
  TransportRateVector();

  // --- Attitude Integration ---
  Eigen::Vector3d dTheta = (wb - C_b_l_.transpose() * (w_ie_n_ + w_en_n_)) * dt;
  double gamma = 0.5 * dTheta.norm();
  double cgamma = std::cos(gamma);
  if (gamma < 1e-5) {
    dTheta *= 0.5;
  } else {
    dTheta *= (0.5 * std::sin(gamma) / gamma);
  }
  Eigen::Matrix4d qdot{
      // clang-format off
      {   cgamma, -dTheta(0), -dTheta(1), -dTheta(2)},
      {dTheta(0),     cgamma,  dTheta(2), -dTheta(1)},
      {dTheta(1), -dTheta(2),     cgamma,  dTheta(0)},
      {dTheta(2),  dTheta(1), -dTheta(0),     cgamma}
      // clang-format on
  };
  q_b_l_ = qdot * q_b_l_;
  q_b_l_ /= q_b_l_.norm();

  // --- Velocity Integration ---
  Eigen::Vector3d fn = C_b_l_ * fb;
  Eigen::Vector3d wn = w_ie_n_ + 2.0 * w_en_n_;
  double dv_n = (fn(0) + g_(0) - (wn(1) * vd_ - wn(2) * ve_)) * dt;
  double dv_e = (fn(1) + g_(1) - (wn(0) * vd_ - wn(2) * vn_)) * dt;
  double dv_d = (fn(2) + g_(2) - (wn(0) * ve_ - wn(1) * vn_)) * dt;

  // --- Position Integration ---
  phi_ += (vn_ + 0.5 * dv_n) / Hn_ * dt;
  lam_ += (ve_ + 0.5 * dv_e) / (He_ * cL_) * dt;
  h_ -= (vd_ + 0.5 * dv_d) * dt;

  // Save integration result
  navtools::quat2dcm<double>(C_b_l_, q_b_l_);
  vn_ += dv_n;
  ve_ += dv_e;
  vd_ += dv_d;
}

// *=== GravityVector ===*
void Strapdown::GravityVector() {
  double hR0sq = h_ / navtools::WGS84_R0<>;
  hR0sq *= hR0sq;
  g_(0) = -8.08e-9 * h_ * std::sin(2.0 * phi_);
  // g_(1) = 0.0;
  g_(2) = g0_ * (1.0 -
                 (2.0 * h_ / navtools::WGS84_R0<>)*(
                     1.0 + navtools::WGS84_F<> * (1.0 - 2.0 * sLsq_) + wR0sqRpMu_) +
                 3.0 * hR0sq);
}

// *=== EarthRateVector ===*
void Strapdown::EarthRateVector() {
  w_ie_n_(0) = navtools::WGS84_OMEGA<> * cL_;
  // w_ie_n_(1) = 0.0;
  w_ie_n_(2) = navtools::WGS84_OMEGA<> * sL_;
}

// *=== TransportRateVector ===*
void Strapdown::TransportRateVector() {
  w_en_n_(0) = ve_ / He_;
  w_en_n_(1) = -vn_ / Hn_;
  w_en_n_(2) = -w_en_n_(0) * tL_;
}

}  // namespace sturdins