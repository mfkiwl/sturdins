/**
 * *kns.cpp*
 *
 * =======  ========================================================================================
 * @file    sturdins/kns.cpp
 * @brief   Kinematic navigation Kalman Filter equations.
 * @date    January 2025
 * @author  Daniel Sturdivant <sturdivant20@gmail.com>
 * @ref     1. "Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems", 2nd
 *              Edition, 2013 - Groves
 * =======  ========================================================================================
 */

#include "sturdins/kns.hpp"

#include <navtools/constants.hpp>

#include "sturdins/least-squares.hpp"

namespace sturdins {

// *=== Kns ===*
Kns::Kns()
    : x_{Eigen::Vector<double, 8>::Zero()},
      P_{Eigen::Matrix<double, 8, 8>::Zero()},
      F_{Eigen::Matrix<double, 8, 8>::Zero()},
      Q_{Eigen::Matrix<double, 8, 8>::Zero()},
      I8_{Eigen::Matrix<double, 8, 8>::Identity()},
      X1ME2_{1.0 - navtools::WGS84_E2<>} {
}
Kns::Kns(
    const double lat,
    const double lon,
    const double alt,
    const double veln,
    const double vele,
    const double veld,
    const double cb,
    const double cd)
    : phi_{lat},
      lam_{lon},
      h_{alt},
      vn_{veln},
      ve_{vele},
      vd_{veld},
      cb_{cb},
      cd_{cd},
      x_{Eigen::Vector<double, 8>::Zero()},
      P_{Eigen::Matrix<double, 8, 8>::Zero()},
      F_{Eigen::Matrix<double, 8, 8>::Identity()},
      Q_{Eigen::Matrix<double, 8, 8>::Zero()},
      I8_{Eigen::Matrix<double, 8, 8>::Identity()} {
}

// *=== SetPosition ===*
void Kns::SetPosition(const double &lat, const double &lon, const double &alt) {
  phi_ = lat;
  lam_ = lon;
  h_ = alt;
}

// *=== SetVelocity ===*
void Kns::SetVelocity(const double &veln, const double &vele, const double &veld) {
  vn_ = veln;
  ve_ = vele;
  vd_ = veld;
}

// *=== SetClock ===*
void Kns::SetClock(const double &cb, const double &cd) {
  cb_ = cb;
  cd_ = cd;
}

// *=== SetClockSpec ===*
void Kns::SetClockSpec(const double &h0, const double &h1, const double &h2) {
  double LS2 = navtools::LIGHT_SPEED<> * navtools::LIGHT_SPEED<>;
  h0_ = LS2 * h0 / 2.0;
  h1_ = LS2 * 2.0 * h1;
  h2_ = LS2 * navtools::PI_SQU<> * h2;
}

// *=== SetProcessNoise ===*
void Kns::SetProcessNoise(const double &Sa) {
  Sa_ = Sa;
  halfSa_ = Sa_ / 2.0;
  thirdSa_ = Sa_ / 3.0;
}

// *=== Propagate ===*
void Kns::Propagate(const double &dt) {
  /**
   * @brief First order F/Phi matrix Groves Ch.9
   * --                      --
   * |  I3   I3*T   Z31   Z31 |
   * |  Z3    I3    Z31   Z31 |
   * | Z13   Z13     1     T  |
   * | Z13   Z13     0     1  |
   * --                      --
   */
  F_(0, 3) = dt;
  F_(1, 4) = dt;
  F_(2, 5) = dt;
  F_(6, 7) = dt;

  /**
   * @brief Process noise matrix Groves Ch.9
   * --                     --
   * | Sa/3*T^3*I3   Sa/2*T^2*I3         Z31           Z31    |
   * | Sa/2*T^2*I3     Sa*T*I3           Z31           Z31    |
   * |     Z13           Z13       Sp*T + Sf/3*T^3   Sf/2*T^2 |
   * |     Z13           Z13           Sf/2*T^2        Sf*T   |
   * --                     --
   */
  double dtsq = dt * dt;
  double dtcb = dtsq * dt;
  Q_(0, 0) = 0.5 * thirdSa_ * dtcb;
  Q_(1, 1) = Q_(0, 0);
  Q_(2, 2) = Q_(0, 0);
  Q_(0, 3) = 0.5 * halfSa_ * dtsq;
  Q_(1, 4) = Q_(0, 3);
  Q_(2, 5) = Q_(0, 3);
  Q_(3, 0) = Q_(0, 3);
  Q_(4, 1) = Q_(0, 3);
  Q_(5, 2) = Q_(0, 3);
  Q_(3, 3) = 0.5 * Sa_ * dt;
  Q_(4, 4) = Q_(3, 3);
  Q_(5, 5) = Q_(3, 3);
  Q_(6, 6) = 0.5 * (h0_ * dt) + (h1_ * dtsq) + (2.0 / 3.0 * h2_ * dtcb);
  Q_(6, 7) = 0.5 * (h1_ * dt) + (h2_ * dtsq);
  Q_(7, 6) = Q_(6, 7);
  Q_(7, 7) = 0.5 * (h0_ / dt) + h1_ + (8.0 / 3.0 * h2_ * dt);

  // Functions of latitude
  sL_ = std::sin(phi_);
  sLsq_ = sL_ * sL_;
  cL_ = std::cos(phi_);

  // radii of curvature
  double t = 1.0 - navtools::WGS84_E2<> * sLsq_;
  double sqt = std::sqrt(t);
  Re_ = navtools::WGS84_R0<> / sqt;
  Rn_ = navtools::WGS84_R0<> * X1ME2_ / (t * t / sqt);
  He_ = Re_ + h_;
  Hn_ = Rn_ + h_;

  // === Kalman Propagation ===
  P_ = F_ * (P_ + Q_) * F_.transpose() + Q_;
  phi_ += vn_ / Hn_ * dt;
  lam_ += ve_ / (cL_ * He_) * dt;
  h_ -= vd_ * dt;
  cb_ += cd_ * dt;
}

// *=== GnssUpdate ===*
void Kns::GnssUpdate(
    const Eigen::MatrixXd &sv_pos,
    const Eigen::MatrixXd &sv_vel,
    const Eigen::VectorXd &psr,
    const Eigen::VectorXd &psrdot,
    const Eigen::VectorXd &psr_var,
    const Eigen::VectorXd &psrdot_var) {
  // Initialize
  const int N = psr.size();
  const int M = 2 * N;
  Eigen::VectorXd dy(M);
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(M, 8);
  for (int i = 0; i < N; i++) {
    H(i, 6) = 1.0;
    H(N + i, 7) = 1.0;
  }
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(M, M);
  R.diagonal() << psr_var, psrdot_var;

  // Functions of current position
  double sLam = std::sin(lam_);
  double cLam = std::cos(lam_);
  sL_ = std::sin(phi_);
  cL_ = std::cos(phi_);
  Eigen::Matrix3d C_l_e{
      {-sL_ * cLam, -sLam, -cL_ * cLam}, {-sL_ * sLam, cLam, -cL_ * sLam}, {cL_, 0.0, -sL_}};
  Eigen::Matrix3d C_e_l = C_l_e.transpose();

  // radii of curvature
  double t = 1.0 - navtools::WGS84_E2<> * sLsq_;
  double sqt = std::sqrt(t);
  Re_ = navtools::WGS84_R0<> / sqt;
  Rn_ = navtools::WGS84_R0<> * X1ME2_ / (t * t / sqt);
  He_ = Re_ + h_;
  Hn_ = Rn_ + h_;

  // Generate observation predictions
  Eigen::Vector3d u, udot;
  Eigen::Vector3d ecef_p{He_ * cL_ * cLam, He_ * cL_ * sLam, (Re_ * X1ME2_ + h_) * sL_};
  Eigen::Vector3d ecef_v{vn_, ve_, vd_};
  ecef_v = C_l_e * ecef_v;
  double pred_psr, pred_psrdot;
  for (int i = 0; i < N; i++) {
    RangeAndRate(
        ecef_p, ecef_v, cb_, cd_, sv_pos.col(i), sv_vel.col(i), u, udot, pred_psr, pred_psrdot);
    u = -C_e_l * u;
    udot = -C_e_l * udot;
    H(i, 0) = u(0);
    H(i, 1) = u(1);
    H(i, 2) = u(2);
    H(N + i, 0) = udot(0);
    H(N + i, 1) = udot(1);
    H(N + i, 2) = udot(2);
    H(i, 3) = u(0);
    H(i, 4) = u(1);
    H(i, 5) = u(2);
    dy(i) = psr(i) - pred_psr;
    dy(N + i) = psrdot(i) - pred_psrdot;
  }

  // TODO: innovation filter

  // === Kalman Update ===
  Eigen::MatrixXd PHt = P_ * H.transpose();
  Eigen::MatrixXd K = PHt * (H * PHt + R).inverse();
  Eigen::MatrixXd L = I8_ - K * H;
  P_ = L * P_ * L.transpose() + K * R + K.transpose();
  x_ += K * dy;

  // Closed loop error corrections
  phi_ -= x_(0) / Hn_;
  lam_ -= x_(1) / (He_ * cL_);
  h_ += x_(2);
  vn_ -= x_(3);
  ve_ -= x_(4);
  vd_ -= x_(5);
  cb_ += x_(6);
  cd_ += x_(7);
  x_.setZero();
}

}  // namespace sturdins