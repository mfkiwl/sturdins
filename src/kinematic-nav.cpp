/**
 * *kinematic-nav.cpp*
 *
 * =======  ========================================================================================
 * @file    sturdins/kinematic-nav.cpp
 * @brief   Kinematic navigation Kalman Filter equations.
 * @date    January 2025
 * @author  Daniel Sturdivant <sturdivant20@gmail.com>
 * @ref     1. "Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems", 2nd
 *              Edition, 2013 - Groves
 * =======  ========================================================================================
 */

#include "sturdins/kinematic-nav.hpp"

#include <iostream>
#include <navtools/attitude.hpp>
#include <navtools/constants.hpp>
#include <navtools/math.hpp>

#include "sturdins/least-squares.hpp"

namespace sturdins {

// *=== KinematicNav ===*
KinematicNav::KinematicNav()
    : q_b_l_{Eigen::Vector4d{1.0, 0.0, 0.0, 0.0}},
      C_b_l_{Eigen::Matrix3d::Identity()},
      P_{Eigen::Matrix<double, 11, 11>::Zero()},
      x_{Eigen::Vector<double, 11>::Zero()},
      F_{Eigen::Matrix<double, 11, 11>::Identity()},
      Q_{Eigen::Matrix<double, 11, 11>::Zero()},
      I11_{Eigen::Matrix<double, 11, 11>::Identity()},
      is_init_{false},
      X1ME2_{1.0 - navtools::WGS84_E2<>},
      LS2_{navtools::LIGHT_SPEED<> * navtools::LIGHT_SPEED<>} {
  P_.diagonal() << 9.0, 9.0, 9.0, 0.05, 0.05, 0.05, 0.01, 0.01, 0.01, 3.0, 0.1;
}
KinematicNav::KinematicNav(
    const double lat,
    const double lon,
    const double alt,
    const double veln,
    const double vele,
    const double veld,
    const double cb,
    const double cd)
    : KinematicNav() {
  phi_ = lat;
  lam_ = lon;
  h_ = alt;
  vn_ = veln;
  ve_ = vele;
  vd_ = veld;
  cb_ = cb;
  cd_ = cd;
  P_.diagonal() << 9.0, 9.0, 9.0, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 3.0, 0.1;
}
KinematicNav::KinematicNav(
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
    const double cd)
    : KinematicNav(lat, lon, alt, veln, vele, veld, cb, cd) {
  SetAttitude(roll, pitch, yaw);
}

// *=== SetPosition ===*
void KinematicNav::SetPosition(const double &lat, const double &lon, const double &alt) {
  phi_ = lat;
  lam_ = lon;
  h_ = alt;
}

// *=== SetVelocity ===*
void KinematicNav::SetVelocity(const double &veln, const double &vele, const double &veld) {
  vn_ = veln;
  ve_ = vele;
  vd_ = veld;
}

// *=== SetAttitude ===*
void KinematicNav::SetAttitude(const double &roll, const double &pitch, const double &yaw) {
  Eigen::Vector3d euler{roll, pitch, yaw};
  navtools::euler2dcm<double>(C_b_l_, euler, true);
  navtools::euler2quat<double>(q_b_l_, euler, true);
}
void KinematicNav::SetAttitude(const Eigen::Ref<const Eigen::Matrix3d> &C) {
  C_b_l_ = C;
  navtools::dcm2quat<double>(q_b_l_, C_b_l_);
}

// *=== SetClock ===*
void KinematicNav::SetClock(const double &cb, const double &cd) {
  cb_ = cb;
  cd_ = cd;
}

// *=== SetClockSpec ===*
void KinematicNav::SetClockSpec(const double &h0, const double &h1, const double &h2) {
  Sb_ = 1.1 * (h0 / 2.0);
  Sd_ = 1.1 * (h2 * 2.0 * navtools::PI_SQU<>);
  Sbd_ = 1.1 * (h1 * 2.0);
}

// *=== SetProcessNoise ===*
void KinematicNav::SetProcessNoise(const double &Svel, const double &Satt) {
  Sv_ = Svel;
  halfSv_ = Sv_ / 2.0;
  thirdSv_ = Sv_ / 3.0;
  Sa_ = Satt;
}

// *=== Propagate ===*
void KinematicNav::Propagate(const double &dt) {
  /**
   * @brief First order F/Phi matrix Groves Ch.9
   * --                      --
   * |  I3   I3*T    Z3   Z31   Z31 |
   * |  Z3    I3     Z3   Z31   Z31 |
   * |  Z3    Z3     I3   Z31   Z31 |
   * | Z13   Z13    Z13    1     T  |
   * | Z13   Z13    Z13    0     1  |
   * --                      --
   */
  F_(0, 3) = dt;
  F_(1, 4) = dt;
  F_(2, 5) = dt;
  F_(9, 10) = dt;

  /**
   * @brief Process noise matrix Groves Ch.9
   * --                                                                --
   * | Sv/3*T^3*I3   Sv/2*T^2*I3      Z3           Z31           Z31    |
   * | Sv/2*T^2*I3     Sv*T*I3        Z3           Z31           Z31    |
   * |      Z3            Z3       Sa*T*I3         Z31           Z31    |
   * |     Z13           Z13         Z13     Sp*T + Sf/3*T^3   Sf/2*T^2 |
   * |     Z13           Z13         Z13         Sf/2*T^2        Sf*T   |
   * --                                                                --
   */
  double dtsq = dt * dt;
  double dtcb = dtsq * dt;
  Q_(0, 0) = thirdSv_ * dtcb;
  Q_(1, 1) = Q_(0, 0);
  Q_(2, 2) = Q_(0, 0);
  Q_(0, 3) = halfSv_ * dtsq;
  Q_(1, 4) = Q_(0, 3);
  Q_(2, 5) = Q_(0, 3);
  Q_(3, 0) = Q_(0, 3);
  Q_(4, 1) = Q_(0, 3);
  Q_(5, 2) = Q_(0, 3);
  Q_(3, 3) = Sv_ * dt;
  Q_(4, 4) = Q_(3, 3);
  Q_(5, 5) = Q_(3, 3);

  Q_(6, 6) = Sa_ * dt;
  Q_(7, 7) = Q_(6, 6);
  Q_(8, 8) = Q_(6, 6);  // * 16.0;

  Q_(9, 9) = LS2_ * ((Sb_ * dt) + (Sbd_ * dtsq) + (Sd_ * dtcb / 3.0));
  Q_(9, 10) = LS2_ * ((Sbd_ * dt) + (Sd_ * dtsq / 2.0));
  Q_(10, 9) = Q_(9, 10);
  Q_(10, 10) = LS2_ * (((Sb_ / dt) + (Sbd_) + (4.0 / 3.0 * Sd_ * dt)));
  // Q_(9, 9) = h0_ * dt + h2_ * dtcb / 3.0;
  // Q_(9, 10) = h2_ * dtsq / 2.0;
  // Q_(10, 9) = Q_(9, 10);
  // Q_(10, 10) = h2_ * dt;

  // std::cout << "F = \n" << F_ << "\n";
  // std::cout << "Q = \n" << Q_ << "\n";

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
  P_ = F_ * P_ * F_.transpose() + Q_;
  phi_ += vn_ / Hn_ * dt;
  lam_ += ve_ / (cL_ * He_) * dt;
  h_ -= vd_ * dt;
  cb_ += cd_ * dt;
}

void KinematicNav::FalsePropagateState(
    Eigen::Ref<Eigen::Vector3d> ecef_p,
    Eigen::Ref<Eigen::Vector3d> ecef_v,
    double &cb,
    double &cd,
    const double &dt) {
  // update ecef position
  double sLam = std::sin(lam_);
  double cLam = std::cos(lam_);
  sL_ = std::sin(phi_);
  cL_ = std::cos(phi_);
  sLsq_ = sL_ * sL_;
  Eigen::Matrix3d C_l_e{
      {-sL_ * cLam, -sLam, -cL_ * cLam}, {-sL_ * sLam, cLam, -cL_ * sLam}, {cL_, 0.0, -sL_}};
  double t = 1.0 - navtools::WGS84_E2<> * sLsq_;
  double sqt = std::sqrt(t);
  Re_ = navtools::WGS84_R0<> / sqt;
  // Rn_ = navtools::WGS84_R0<> * X1ME2_ / (t * t / sqt);
  He_ = Re_ + h_;
  // Hn_ = Rn_ + h_;
  ecef_p_ << He_ * cL_ * cLam, He_ * cL_ * sLam, (Re_ * X1ME2_ + h_) * sL_;
  ecef_v_ << vn_, ve_, vd_;
  ecef_v_ = C_l_e * ecef_v_;

  // fake-propagate state into provided vectors
  ecef_p = ecef_p_ + ecef_v_ * dt;
  ecef_v = ecef_v_;
  cb = cb_ + cd_ * dt;
  cd = cd_;
}

// *=== GnssUpdate ===*
void KinematicNav::GnssUpdate(
    const Eigen::Ref<const Eigen::Matrix3Xd> &sv_pos,
    const Eigen::Ref<const Eigen::Matrix3Xd> &sv_vel,
    const Eigen::Ref<const Eigen::VectorXd> &psr,
    const Eigen::Ref<const Eigen::VectorXd> &psrdot,
    const Eigen::Ref<const Eigen::VectorXd> &psr_var,
    const Eigen::Ref<const Eigen::VectorXd> &psrdot_var) {
  // Initialize
  const int N = psr.size();
  const int M = 2 * N;
  Eigen::VectorXd dy(M);
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(M, 11);
  for (int i = 0; i < N; i++) {
    H(i, 9) = 1.0;
    H(N + i, 10) = 1.0;
  }
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(M, M);
  R.diagonal() << psr_var, psrdot_var;

  // Functions of current position
  double sLam = std::sin(lam_);
  double cLam = std::cos(lam_);
  sL_ = std::sin(phi_);
  cL_ = std::cos(phi_);
  sLsq_ = sL_ * sL_;
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
  ecef_p_ << He_ * cL_ * cLam, He_ * cL_ * sLam, (Re_ * X1ME2_ + h_) * sL_;
  ecef_v_ << vn_, ve_, vd_;
  ecef_v_ = C_l_e * ecef_v_;
  double pred_psr, pred_psrdot;
  for (int i = 0; i < N; i++) {
    RangeAndRate(
        ecef_p_, ecef_v_, cb_, cd_, sv_pos.col(i), sv_vel.col(i), u, udot, pred_psr, pred_psrdot);
    u = C_e_l * u;
    udot = C_e_l * udot;
    H(i, 0) = u(0);
    H(i, 1) = u(1);
    H(i, 2) = u(2);
    H(N + i, 0) = udot(0);
    H(N + i, 1) = udot(1);
    H(N + i, 2) = udot(2);
    H(N + i, 3) = u(0);
    H(N + i, 4) = u(1);
    H(N + i, 5) = u(2);
    dy(i) = psr(i) - pred_psr;
    dy(N + i) = psrdot(i) - pred_psrdot;
  }

  // Kalman Update
  KalmanUpdate(R, H, dy);
}

// *=== PhasedArrayUpdate ===*
void KinematicNav::PhasedArrayUpdate(
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
    const double &lamb) {
  // Initialize
  const int N = psr.size();
  const int M = 2 * N;
  const int MM = M + (n_ant - 1) * N;
  Eigen::VectorXd dy(MM);
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(MM, 11);
  for (int ii = 0; ii < N; ii++) {
    H(ii, 9) = 1.0;
    H(N + ii, 10) = 1.0;
  }
  // std::cout << "psr_var = " << psr_var(0) << ", psrdot_var = " << psrdot_var(0) << "\n";

  // measurement variance
  int k1, k2;
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(MM, MM);
  for (int ii = 0; ii < N; ii++) {
    k1 = N + ii;
    R(ii, ii) = psr_var(ii);
    R(k1, k1) = psrdot_var(ii);
    for (int jj = 1; jj < n_ant; jj++) {
      k2 = M + (n_ant - 1) * ii + jj - 1;
      R(k2, k2) = phase_var(jj, ii);
    }
  }

  // Functions of current position
  double sLam = std::sin(lam_);
  double cLam = std::cos(lam_);
  sL_ = std::sin(phi_);
  cL_ = std::cos(phi_);
  sLsq_ = sL_ * sL_;
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
  Eigen::Vector3d u, udot, hp;
  ecef_p_ << He_ * cL_ * cLam, He_ * cL_ * sLam, (Re_ * X1ME2_ + h_) * sL_;
  ecef_v_ << vn_, ve_, vd_;
  ecef_v_ = C_l_e * ecef_v_;
  double pred_psr, pred_psrdot, pred_phase;
  Eigen::Matrix3Xd ant_ned = C_b_l_ * ant_xyz;
  // std::cout << "C_b_l = \n" << C_b_l_ << "\n";
  for (int ii = 0; ii < N; ii++) {
    RangeAndRate(
        ecef_p_, ecef_v_, cb_, cd_, sv_pos.col(ii), sv_vel.col(ii), u, udot, pred_psr, pred_psrdot);
    u = C_e_l * u;
    udot = C_e_l * udot;
    H(ii, 0) = u(0);
    H(ii, 1) = u(1);
    H(ii, 2) = u(2);
    H(N + ii, 0) = udot(0);
    H(N + ii, 1) = udot(1);
    H(N + ii, 2) = udot(2);
    H(N + ii, 3) = u(0);
    H(N + ii, 4) = u(1);
    H(N + ii, 5) = u(2);
    dy(ii) = psr(ii) - pred_psr;
    dy(N + ii) = psrdot(ii) - pred_psrdot;

    for (int jj = 1; jj < n_ant; jj++) {
      k2 = M + (n_ant - 1) * ii + jj - 1;
      pred_phase = -(ant_ned.col(jj).transpose() * u)(0) / lamb;
      pred_phase = std::fmod(pred_phase + navtools::PI<>, navtools::TWO_PI<>) - navtools::PI<>;
      if (pred_phase < -navtools::PI<>) {
        pred_phase += navtools::TWO_PI<>;
      } else if (pred_phase > navtools::PI<>) {
        pred_phase -= navtools::TWO_PI<>;
      }
      // hp = -(navtools::Skew<double>(ant_ned.col(j)).transpose() * u) / lamb;
      hp = (navtools::Skew(ant_ned.col(jj)).transpose() * u) / lamb;

      H(k2, 6) = hp(0);
      H(k2, 7) = hp(1);
      H(k2, 8) = hp(2);
      dy(k2) = phase(jj, ii) - pred_phase;
      dy(k2) = std::fmod(dy(k2) + navtools::PI<>, navtools::TWO_PI<>) - navtools::PI<>;
      if (dy(k2) < -navtools::PI<>) {
        dy(k2) += navtools::TWO_PI<>;
      } else if (dy(k2) > navtools::PI<>) {
        dy(k2) -= navtools::TWO_PI<>;
      }
      // std::cout << "meas_phase(" << k2 << "): " << phase(jj, ii) << " | est_phase(" << k2
      //           << "): " << pred_phase << " | dy(" << k2 << "): " << dy(k2) << "\n";
    }
  }

  // std::cout << "H: \n" << H << "\n";
  // std::cout << "dy: \n" << dy.segment(M - 1, MM - M).transpose() << "\n";

  // === Kalman Update ===
  KalmanUpdate(R, H, dy);
}

void KinematicNav::KalmanUpdate(
    const Eigen::Ref<const Eigen::MatrixXd> &R,
    const Eigen::Ref<const Eigen::MatrixXd> &H,
    const Eigen::Ref<const Eigen::VectorXd> &dy) {
  // // innovation filter
  // Eigen::MatrixXd PHt = P_ * H.transpose();
  // Eigen::VectorXd sqrt_S = (H * PHt + R).diagonal().cwiseSqrt();
  // Eigen::VectorX<bool> mask = (dy.array().abs() / sqrt_S.array()) < 3.0;
  // const int new_M = (mask).count();
  // if (new_M < M) {
  //   if (new_M > 0) {
  //     Eigen::MatrixXd new_H(new_M, 8);
  //     Eigen::MatrixXd new_R(new_M, new_M);
  //     Eigen::VectorXd new_dy(new_M);
  //     int k = 0;
  //     for (int i = 0; i < M; i++) {
  //       if (mask(i)) {
  //         new_H.row(k) = H.row(i);
  //         new_R(k, k) = R(i, i);
  //         new_dy(k) = dy(i);
  //         k++;
  //       }
  //     }

  //     // === Kalman Update ===
  //     Eigen::MatrixXd PHt = P_ * new_H.transpose();
  //     Eigen::MatrixXd K = PHt * (new_H * PHt + new_R).inverse();
  //     Eigen::MatrixXd L = I11_ - K * new_H;
  //     P_ = L * P_ * L.transpose() + K * new_R * K.transpose();
  //     x_ += K * new_dy;
  //   }
  // } else {
  //   // === Kalman Update ===
  //   Eigen::MatrixXd PHt = P_ * H.transpose();
  //   Eigen::MatrixXd K = PHt * (H * PHt + R).inverse();
  //   Eigen::MatrixXd L = I11_ - K * H;
  //   P_ = L * P_ * L.transpose() + K * R * K.transpose();
  //   x_ += K * dy;
  // }

  Eigen::MatrixXd PHt = P_ * H.transpose();
  Eigen::MatrixXd K = PHt * (H * PHt + R).inverse();
  Eigen::MatrixXd L = I11_ - K * H;
  // if (!is_init_) {
  //   for (int i = 0; i < 100; i++) {
  //     P_ = L * P_ * L.transpose() + K * R * K.transpose();
  //     P_ = F_ * P_ * F_.transpose() + Q_;
  //     PHt = P_ * H.transpose();
  //     K = PHt * (H * PHt + R).inverse();
  //     L = I11_ - K * H;
  //   }
  //   is_init_ = true;
  // } else {
  P_ = L * P_ * L.transpose() + K * R * K.transpose();
  // }
  x_ += K * dy;

  // std::cout << "dy = \n" << dy.transpose() << "\ndx = \n" << x_.transpose() << "\n";
  // std::cout << "H = \n" << H << "\n";

  // Closed loop error corrections
  phi_ += x_(0) / Hn_;
  lam_ += x_(1) / (He_ * cL_);
  h_ -= x_(2);
  vn_ += x_(3);
  ve_ += x_(4);
  vd_ += x_(5);
  Eigen::Vector4d q_err{1.0, 0.5 * x_(6), 0.5 * x_(7), 0.5 * x_(8)};
  q_b_l_ = navtools::quatdot<double>(q_err, q_b_l_);
  q_b_l_ /= q_b_l_.norm();
  navtools::quat2dcm<double>(C_b_l_, q_b_l_);
  cb_ += x_(9);
  cd_ += x_(10);
  x_.setZero();

  // std::cout << "P = \n" << P_.diagonal().transpose() << "\n";
}

// //! ITERATIVE EKF
// double sLam, cLam, t, sqt, pred_psr, pred_psrdot, pred_phase;
// double lat = phi_, lon = lam_, h = h_, vn = vn_, ve = ve_, vd = vd_, cb = cb_, cd = cd_;
// Eigen::Matrix3d C_l_e, C_e_l, C_b_l = C_b_l_;
// Eigen::Vector4d q_err, q_b_l = q_b_l_;
// Eigen::Vector3d u, udot, hp, ecef_p, ecef_v;
// Eigen::Matrix3Xd ant_ned(3, n_ant);
// Eigen::MatrixXd PHt(11, MM), K(11, MM);
// Eigen::VectorXd x_k = x_;
// Eigen::VectorXd delta(11);
// for (int kk = 0; kk < 10; kk++) {
//   sLam = std::sin(lon);
//   cLam = std::cos(lon);
//   sL_ = std::sin(lat);
//   cL_ = std::cos(lat);
//   sLsq_ = sL_ * sL_;
//   C_l_e << -sL_ * cLam, -sLam, -cL_ * cLam, -sL_ * sLam, cLam, -cL_ * sLam, cL_, 0.0, -sL_;
//   C_e_l = C_l_e.transpose();

//   t = 1.0 - navtools::WGS84_E2<> * sLsq_;
//   sqt = std::sqrt(t);
//   Re_ = navtools::WGS84_R0<> / sqt;
//   Rn_ = navtools::WGS84_R0<> * X1ME2_ / (t * t / sqt);
//   He_ = Re_ + h;
//   Hn_ = Rn_ + h;

//   ecef_p << He_ * cL_ * cLam, He_ * cL_ * sLam, (Re_ * X1ME2_ + h) * sL_;
//   ecef_v << vn, ve, vd;
//   ecef_v = C_l_e * ecef_v;
//   ant_ned = C_b_l * ant_xyz;

//   for (int ii = 0; ii < N; ii++) {
//     // range and range rate
//     RangeAndRate(
//         ecef_p, ecef_v, cb, cd, sv_pos.col(ii), sv_vel.col(ii), u, udot, pred_psr, pred_psrdot);
//     u = C_e_l * u;
//     udot = C_e_l * udot;
//     H(ii, 0) = u(0);
//     H(ii, 1) = u(1);
//     H(ii, 2) = u(2);
//     H(N + ii, 0) = udot(0);
//     H(N + ii, 1) = udot(1);
//     H(N + ii, 2) = udot(2);
//     H(N + ii, 3) = u(0);
//     H(N + ii, 4) = u(1);
//     H(N + ii, 5) = u(2);
//     dy(ii) = psr(ii) - pred_psr;
//     dy(N + ii) = psrdot(ii) - pred_psrdot;

//     for (int jj = 1; jj < n_ant; jj++) {
//       k2 = M + (n_ant - 1) * ii + jj - 1;
//       pred_phase = -(ant_ned.col(jj).transpose() * u)(0) / lamb;
//       pred_phase = std::fmod(pred_phase + navtools::PI<>, navtools::TWO_PI<>) - navtools::PI<>;
//       if (pred_phase < -navtools::PI<>) {
//         pred_phase += navtools::TWO_PI<>;
//       } else if (pred_phase > navtools::PI<>) {
//         pred_phase -= navtools::TWO_PI<>;
//       }
//       // hp = -(navtools::Skew<double>(ant_ned.col(jj)).transpose() * u) / lamb;
//       hp = (navtools::Skew(ant_ned.col(jj)).transpose() * u) / lamb;

//       H(k2, 6) = hp(0);
//       H(k2, 7) = hp(1);
//       H(k2, 8) = hp(2);
//       dy(k2) = phase(jj, ii) - pred_phase;
//       dy(k2) = std::fmod(dy(k2) + navtools::PI<>, navtools::TWO_PI<>) - navtools::PI<>;
//       if (dy(k2) < -navtools::PI<>) {
//         dy(k2) += navtools::TWO_PI<>;
//       } else if (dy(k2) > navtools::PI<>) {
//         dy(k2) -= navtools::TWO_PI<>;
//       }
//       std::cout << "meas_phase(" << k2 << "): " << phase(jj, ii) << " | est_phase(" << k2
//                 << "): " << pred_phase << " | dy(" << k2 << "): " << dy(k2) << "\n";
//     }
//   }

//   // iterative corrections
//   PHt = P_ * H.transpose();
//   K = PHt * (H * PHt + R).inverse();
//   delta = K * (H * x_k + dy) - x_k;  // x_ is zero
//   x_k += delta;

//   // Closed loop error corrections
//   lat = phi_ + x_k(0) / Hn_;
//   lon = lam_ + x_k(1) / (He_ * cL_);
//   h = h_ - x_k(2);
//   vn = vn_ + x_k(3);
//   ve = ve_ + x_k(4);
//   vd = vd_ + x_k(5);
//   q_err << 1.0, 0.5 * x_k(6), 0.5 * x_k(7), 0.5 * x_k(8);
//   q_b_l = navtools::quatdot<double>(q_err, q_b_l_);
//   q_b_l /= q_b_l.norm();
//   navtools::quat2dcm<double>(C_b_l, q_b_l);
//   cb = cb_ + x_k(9);
//   cd = cd_ + x_k(10);

//   if (delta.norm() < 1e-6) break;
// }

// // finish Kalman update
// Eigen::MatrixXd L = I11_ - K * H;
// P_ = L * P_ * L.transpose() + K * R * K.transpose();
// x_.setZero();
// phi_ = lat;
// lam_ = lon;
// h_ = h;
// vn_ = vn;
// ve_ = ve;
// vd_ = vd;
// cb_ = cb;
// cd_ = cd;
// q_b_l_ = q_b_l;
// C_b_l_ = C_b_l;

}  // namespace sturdins