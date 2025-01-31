/**
 * *inertial-filter.cpp*
 *
 * =======  ========================================================================================
 * @file    sturdins/inertial-filter.cpp
 * @brief   Inertial navigation Kalman Filter equations.
 * @date    January 2025
 * @author  Daniel Sturdivant <sturdivant20@gmail.com>
 * @ref     1. "Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems", 2nd
 *              Edition, 2013 - Groves
 *          2. "Fundamentals of Inertial Navigation, Satellite-based Positioning and their
 *              Integration" - Noureldin, Karamat, & Georgy
 * =======  ========================================================================================
 */

#include "sturdins/ins.hpp"

#include <Eigen/src/Core/Matrix.h>

#include <navtools/attitude.hpp>
#include <navtools/constants.hpp>
#include <navtools/math.hpp>

#include "sturdins/least-squares.hpp"
#include "sturdins/strapdown.hpp"

namespace sturdins {

// *=== Ins ===*
Ins::Ins()
    : Strapdown(),
      bg_{Eigen::Vector<double, 3>::Zero()},
      ba_{Eigen::Vector<double, 3>::Zero()},
      x_{Eigen::Vector<double, 17>::Zero()},
      P_{Eigen::Matrix<double, 17, 17>::Zero()},
      F_{Eigen::Matrix<double, 17, 17>::Zero()},
      Q_{Eigen::Matrix<double, 17, 17>::Zero()},
      I17_{Eigen::Matrix<double, 17, 17>::Identity()} {
}
Ins::Ins(
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
    : Strapdown(lat, lon, alt, veln, vele, veld, roll, pitch, yaw),
      bg_{Eigen::Vector<double, 3>::Zero()},
      ba_{Eigen::Vector<double, 3>::Zero()},
      cb_{cb},
      cd_{cd},
      x_{Eigen::Vector<double, 17>::Zero()},
      P_{Eigen::Matrix<double, 17, 17>::Zero()},
      F_{Eigen::Matrix<double, 17, 17>::Zero()},
      Q_{Eigen::Matrix<double, 17, 17>::Zero()},
      I17_{Eigen::Matrix<double, 17, 17>::Identity()} {
}

// *=== ~Ins ===*
Ins::~Ins() {
}

// *=== SetImuSpec ===*
void Ins::SetImuSpec(const double &Ba, const double &Na, const double &Bg, const double &Ng) {
  double B_acc = Ba * 9.80665 / 1000.0;                // [mg] -> [(m/s)/s]
  double N_acc = Na / 60.0;                            // [(m/s)/√hr] -> [(m/s)/√s]
  double B_gyr = navtools::DEG2RAD<> * (Bg / 3600.0);  // [deg/hr] ->  [rad/s]
  double N_gyr = navtools::DEG2RAD<> * (Ng / 60.0);    // [deg/√hr]  ->  [rad/√s]

  // 1.21 to account for 10% stochastic noise
  Sra_ = 1.21 * N_acc * N_acc;
  Sbad_ = 1.21 * B_acc * B_acc;
  Srg_ = 1.21 * N_gyr * N_gyr;
  Sbgd_ = 1.21 * B_gyr * B_gyr;

  // FOR DYNAMIC BIASES
  // double K_acc = Ka / std::pow(3600, 1.5);             // [(m/s)/(hr*√hr] -> [(m/s)/(s*√s)]
  // double K_gyr =
  //     navtools::DEG2RAD<> * (Kg / std::pow(3600.0, 1.5));  // [deg/(hr*√hr)] ->  [rad/(s*√s)]
}

// *=== SetClockSpec ===*
void Ins::SetClockSpec(const double &h0, const double &h1, const double &h2) {
  double LS2 = navtools::LIGHT_SPEED<> * navtools::LIGHT_SPEED<>;
  h0_ = LS2 * h0 / 2.0;
  h1_ = LS2 * 2.0 * h1;
  h2_ = LS2 * navtools::PI_SQU<> * h2;
}

// *=== SetClock ===*
void Ins::SetClock(const double &cb, const double &cd) {
  cb_ = cb;
  cd_ = cd;
}

// *=== Propagate ===*
void Ins::Propagate(const Eigen::Vector3d &wb, const Eigen::Vector3d &fb, const double &dt) {
  Rg_ = Re_ * std::sqrt(1.0 - sLsq_ * navtools::WGS84_E<> * (2.0 - navtools::WGS84_E<>));
  cLsq_ = cL_ * cL_;
  Hnsq_ = Hn_ * Hn_;
  Hesq_ = He_ * He_;
  double vnve_ = vn_ * ve_;
  double vn_Hn_ = -w_en_n_(1);
  double ve_He_ = w_en_n_(0);
  double ve_Hn_ = ve_ / Hn_;
  double vnsq_ = vn_ * vn_;
  double vesq_ = ve_ * ve_;
  double HnHe_ = Hn_ * He_;
  double vewie = ve_ * navtools::WGS84_OMEGA<>;
  double vnwie = vn_ * navtools::WGS84_OMEGA<>;
  double two_wie = 2.0 * navtools::WGS84_OMEGA<>;

  /**
   * @brief First order F/Phi matrix from Groves Ch.14 and appendix I
   * --                                            --
   * | I3 + F33*T     F32*T       Z3        Z3   Z3 |
   * |    F23*T    I3 + F22*T    F21*T     C*T   Z3 |
   * |    F13*T       F12*T    I3 + F11*T   Z3  C*T |
   * |      Z3          Z3         Z3       I3   Z3 |
   * |      Z3          Z3         Z3       Z3   I3 |
   * --                                            --
   */
  Eigen::Vector3d wn = C_b_l_ * wb - (w_ie_n_ + w_en_n_);
  Eigen::Vector3d fn = C_b_l_ * fb;

  // F33
  F_(0, 0) = 1.0;
  F_(0, 2) = vn_Hn_ * dt;
  F_(1, 0) = ve_Hn_ * tL_ * dt;
  F_(1, 1) = 1.0;
  F_(1, 2) = ve_He_ * dt;
  F_(2, 2) = 1.0;

  // F32*T
  F_(0, 3) = dt;
  F_(1, 4) = dt;
  F_(2, 5) = dt;

  // F23*T
  F_(3, 0) = -(vesq_ / HnHe_ / cLsq_ + 2.0 * vewie * cL_ / Hn_) * dt;
  F_(3, 2) = (-vesq_ * tL_ / Hesq_ + vn_ * vd_ / Hnsq_) * dt;
  F_(4, 0) = (vnve_ / HnHe_ / cLsq_ + 2.0 * (vnwie * cL_ - vd_ * sL_) / Hn_) * dt;
  F_(4, 2) = ((vnve_ * tL_ + ve_ * vd_) / Hesq_) * dt;
  F_(5, 0) = (2.0 * vewie * sL_ / Hn_) * dt;
  F_(5, 2) = (-vesq_ / Hesq_ - vnsq_ / Hnsq_ + 2.0 * g0_ / Rg_) * dt;

  // I3 + F22 * T
  F_(3, 3) = 1.0 + (vd_ / Hn_) * dt;
  F_(3, 4) = -(2.0 * ve_He_ * tL_ + two_wie * sL_) * dt;
  F_(3, 5) = vn_Hn_ * dt;
  F_(4, 3) = (ve_He_ * tL_ + 2.0 * two_wie * sL_) * dt;
  F_(4, 4) = 1.0 + ((vn_ * tL_ + vd_) / He_) * dt;
  F_(4, 5) = (ve_He_ + two_wie * cL_) * dt;
  F_(5, 3) = -2.0 * vn_Hn_ * dt;
  F_(5, 4) = -(2.0 * ve_He_ + two_wie * cL_) * dt;
  F_(5, 5) = 1.0;

  // F21*T
  F_(3, 7) = fn(2) * dt;
  F_(3, 8) = -fn(1) * dt;
  F_(4, 6) = -fn(2) * dt;
  F_(4, 8) = fn(0) * dt;
  F_(5, 6) = fn(1) * dt;
  F_(5, 7) = -fn(0) * dt;

  // C*T
  F_(3, 9) = C_b_l_(0, 0) * dt;
  F_(3, 10) = C_b_l_(0, 1) * dt;
  F_(3, 11) = C_b_l_(0, 2) * dt;
  F_(4, 9) = C_b_l_(1, 0) * dt;
  F_(4, 10) = C_b_l_(1, 1) * dt;
  F_(4, 11) = C_b_l_(1, 2) * dt;
  F_(5, 9) = C_b_l_(2, 0) * dt;
  F_(5, 10) = C_b_l_(2, 1) * dt;
  F_(5, 11) = C_b_l_(2, 2) * dt;

  // F13*T
  F_(6, 0) = (navtools::WGS84_OMEGA<> * sL_ / Hn_) * dt;
  F_(6, 2) = (-ve_ / Hesq_) * dt;
  F_(7, 2) = (vn_ / Hnsq_) * dt;
  F_(8, 0) = ((navtools::WGS84_OMEGA<> * cL_ + ve_He_ / cLsq_) / Hn_) * dt;
  F_(8, 2) = (ve_ * tL_ / Hesq_) * dt;

  // F12*T
  F_(6, 4) = -dt / He_;
  F_(7, 3) = dt / Hn_;
  F_(8, 4) = dt * tL_ / He_;

  // I3 + F11*T
  F_(6, 6) = 1.0;
  F_(6, 7) = wn(2) * dt;
  F_(6, 8) = -wn(1) * dt;
  F_(7, 6) = -wn(2) * dt;
  F_(7, 7) = 1.0;
  F_(7, 8) = wn(0) * dt;
  F_(8, 6) = wn(1) * dt;
  F_(8, 7) = -wn(0) * dt;
  F_(8, 8) = 1.0;

  // C*T
  F_(6, 12) = F_(3, 9);
  F_(6, 13) = F_(3, 10);
  F_(6, 14) = F_(3, 11);
  F_(7, 12) = F_(4, 9);
  F_(7, 13) = F_(4, 10);
  F_(7, 14) = F_(4, 11);
  F_(8, 12) = F_(5, 9);
  F_(8, 13) = F_(5, 10);
  F_(8, 14) = F_(5, 11);

  // I3
  F_(9, 9) = 1.0;
  F_(10, 10) = 1.0;
  F_(11, 11) = 1.0;

  // I3
  F_(12, 12) = 1.0;
  F_(13, 13) = 1.0;
  F_(14, 14) = 1.0;

  // Clock F
  F_(15, 15) = 1.0;
  F_(15, 16) = dt;
  F_(16, 16) = 1.0;

  /**
   * @brief First order simplified Q matrix from Groves Ch.14
   * --                                   --
   * | Srg*I3    Z3    Z3     Z3      Z3   |
   * |  Z3     Sra*I3  Z3     Z3      Z3   |
   * |  Z3       Z3    Z3     Z3      Z3   |
   * |  Z3       Z3    Z3  Sbad*I3    Z3   |
   * |  Z3       Z3    Z3     Z3    Sgd*I3 |
   * --                                   --
   */
  double dtsq = dt * dt;
  double dtcb = dtsq * dt;
  Q_(0, 0) = 0.5 * Srg_ * dt;
  Q_(1, 1) = Q_(0, 0);
  Q_(2, 2) = Q_(0, 0);
  Q_(3, 3) = 0.5 * Sra_ * dt;
  Q_(4, 4) = Q_(3, 3);
  Q_(5, 5) = Q_(3, 3);
  Q_(9, 9) = 0.5 * Sbad_ * dt;
  Q_(10, 10) = Q_(9, 9);
  Q_(11, 11) = Q_(9, 9);
  Q_(12, 12) = 0.5 * Sbgd_ * dt;
  Q_(13, 13) = Q_(12, 12);
  Q_(14, 14) = Q_(12, 12);

  // Clock Q
  Q_(15, 15) = 0.5 * (h0_ * dt) + (h1_ * dtsq) + (2.0 / 3.0 * h2_ * dtcb);
  Q_(15, 16) = 0.5 * (h1_ * dt) + (h2_ * dtsq);
  Q_(16, 15) = Q_(6, 7);
  Q_(16, 16) = 0.5 * (h0_ / dt) + h1_ + (8.0 / 3.0 * h2_ * dt);

  // === Kalman Propagation ===
  P_ = F_ * (P_ + Q_) * F_.transpose() + Q_;
  // x_ = F_ * x_;
}

// *=== GnssUpdate ===*
void Ins::GnssUpdate(
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
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(M, 17);
  for (int i = 0; i < N; i++) {
    H(i, 15) = 1.0;
    H(N + i, 16) = 1.0;
  }
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(M, M);
  R.diagonal() << psr_var, psrdot_var;

  // //! Iterative EKF
  // Eigen::Vector3d u, udot;
  // Eigen::Vector3d ecef_p, ecef_v;
  // double pred_psr, pred_psrdot, cLam, sLam, t, sqt;
  // double lat = phi_, lon = lam_, alt = h_, veln = vn_, vele = ve_, veld = vd_, clkb = cb_,
  //        clkd = cd_;
  // Eigen::Matrix3d C_l_e, C_e_l;
  // Eigen::MatrixXd PHt(8, M), K(8, M);
  // Eigen::VectorXd x_k = x_;
  // Eigen::VectorXd delta(8);
  // for (int k = 0; k < 10; k++) {
  //   // update user position
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
  //   He_ = Re_ + alt;
  //   Hn_ = Rn_ + alt;
  //   ecef_p << He_ * cL_ * cLam, He_ * cL_ * sLam, (Re_ * X1ME2_ + alt) * sL_;
  //   ecef_v << veln, vele, veld;
  //   ecef_v = C_l_e * ecef_v;

  //   // iterate
  //   for (int i = 0; i < N; i++) {
  //     RangeAndRate(
  //         ecef_p, ecef_v, clkb, clkd, sv_pos.col(i), sv_vel.col(i), u, udot, pred_psr,
  //         pred_psrdot);
  //     u = -C_e_l * u;
  //     udot = -C_e_l * udot;
  //     H(i, 0) = u(0);
  //     H(i, 1) = u(1);
  //     H(i, 2) = u(2);
  //     H(N + i, 0) = udot(0);
  //     H(N + i, 1) = udot(1);
  //     H(N + i, 2) = udot(2);
  //     H(N + i, 3) = u(0);
  //     H(N + i, 4) = u(1);
  //     H(N + i, 5) = u(2);
  //     dy(i) = psr(i) - pred_psr;
  //     dy(N + i) = psrdot(i) - pred_psrdot;
  //   }

  //   // iterative corrections
  //   PHt = P_ * H.transpose();
  //   K = PHt * (H * PHt + R).inverse();
  //   delta = K * (H * x_k + dy) - x_k;  // x_ is zero
  //   x_k += delta;

  //   // Closed loop error corrections
  //   lat = phi_ - x_k(0) / Hn_;
  //   lon = lam_ - x_k(1) / (He_ * cL_);
  //   alt = h_ + x_k(2);
  //   veln = vn_ - x_k(3);
  //   vele = ve_ - x_k(4);
  //   veld = vd_ - x_k(5);
  //   clkb = cb_ + x_k(6);
  //   clkd = cd_ + x_k(7);

  //   if (delta.squaredNorm() < 1e-10) break;
  // }

  // // finish Kalman update
  // Eigen::MatrixXd L = I17_ - K * H;
  // P_ = L * P_ * L.transpose() + K * R * K.transpose();
  // x_.setZero();
  // phi_ = lat;
  // lam_ = lon;
  // h_ = alt;
  // vn_ = veln;
  // ve_ = vele;
  // vd_ = veld;
  // cb_ = clkb;
  // cd_ = clkd;
  // Eigen::Vector4d q_err{1.0, -0.5 * x_k(6), -0.5 * x_k(7), -0.5 * x_k(8)};
  // q_b_l_ = navtools::quatdot<double>(q_err, q_b_l_);
  // navtools::quat2dcm<double>(C_b_l_, q_b_l_);
  // ba_(0) += x_k(9);
  // ba_(1) += x_k(10);
  // ba_(2) += x_k(11);
  // bg_(0) += x_k(12);
  // bg_(1) += x_k(13);
  // bg_(2) += x_k(14);

  //! Regular EKF
  // Functions of current position
  double sLam = std::sin(lam_);
  double cLam = std::cos(lam_);
  sL_ = std::sin(phi_);
  cL_ = std::cos(phi_);
  RadiiOfCurvature<false>();  // Re_ and Rn_
  Hn_ = Rn_ + h_;
  He_ = Re_ + h_;
  Eigen::Matrix3d C_l_e{
      {-sL_ * cLam, -sLam, -cL_ * cLam}, {-sL_ * sLam, cLam, -cL_ * sLam}, {cL_, 0.0, -sL_}};
  Eigen::Matrix3d C_e_l = C_l_e.transpose();

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
    H(N + i, 3) = u(0);
    H(N + i, 4) = u(1);
    H(N + i, 5) = u(2);
    dy(i) = psr(i) - pred_psr;
    dy(N + i) = psrdot(i) - pred_psrdot;
  }

  // innovation filter
  Eigen::MatrixXd PHt = P_ * H.transpose();
  Eigen::VectorXd sqrt_S = (H * PHt + R).diagonal().cwiseSqrt();
  Eigen::VectorX<bool> mask = (dy.array().abs() / sqrt_S.array()) < 3.0;
  const int new_M = (mask).count();
  if (new_M < M) {
    Eigen::MatrixXd new_H(new_M, 8);
    Eigen::MatrixXd new_R(new_M, new_M);
    Eigen::VectorXd new_dy(new_M);
    int k = 0;
    for (int i = 0; i < M; i++) {
      if (mask(i)) {
        new_H.row(k) = H.row(i);
        new_R(k, k) = R(i, i);
        new_dy(k) = dy(i);
        k++;
      }
    }

    // === Kalman Update ===
    Eigen::MatrixXd PHt = P_ * new_H.transpose();
    Eigen::MatrixXd K = PHt * (new_H * PHt + new_R).inverse();
    Eigen::MatrixXd L = I17_ - K * new_H;
    P_ = L * P_ * L.transpose() + K * new_R * K.transpose();
    x_ += K * new_dy;
  } else {
    // === Kalman Update ===
    Eigen::MatrixXd PHt = P_ * H.transpose();
    Eigen::MatrixXd K = PHt * (H * PHt + R).inverse();
    Eigen::MatrixXd L = I17_ - K * H;
    P_ = L * P_ * L.transpose() + K * R * K.transpose();
    x_ += K * dy;
  }

  // Closed loop error corrections
  phi_ -= x_(0) / Hn_;
  lam_ -= x_(1) / (He_ * cL_);
  h_ -= -x_(2);
  vn_ -= x_(3);
  ve_ -= x_(4);
  vd_ -= x_(5);
  Eigen::Vector4d q_err{1.0, -0.5 * x_(6), -0.5 * x_(7), -0.5 * x_(8)};
  q_b_l_ = navtools::quatdot<double>(q_err, q_b_l_);
  navtools::quat2dcm<double>(C_b_l_, q_b_l_);
  ba_(0) += x_(9);
  ba_(1) += x_(10);
  ba_(2) += x_(11);
  bg_(0) += x_(12);
  bg_(1) += x_(13);
  bg_(2) += x_(14);
  cb_ += x_(15);
  cd_ += x_(16);
  x_.setZero();
}

}  // namespace sturdins