/**
 * *least-squares.cpp*
 *
 * =======  ========================================================================================
 * @file    sturdins/least-squares.cpp
 * @brief   Least squares algorithms for initializing the GNSS-INS filter.
 * @date    January 2025
 * @author  Daniel Sturdivant <sturdivant20@gmail.com>
 * @ref     1. "Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems", 2nd
 *              Edition, 2013 - Groves
 * =======  ========================================================================================
 */

#include "sturdins/least-squares.hpp"

#include <navtools/constants.hpp>

namespace sturdins {

// *=== RangeAndRate ===*
void RangeAndRate(
    const Eigen::Ref<const Eigen::Vector3d> &pos,
    const Eigen::Ref<const Eigen::Vector3d> &vel,
    const double &cb,
    const double &cd,
    const Eigen::Ref<const Eigen::Vector3d> &sv_pos,
    const Eigen::Ref<const Eigen::Vector3d> &sv_vel,
    Eigen::Ref<Eigen::Vector3d> u,
    Eigen::Ref<Eigen::Vector3d> udot,
    double &pred_psr,
    double &pred_psrdot) {
  // predict approximate range and account for earth's rotation (Groves 8.34)
  Eigen::Vector3d dr = pos - sv_pos;
  double r = dr.norm();
  double wtau = navtools::WGS84_OMEGA<double> * r / navtools::LIGHT_SPEED<double>;
  double swtau = std::sin(wtau);
  double cwtau = std::cos(wtau);
  Eigen::Matrix3d C_i_e{{cwtau, swtau, 0.0}, {-swtau, cwtau, 0.0}, {0.0, 0.0, 1.0}};

  // predict pseudorange (Groves 8.35, 9.165)
  dr = pos - C_i_e * sv_pos;
  r = dr.norm();
  u = dr.array() / r;
  pred_psr = r + cb;

  // predict pseudorange-rate (Groves 8.44, 9.165)
  double r2 = r * r;
  Eigen::Vector3d dv = (vel + navtools::WGS84_OMEGA_SKEW<double> * pos) -
                       C_i_e * (sv_vel + navtools::WGS84_OMEGA_SKEW<double> * sv_pos);
  udot = (dv * r2 - dr * dv.dot(dr)) / (r2 * r);
  pred_psrdot = u.dot(dv) + cd;
}

// *=== GaussNewton ===*
bool GaussNewton(
    Eigen::Ref<Eigen::VectorXd> x,
    Eigen::Ref<Eigen::MatrixXd> P,
    const Eigen::Ref<const Eigen::MatrixXd> &sv_pos,
    const Eigen::Ref<const Eigen::MatrixXd> &sv_vel,
    const Eigen::Ref<const Eigen::VectorXd> &psr,
    const Eigen::Ref<const Eigen::VectorXd> &psrdot,
    const Eigen::Ref<const Eigen::VectorXd> &psr_var,
    const Eigen::Ref<const Eigen::VectorXd> &psrdot_var) {
  // Initialize
  const int N = psr.size();
  const int M = 2 * N;
  Eigen::VectorXd dy(M);
  Eigen::MatrixXd H{Eigen::MatrixXd::Zero(M, 8)};
  for (int i = 0; i < N; i++) {
    H(i, 6) = 1.0;
    H(N + i, 7) = 1.0;
  }
  Eigen::MatrixXd Ht(H.rows(), H.cols());
  Eigen::MatrixXd W = Eigen::MatrixXd::Zero(M, M);
  W.diagonal() << 1.0 / psr_var.array(), 1.0 / psrdot_var.array();

  // Recursive Estimation
  Eigen::VectorXd dx(8);
  Eigen::Vector3d u, udot;
  double pred_psr, pred_psrdot;
  for (int k = 0; k < 10; k++) {  // should converge within 5 iterations

    // update predicted measurements based on updated state
    for (int i = 0; i < N; i++) {
      RangeAndRate(
          x.segment(0, 3),
          x.segment(3, 3),
          x(6),
          x(7),
          sv_pos.col(i),
          sv_vel.col(i),
          u,
          udot,
          pred_psr,
          pred_psrdot);
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

    // Gauss-Newton weighted least squares formula
    Ht = H.transpose();
    P = (Ht * W * H).inverse();
    dx = P * Ht * W * dy;
    x += dx;
    if ((dx(0) * dx(0) + dx(1) * dx(1) + dx(2) * dx(2)) < 1e-3) {
      return true;
    }
  }

  // failed to converge in time!
  return false;
}

}  // namespace sturdins