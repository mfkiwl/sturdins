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

#include <unsupported/Eigen/MatrixFunctions>

// #include <Eigen/Eigenvalues>
#include <complex>
#include <iostream>
#include <navtools/constants.hpp>
#include <navtools/math.hpp>

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

// *=== GnssPVT ===*
bool GnssPVT(
    Eigen::Ref<Eigen::VectorXd> x,
    Eigen::Ref<Eigen::MatrixXd> P,
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
  Eigen::MatrixXd H{Eigen::MatrixXd::Zero(M, 8)};
  for (int i = 0; i < N; i++) {
    H(i, 6) = 1.0;
    H(N + i, 7) = 1.0;
  }
  Eigen::MatrixXd Ht(H.rows(), H.cols());
  Eigen::MatrixXd W{Eigen::MatrixXd::Zero(M, M)};
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
    if (dx.squaredNorm() < 1e-6) {
      return true;
    }
  }

  // failed to converge in time!
  return false;
}

// *=== PhasedArrayAttitude ===*
bool PhasedArrayAttitude(
    Eigen::Ref<Eigen::Matrix3d> C_b_l,
    const Eigen::Ref<const Eigen::Matrix3Xd> &u_ned,
    const Eigen::Ref<const Eigen::MatrixXd> &meas_phase,
    const Eigen::Ref<const Eigen::MatrixXd> &meas_phase_var,
    const Eigen::Ref<const Eigen::MatrixXd> &ant_xyz,
    const int &n_ant,
    const double &lambda,
    const double &thresh) {
  const int N = meas_phase.cols();  // number of satellites
  const int M = N * n_ant;
  // std::cout << "M: " << M << "\n";

  int k;
  Eigen::MatrixXd R{Eigen::MatrixXd::Zero(M, M)};
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < n_ant; j++) {
      k = n_ant * i + j;
      R(k, k) = 1.0 / meas_phase_var(j, i);
    }
  }

  Eigen::Vector3d dx{Eigen::Vector3d::Ones()};
  Eigen::MatrixXd H{Eigen::MatrixXd::Zero(M, 3)};
  Eigen::VectorXd dy{Eigen::VectorXd::Zero(M)};
  Eigen::MatrixXd ant_ned(3, n_ant);
  double est_phase;
  for (int z = 0; z < 10; z++) {
    // construct innovation and observation matrix
    ant_ned = C_b_l * ant_xyz;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < n_ant; j++) {
        k = n_ant * i + j;
        // std::cout << "k: " << k << "\n";
        est_phase = (ant_ned.col(j).transpose() * u_ned.col(i))(0) / lambda;
        // std::cout << "meas_phase(" << k << "): " << meas_phase(j, i) << " | est_phase(" << k
        //           << "): " << est_phase << "\n";
        dy(k) = meas_phase(j, i) - est_phase;
        navtools::WrapPiToPi<double>(dy(k));
        H.row(k) = (navtools::Skew<double>(ant_ned.col(j)).transpose() * u_ned.col(i)) / lambda;
      }
    }

    // recursive least squares
    // std::cout << "\n";
    // std::cout << "dy: " << dy.transpose() << "\n";
    // std::cout << "H: \n" << H << "\n";
    dx = (H.transpose() * R * H).inverse() * H.transpose() * R * dy;
    // std::cout << "dx: " << dx.transpose() << "\n";
    C_b_l = navtools::Skew<double>(dx).exp() * C_b_l;

    if (dx.squaredNorm() < thresh) {
      std::cout << "z = " << z << "\n";
      return true;
    }
  }
  return false;
}

// *=== Wahba ===*
void Wahba(
    Eigen::Ref<Eigen::Matrix3d> C_l_b,
    const Eigen::Ref<const Eigen::Matrix3Xd> &u_body,
    const Eigen::Ref<const Eigen::Matrix3Xd> &u_ned,
    const Eigen::Ref<const Eigen::VectorXd> &u_body_var) {
  // const int N = u_body.cols();
  // Eigen::MatrixXd R(N, N);
  // R.diagonal() << 1.0 / u_body_var.array();
  // Eigen::MatrixXd u_ned_T = u_ned.transpose();
  // C_l_b = u_body * R * u_ned_T * (u_ned * R * u_ned_T).inverse();

  // F. Landis Markley - Attitude Determination using Vector Observations and Singular Value
  // Decomposition (1968)
  const int N = u_body.cols();
  Eigen::MatrixXd R{Eigen::MatrixXd::Zero(N, N)};
  R.diagonal() << 1.0 / u_body_var.array();
  Eigen::Matrix3d B = u_body * R * u_ned.transpose();
  // Eigen::Matrix3d B = u_body * u_ned.transpose();
  // std::cout << "B: \n" << B << "\n";
  Eigen::JacobiSVD<Eigen::MatrixXd> svd;
  svd.compute(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::Matrix3d U = svd.matrixU();
  // std::cout << "U: \n" << U << "\n";
  Eigen::Matrix3d V = svd.matrixV();
  // std::cout << "V: \n" << V << "\n";
  double d = U.determinant() * V.determinant();
  // std::cout << "d: \n" << d << "\n";
  Eigen::Matrix3d S{Eigen::Matrix3d::Zero()};
  S.diagonal() << 1.0, 1.0, d;
  // std::cout << "S: \n" << S << "\n";
  C_l_b = (U * S * V.transpose());
  // std::cout << "C_l_b: \n" << C_l_b << "\n";
}

// *=== MUSIC ===*
void MUSIC(
    double &az_mean,
    double &el_mean,
    const Eigen::Ref<const Eigen::VectorXcd> &P,
    const Eigen::Ref<const Eigen::Matrix3Xd> &ant_xyz,
    const int &n_ant,
    const double &lambda,
    const double &thresh) {
  // std::cout << "MUSIC inputs: \n\taz_mean: " << az_mean << "\n\tel_mean: " << el_mean
  //           << "\n\tP: " << P.transpose() << "\n\tant_xyz: \n"
  //           << ant_xyz << "\n\tn_ant: " << n_ant << "\n\tlambda: " << lambda
  //           << "\n\tthresh: " << thresh << "\n\n";

  // 1) calculate correlation/covariance
  Eigen::MatrixXcd S = (P * P.conjugate().transpose());
  // std::cout << "S: \n" << S << "\n";

  // 2) calculate eigen-structure of S
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigen_solver(S, true);
  Eigen::VectorXcd e = eigen_solver.eigenvalues();
  Eigen::MatrixXcd v = eigen_solver.eigenvectors();
  // std::cout << "e: \n" << e.transpose() << "\n";
  // std::cout << "v: \n" << v << "\n";

  // 3) determine number of incident signals
  Eigen::VectorX<bool> idx = (e.array().abs()) < 10.0;
  const int N = idx.count();
  // const int D = n_ant - N
  Eigen::MatrixXcd U(n_ant, N);
  for (int i = 0; i < n_ant; i++) {
    if (idx(i)) {
      U.col(i) = v.col(i);
    }
  }
  // std::cout << "U: \n" << U << "\n";

  // 4) Find max power given arrival angles
  double az_span = navtools::DEG2RAD<> * 180.0;
  double el_span = navtools::DEG2RAD<> * 45.0;
  double res = navtools::DEG2RAD<> * 1.0;
  az_mean = navtools::DEG2RAD<> * 0.0;
  el_mean = navtools::DEG2RAD<> * -45.0;
  Eigen::VectorXcd tmp(n_ant);
  Eigen::Vector3d u;
  double sin_az, cos_az, sin_el, cos_el;
  int max_az_idx, max_el_idx;
  // std::complex<double> xyz2phase = -navtools::COMPLEX_I<> * navtools::TWO_PI<> / lambda;
  std::complex<double> xyz2phase = -navtools::COMPLEX_I<> / lambda;
  while (res >= thresh) {
    Eigen::VectorXd az{Eigen::VectorXd::LinSpaced(
        2 * static_cast<int>(az_span / res) + 1, az_mean - az_span, az_mean + az_span)};
    Eigen::VectorXd el{Eigen::VectorXd::LinSpaced(
        2 * static_cast<int>(el_span / res) + 1, el_mean - el_span, el_mean + el_span)};
    Eigen::MatrixXd P_music(az.size(), el.size());

    for (int i = 0; i < az.size(); i++) {
      // loop through azimuths
      sin_az = std::sin(az(i));
      cos_az = std::cos(az(i));

      for (int j = 0; j < el.size(); j++) {
        // loop through elevations
        sin_el = std::sin(el(j));
        cos_el = std::cos(el(j));

        // calculate spatial phase
        u << cos_az * cos_el, sin_az * cos_el, -sin_el;
        for (int k = 0; k < n_ant; k++) {
          tmp(k) = std::exp(xyz2phase * u.dot(ant_xyz.col(k)));
        }

        // determine power based on input az & el
        P_music(i, j) =
            1.0 / std::abs((tmp.conjugate().transpose() * U * U.conjugate().transpose() * tmp)(0));
      }
    }

    // find peak, recenter, and increase resolution
    P_music.maxCoeff(&max_az_idx, &max_el_idx);
    az_mean = az(max_az_idx);
    el_mean = el(max_el_idx);
    az_span = res / 2.0;
    el_span = res / 2.0;
    res /= 10.0;
    // std::cout << "az & el mean: " << az_mean * navtools::RAD2DEG<> << ", "
    //           << el_mean * navtools::RAD2DEG<> << "\n";
  }
}

}  // namespace sturdins