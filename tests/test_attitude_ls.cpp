
#include <cmath>
#include <iostream>
#include <navtools/attitude.hpp>
#include <navtools/constants.hpp>
#include <navtools/frames.hpp>
#include <navtools/math.hpp>

#include "sturdins/least-squares.hpp"
#include "test_common.hpp"

int main() {
  // parse ephemeris
  std::cout << "reading ephemerides\n";
  std::vector<satutils::KeplerEphem<double>> eph =
      ParseEphemeris<double>("src/sturdins/tests/sv_ephem.bin");

  // simplified initialization
  std::cout << "initializing\n";
  double lamb = navtools::LIGHT_SPEED<> / 1575.42e6;
  double beta = navtools::LIGHT_SPEED<> / 1.023e6;
  double cno = std::pow(10.0, 4.4);
  double T = 0.02;
  Eigen::Vector3d lla{32.5863 * navtools::DEG2RAD<>, -85.4943 * navtools::DEG2RAD<>, 213.0};
  Eigen::Vector3d nedv{0.0, 0.0, 0.0};
  Eigen::Vector3d rpy{
      2.5 * navtools::DEG2RAD<>, -15.1 * navtools::DEG2RAD<>, -67.9 * navtools::DEG2RAD<>};
  // Eigen::Vector3d rpy{-0.5, 0.5, 3.0};
  // Eigen::Vector3d rpy{0.0, 0.0, 0.0};
  Eigen::MatrixXd ant_xyz{
      {0.0, 0.09514, 0.0, 0.09514}, {0.0, 0.0, -0.09514, -0.09514}, {0.0, 0.0, 0.0, 0.0}};
  const int N = eph.size();
  const int n_ant = 4;

  // get user ecef states
  std::cout << "rotating user position to ecef\n";
  Eigen::Matrix3d C_b_l = navtools::euler2dcm<double>(rpy, true);
  Eigen::Matrix3d C_e_l = navtools::ecef2nedDcm<double>(lla);
  Eigen::MatrixXd ant_xyz_ecef(3, n_ant);
  Eigen::Vector3d ecefv = navtools::ned2ecefv<double>(nedv, lla);
  for (int j = 0; j < n_ant; j++) {
    ant_xyz_ecef.col(j) = navtools::ned2ecef<double>(C_b_l * ant_xyz.col(j), lla);
  }

  // get satellite ecef positions and velocities (simplified)
  std::cout << "calculating satellite positions\n";
  double ToW = 521400;
  Eigen::Vector3d sv_clk, sv_acc, u, udot;
  Eigen::MatrixXd sv_pos(3, N), sv_vel(3, N);
  Eigen::MatrixXd true_psr(n_ant, N), true_psrdot(n_ant, N);
  Eigen::VectorXd true_az(N), true_el(N);
  Eigen::MatrixXd u_ned(3, N);
  double cb = 0, cd = 0;
  for (int i = 0; i < N; i++) {
    eph[i].CalcNavStates<false>(sv_clk, sv_pos.col(i), sv_vel.col(i), sv_acc, ToW);

    // get range and rate measurements
    for (int j = 0; j < n_ant; j++) {
      sturdins::RangeAndRate(
          ant_xyz_ecef.col(j),
          ecefv,
          cb,
          cd,
          sv_pos.col(i),
          sv_vel.col(i),
          u,
          udot,
          true_psr(j, i),
          true_psrdot(j, i));
      if (j == 0) {
        // Eigen::Vector3d aer = navtools::ecef2aer<double>(ant_xyz_ecef.col(j), sv_pos.col(i));
        // true_az(i) = aer(0);
        // true_el(i) = aer(1);
        u_ned.col(i) = C_e_l * u;
        // u_ned.col(i) = -C_e_l * u;
        true_az(i) = std::atan2(u_ned(1, i), u_ned(0, i));
        true_el(i) = -std::asin(u_ned(2, i));
      }
    }
  }

  // add noise to measurements
  std::cout << "adding noise to measurements\n";
  Eigen::MatrixXd code_range(n_ant, N), phase_range(n_ant, N), range_rate(n_ant, N);
  for (int i = 0; i < N; i++) {
    code_range.col(i) = true_psr.col(i).array() + 5.0 * noise_dist(noise_gen);
    phase_range.col(i) = true_psr.col(i).array() + 20.0 * noise_dist(noise_gen);
    range_rate.col(i) = true_psrdot.col(i).array() + 0.1 * noise_dist(noise_gen);
  }

  // simulate correlators
  std::cout << "simulating correlators\n";
  // Eigen::MatrixXd chip_err = ((-code_range).rowwise() + true_psr.row(0)).array() / beta;
  // Eigen::MatrixXd phase_err = ((-phase_range).rowwise() + true_psr.row(0)).array() / lamb;
  // Eigen::MatrixXd freq_err = ((-range_rate).rowwise() + true_psrdot.row(0)).array() / lamb;
  Eigen::MatrixXd chip_err = (true_psr.rowwise() - code_range.row(0)).array() / beta;
  Eigen::MatrixXd phase_err = -(true_psr.rowwise() - phase_range.row(0)).array() / lamb;
  Eigen::MatrixXd freq_err = (true_psrdot.rowwise() - range_rate.row(0)).array() / lamb;
  double A = std::sqrt(2.0 * cno * T);
  Eigen::MatrixXd F = navtools::PI<> * freq_err * T;
  F = F.array().sin() / F.array();
  Eigen::MatrixXd R = 1.0 - chip_err.array().abs();
  Eigen::MatrixXcd P = (navtools::COMPLEX_I<> * navtools::TWO_PI<> * phase_err.array()).exp();
  Eigen::MatrixXcd prompt = A * R.array() * F.array() * P.array();

  // add noise to correlators
  std::cout << "adding noise to correlators\n";
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < n_ant; j++) {
      prompt(j, i) += std::complex<double>(noise_dist(noise_gen), noise_dist(noise_gen));
    }
  }
  std::cout << "PromptCorrelators: \n" << prompt.transpose() << "\n";

  // run music to predict azimuth and elevation
  std::cout << "\n--- starting music & solving wabhas problem ----------------------------------\n";
  Eigen::VectorXd est_az(N), est_el(N);
  for (int i = 0; i < N; i++) {
    sturdins::MUSIC(
        est_az(i), est_el(i), prompt.col(i), ant_xyz, n_ant, lamb / navtools::TWO_PI<>, 1e-4);
  }

  // run wahba to predict user attitude
  Eigen::Matrix3d C_l_b_est;
  Eigen::MatrixXd u_body_est(3, N);
  Eigen::VectorXd u_body_var{Eigen::VectorXd::Ones(N)};
  u_body_est.row(0) = est_az.array().cos() * est_el.array().cos();
  u_body_est.row(1) = est_az.array().sin() * est_el.array().cos();
  u_body_est.row(2) = -est_el.array().sin();
  sturdins::Wahba(C_l_b_est, u_body_est, u_ned, u_body_var);
  Eigen::Vector3d rpy_est = navtools::dcm2euler<double>(C_l_b_est.transpose(), true);
  // Eigen::Vector3d rpy_est = navtools::dcm2euler<true, double>(C_l_b_est);

  std::cout << "true_az: \n\t" << true_az.array().transpose() * navtools::RAD2DEG<> << "\n";
  std::cout << "true_el: \n\t" << true_el.array().transpose() * navtools::RAD2DEG<> << "\n";
  std::cout << "est_az: \n\t" << est_az.array().transpose() * navtools::RAD2DEG<> << "\n";
  std::cout << "est_el: \n\t" << est_el.array().transpose() * navtools::RAD2DEG<> << "\n";
  std::cout << "u_from_sv: \n" << u_ned.transpose() << "\n";
  std::cout << "u_from_music: \n" << u_body_est.transpose() << "\n";
  std::cout << "true_rpy: \n" << rpy * navtools::RAD2DEG<> << "\n";
  std::cout << "est_rpy: \n" << rpy_est * navtools::RAD2DEG<> << "\n";

  std::cout << "\n--- running alternative attitude estimator -----------------------------------\n";

  // calculate phase discriminators
  Eigen::MatrixXd phase_disc{Eigen::MatrixXd::Zero(n_ant, N)};
  // phase_disc = navtools::TWO_PI<> / lamb * true_psr;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < n_ant; j++) {
      phase_disc(j, i) = std::atan2(prompt(j, i).imag(), prompt(j, i).real());
    }
  }
  std::cout << "phase_disc: \n" << phase_disc.transpose() << "\n";

  Eigen::RowVectorXd phase_disc_0 = phase_disc.row(0);
  phase_disc = (-phase_disc).rowwise() + phase_disc_0;
  // phase_disc = (phase_disc).rowwise() - phase_disc_0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < n_ant; j++) {
      navtools::WrapPiToPi<double>(phase_disc(j, i));
    }
  }
  std::cout << "phase_disc: \n" << phase_disc.transpose() << "\n";

  // Eigen::Matrix3d C_b_l_est = C_l_b_est.transpose();
  Eigen::Matrix3d C_b_l_est{Eigen::Matrix3d::Identity()};
  Eigen::MatrixXd phase_disc_var{Eigen::MatrixXd::Ones(n_ant, N)};
  sturdins::PhasedArrayAttitude(
      C_b_l_est,
      u_ned,
      phase_disc,
      phase_disc_var,
      ant_xyz,
      n_ant,
      lamb / navtools::TWO_PI<>,
      1e-9);
  Eigen::Vector3d rpy_est2 = navtools::dcm2euler<double>(C_b_l_est, true);
  std::cout << "est_rpy2: \n" << rpy_est2 * navtools::RAD2DEG<> << "\n";

  return 0;
}