
#include <cmath>
#include <iostream>
#include <navtools/attitude.hpp>
#include <navtools/constants.hpp>
#include <navtools/frames.hpp>

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
      2.5 * navtools::DEG2RAD<>, -14.1 * navtools::DEG2RAD<>, 65.0 * navtools::DEG2RAD<>};
  Eigen::MatrixXd ant_xyz{
      {0.0, 0.09514, 0.0, 0.09514}, {0.0, 0.0, -0.09514, -0.09514}, {0.0, 0.0, 0.0, 0.0}};
  const int N = eph.size();
  const int n_ant = 4;

  // get user ecef states
  std::cout << "rotating user position to ecef\n";
  Eigen::Matrix3d C_b_l = navtools::euler2dcm<true, double>(rpy);
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
        Eigen::Vector3d aer = navtools::ecef2aer<double>(ant_xyz_ecef.col(j), sv_pos.col(i));
        true_az(i) = aer(0);
        true_el(i) = aer(1);
      }
    }
  }
  std::cout << "true_az: \n\t" << true_az.array().transpose() * navtools::RAD2DEG<> << "\n";
  std::cout << "true_el: \n\t" << true_el.array().transpose() * navtools::RAD2DEG<> << "\n";

  // add noise to measurements
  std::cout << "adding noise to measurements\n";
  Eigen::MatrixXd code_range(n_ant, N), phase_range(n_ant, N), range_rate(n_ant, N);
  for (int i = 0; i < N; i++) {
    code_range.col(i) = true_psr.col(i).array() + 5.0 * noise_dist(noise_gen);
    phase_range.col(i) = true_psr.col(i).array() + 0.1 * noise_dist(noise_gen);
    range_rate.col(i) = true_psrdot.col(i).array() + 0.1 * noise_dist(noise_gen);
  }

  // simulate correlators
  std::cout << "simulating correlators\n";
  Eigen::MatrixXd chip_err = ((-code_range).rowwise() + true_psr.row(0)).array() / beta;
  Eigen::MatrixXd phase_err = ((-phase_range).rowwise() + true_psr.row(0)).array() / lamb;
  Eigen::MatrixXd freq_err = ((-range_rate).rowwise() + true_psrdot.row(0)).array() / lamb;
  double A = std::sqrt(2.0 * cno * T);
  Eigen::MatrixXd F = navtools::PI<> * freq_err * T;
  F = F.array().sin() / F.array();
  Eigen::MatrixXd R = 1.0 - chip_err.array().abs();
  Eigen::MatrixXcd P = (navtools::COMPLEX_I<> * 2.0 * navtools::PI<> * phase_err.array()).exp();
  Eigen::MatrixXcd prompt = A * R.array() * F.array() * P.array();

  // add noise to correlators
  std::cout << "adding noise to correlators\n";
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < n_ant; j++) {
      prompt(j, i) += std::complex<double>(noise_dist(noise_gen), noise_dist(noise_gen));
    }
  }

  // run music to predict azimuth and elevation
  std::cout << "starting music\n";
  Eigen::VectorXd est_az(N), est_el(N);
  for (int i = 0; i < N; i++) {
    sturdins::MUSIC(est_az(i), est_el(i), prompt.col(i), ant_xyz, n_ant, lamb, 1e-4);
  }
  std::cout << "est_az: \n" << est_az.array().transpose() * navtools::RAD2DEG<> << "\n";
  std::cout << "est_el: \n" << est_el.array().transpose() * navtools::RAD2DEG<> << "\n";

  return 0;
}