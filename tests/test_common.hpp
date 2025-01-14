
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <navtools/types.hpp>
#include <random>
#include <satutils/ephemeris.hpp>

#include "navtools/constants.hpp"
#include "sturdins/least-squares.hpp"

// normal distribution random number generator
std::default_random_engine noise_gen;
std::normal_distribution<double> noise_dist(0.0, 1.0);

// automotive grade IMU
inline constexpr double Ba = 1.2;
inline constexpr double Ka = 0.0;
inline constexpr double Na = 0.5884;
inline constexpr double Bg = 180.0;
inline constexpr double Kg = 0.0;
inline constexpr double Ng = 3.0;
inline constexpr double dt = 0.01;
inline constexpr double B_acc = Ba * 9.80665 / 1000.0;  // [mg] -> [(m/s)/s]
inline constexpr double K_acc = Ka / 216000.0;          // [(m/s)/(hr*√hr] -> [(m/s)/(s*√s)]
inline constexpr double N_acc = Na / 60.0;              // [(m/s)/√hr] -> [(m/s)/√s]
inline constexpr double Tc_acc = 100.0;                 // accel correlation times [s]
inline constexpr double B_gyr = navtools::DEG2RAD<> * (Bg / 3600.0);  // [deg/hr] ->  [rad/s]
inline constexpr double K_gyr =
    navtools::DEG2RAD<> * (Kg / 216000.0);  // [deg/(hr*√hr)] ->  [rad/(s*√s)]
inline constexpr double N_gyr = navtools::DEG2RAD<> * (Ng / 60.0);  // [deg/√hr]  ->  [rad/√s]
inline constexpr double Tc_gyr = 300.0;                             // gyro correlation times [s]
inline constexpr double alpha_a = 1.0 - dt / Tc_acc;
inline constexpr double alpha_g = 1.0 - dt / Tc_gyr;
inline const double beta_a = B_acc * std::sqrt(1.0 - std::exp(-2.0 * dt / Tc_acc));
inline const double beta_g = B_gyr * std::sqrt(1.0 - std::exp(-2.0 * dt / Tc_gyr));
inline const double gamma_a = N_acc / std::sqrt(dt);
inline const double gamma_g = N_gyr / std::sqrt(dt);

// high quality tcxo clock model
inline constexpr double LS2 = navtools::LIGHT_SPEED<> * navtools::LIGHT_SPEED<>;
inline constexpr double h0 = 2e-21;
inline constexpr double h1 = 1e-22;
inline constexpr double h2 = 2e-20;

// navigation data structure of binary files
template <typename T = double>
struct NavData {
  T t;
  T lat;
  T lon;
  T h;
  T vn;
  T ve;
  T vd;
  T roll;
  T pitch;
  T yaw;
  T fx;
  T fy;
  T fz;
  T wx;
  T wy;
  T wz;
};

template <typename T = double>
struct NavResult {
  T t;
  T lat;
  T lon;
  T h;
  T vn;
  T ve;
  T vd;
  T roll;
  T pitch;
  T yaw;
  T cb;
  T cd;
};

// function to parse ephemeris binary file
template <typename T = double>
std::vector<satutils::KeplerEphem<T>> ParseEphemeris(std::string filename) {
  // read ephemeris
  std::ifstream fid(filename, std::ios::binary);
  if (!fid) {
    std::cerr << "Error opening file!\n";
  }
  std::vector<satutils::KeplerEphem<T>> eph;
  satutils::KeplerElements<T> elem;
  while (fid.read(reinterpret_cast<char *>(&elem), sizeof(elem))) {
    std::cout << "Next ephemerides: ";
    std::cout << "\n\tiode     : " << elem.iode << "\n\tiodc     : " << elem.iodc
              << "\n\ttoe      : " << elem.toe << "\n\ttoc      : " << elem.toc
              << "\n\ttgd      : " << elem.tgd << "\n\toaf2     : " << elem.af2
              << "\n\taf1      : " << elem.af1 << "\n\taf0      : " << elem.af0
              << "\n\te        : " << elem.e << "\n\tsqrtA    : " << elem.sqrtA
              << "\n\tdeltan   : " << elem.deltan << "\n\tm0       : " << elem.m0
              << "\n\tomega0   : " << elem.omega0 << "\n\tomega    : " << elem.omega
              << "\n\tomegaDot : " << elem.omegaDot << "\n\ti0       : " << elem.i0
              << "\n\tiDot     : " << elem.iDot << "\n\tcuc      : " << elem.cuc
              << "\n\tcus      : " << elem.cus << "\n\tcic      : " << elem.cic
              << "\n\tcis      : " << elem.cis << "\n\tcrc      : " << elem.crc
              << "\n\tcrs      : " << elem.crs << "\n\tura      : " << elem.ura
              << "\n\thealth   : " << elem.health << "\n\n";
    eph.push_back(satutils::KeplerEphem<T>(elem));
  }
  fid.close();

  return eph;
};

// measurement data structure
struct MeasurementData {
  Eigen::MatrixXd sv_pos;
  Eigen::MatrixXd sv_vel;
  Eigen::VectorXd psr;
  Eigen::VectorXd psrdot;

  MeasurementData(int size) : sv_pos(3, size), sv_vel(3, size), psr(size), psrdot(size){};
};

// measurement model and awgn noise
MeasurementData MeasurementModel(
    const double &transmit_time,
    const double &psr_std,
    const double &psrdot_std,
    const Eigen::Vector3d &pos,
    const Eigen::Vector3d &vel,
    const double &cb,
    const double &cd,
    std::vector<satutils::KeplerEphem<double>> &eph) {
  const int N = eph.size();
  MeasurementData data(N);
  Eigen::MatrixXd sv_clk(3, N);
  Eigen::MatrixXd sv_acc(3, N);
  Eigen::Vector3d u, udot;
  for (int i = 0; i < N; i++) {
    // get satellite positions
    eph[i].CalcNavStates<false>(
        sv_clk.col(i), data.sv_pos.col(i), data.sv_vel.col(i), sv_acc.col(i), transmit_time);

    // get range and rate measurements
    sturdins::RangeAndRate(
        pos,
        vel,
        cb,
        cd,
        data.sv_pos.col(i),
        data.sv_vel.col(i),
        u,
        udot,
        data.psr(i),
        data.psrdot(i));

    // add noise to measurement
    data.psr(i) += psr_std * noise_dist(noise_gen);
    data.psrdot(i) += psrdot_std * noise_dist(noise_gen);
  }

  return data;
};

// clock awgn model
void ClockModel(Eigen::Vector2d &x, const double &T) {
  double T2 = T * T;
  double T3 = T2 * T;
  double q_bb =
      std::sqrt(LS2 * (h0 / 2 * T + 2 * h1 * T2 + (2.0 / 3.0) * navtools::PI_SQU<> * h2 * T3));
  double q_dd =
      std::sqrt(LS2 * (h0 / (2 * T) + 2 * h1 + (8.0 / 3.0) * navtools::PI_SQU<> * h2 * T));
  double q_bd = std::sqrt(LS2 * (2 * h1 * T + navtools::PI_SQU<> * h2 * T2));

  Eigen::Matrix2d A{{1.0, T}, {0.0, 1.0}};
  Eigen::Matrix2d Q{{q_bb, q_bd}, {q_bd, q_dd}};
  Eigen::Vector2d n{noise_dist(noise_gen), noise_dist(noise_gen)};
  x = A * x + Q * n;
}

// imu awgn noise model
void ImuModel(
    Eigen::Vector3d &w, Eigen::Vector3d &f, Eigen::Vector3d &drift_g, Eigen::Vector3d &drift_a) {
  double wn_a, wn_g;
  for (int i = 0; i < 3; i++) {
    // white noise
    wn_a = gamma_a * noise_dist(noise_gen);
    wn_g = gamma_g * noise_dist(noise_gen);

    // random walk
    drift_a(i) = (alpha_a * drift_a(i)) + (beta_a * noise_dist(noise_gen));
    drift_g(i) = (alpha_g * drift_g(i)) + (beta_g * noise_dist(noise_gen));

    // misalignment and scale factor ...

    // combine
    f(i) += (drift_a(i) + wn_a);
    w(i) += (drift_g(i) + wn_g);
  }
};
