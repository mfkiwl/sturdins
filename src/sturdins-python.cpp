/**
 * *sturdins-python.cpp*
 *
 * =======  ========================================================================================
 * @file    src/sturdins-python.cpp
 * @brief   PyBind11 wrapper for using SturdINS in python!
 * @date    March 2025
 * =======  ========================================================================================
 */

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "sturdins/inertial-nav.hpp"
#include "sturdins/kinematic-nav.hpp"
#include "sturdins/least-squares.hpp"
#include "sturdins/nav-clock.hpp"
#include "sturdins/nav-imu.hpp"
#include "sturdins/strapdown.hpp"

namespace py = pybind11;
using namespace sturdins;

PYBIND11_MODULE(_sturdins_core, h) {
  h.doc() = R"pbdoc(
    SturdINS
    ========
    
    Daniel Sturdivant's GNSS-INS Sofware.
    
    Contains the following classes:
    
    1. `InertialNav`
    2. `KinematicNav`
    3. `Strapdown`

    Contains the following modules:

    1. `leastsquares`
    2. `navsense`

    )pbdoc";
  h.attr("__version__") = "1.0.0";

  // Strapdown
  py::class_<Strapdown>(h, "Strapdown")
      .def(py::init<>())
      .def(
          py::init<
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double>(),
          py::arg("lat"),
          py::arg("lon"),
          py::arg("alt"),
          py::arg("vn"),
          py::arg("ve"),
          py::arg("vd"),
          py::arg("r"),
          py::arg("p"),
          py::arg("y"))
      .def(
          "SetPosition",
          &Strapdown::SetPosition,
          py::arg("lat"),
          py::arg("lon"),
          py::arg("alt"),
          R"pbdoc(
          SetPosition
          ===========
          
          Set the position of the navigator
          
          Parameters
          ----------
          
          lat : double
          
              Latitude [rad]
          
          lon : double
          
              Longitude [rad]
              
          alt : double
          
              Altitude/Height [m]
          )pbdoc")
      .def(
          "SetVelocity",
          &Strapdown::SetVelocity,
          py::arg("vn"),
          py::arg("ve"),
          py::arg("vd"),
          R"pbdoc(
          SetVelocity
          ===========
          
          Set the velocity of the navigator
          
          Parameters
          ----------
          
          vn : double
          
              north velocity [m/s]
          
          ve : double
          
              east velocity [m/s]
              
          vd : double
          
              down velocity [m/s]
          )pbdoc")
      .def(
          "SetAttitude",
          py::overload_cast<const Eigen::Ref<const Eigen::Matrix3d> &>(&Strapdown::SetAttitude),
          py::arg("C"),
          R"pbdoc(
          SetAttitude
          ===========
          
          Set the attitude of the navigator
          
          Parameters
          ----------
          
          C : np.ndarray
          
              body-to-ned rotation matrix
          )pbdoc")
      .def(
          "SetAttitude",
          py::overload_cast<const double &, const double &, const double &>(
              &Strapdown::SetAttitude),
          py::arg("r"),
          py::arg("p"),
          py::arg("y"),
          R"pbdoc(
          SetAttitude
          ===========
          
          Set the attitude of the navigator
          
          Parameters
          ----------
          
          r : double
          
              roll angle [rad]
          
          p : double
          
              pitch angle [rad]
              
          y : double
          
              yaw angle [rad]
          )pbdoc")
      .def(
          "Mechanize",
          &Strapdown::Mechanize,
          py::arg("w_ib_b"),
          py::arg("f_ib_b"),
          py::arg("dt"),
          R"pbdoc(
          Mechanize
          =========
          
          Integrate measured angular rates (delta thetas) and specific forces (delta velocities)
          
          Parameters
          ----------
          
          w_ib_b : np.ndarray
          
              Measured angular rates (delta thetas) in the body frame [rad/s]
          
          f_ib_b : np.ndarray
          
              Measured specific forces (delta velocities) in the body frame [m/s^2]
              
          dt : double
          
              Integration time [s]
          )pbdoc")
      .def_readwrite("phi_", &Strapdown::phi_)
      .def_readwrite("lam_", &Strapdown::lam_)
      .def_readwrite("h_", &Strapdown::h_)
      .def_readwrite("vn_", &Strapdown::vn_)
      .def_readwrite("ve_", &Strapdown::ve_)
      .def_readwrite("vd_", &Strapdown::vd_)
      .def_readwrite("q_b_l_", &Strapdown::q_b_l_)
      .def_readwrite("C_b_l_", &Strapdown::C_b_l_)
      .doc() = R"pbdoc(
               Strapdown
               ========= 

               Inertial navigation strapdown integration equations.
               )pbdoc";

  // InertialNav
  py::class_<InertialNav>(h, "InertialNav")
      .def(py::init<>())
      .def(
          py::init<
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double>(),
          py::arg("lat"),
          py::arg("lon"),
          py::arg("alt"),
          py::arg("vn"),
          py::arg("ve"),
          py::arg("vd"),
          py::arg("r"),
          py::arg("p"),
          py::arg("y"),
          py::arg("cb"),
          py::arg("cd"))
      .def(
          "SetImuSpec",
          &InertialNav::SetImuSpec,
          py::arg("Ba"),
          py::arg("Na"),
          py::arg("Bg"),
          py::arg("Ng"),
          R"pbdoc(
          SetImuSpec
          ==========
          
          Set the noise parameters of the IMU
          
          Parameters
          ----------
          
          Ba : double
          
              Accelerometer bias InertialNavtability [mg]
          
          Na : double
          
              Accelerometer random walk [m/s/sqrt(hr)]

          Bg : double
          
              Gyroscope bias InertialNavtability [deg/hr]

          Ng : double
          
              Gyroscope random walk [deg/sqrt(hr)]
          )pbdoc")
      .def(
          "SetClockSpec",
          &InertialNav::SetClockSpec,
          py::arg("h0"),
          py::arg("h1"),
          py::arg("h2"),
          R"pbdoc(
          SetClockSpec
          ============
          
          Set the noise parameters of the Clock
          
          Parameters
          ----------
          
          h0 : double
          
              white frequency modulation
          
          h1 : double
          
              flicker frequency modulation

          h2 : double
          
              random walk frequency modulation
          )pbdoc")
      .def(
          "SetClock",
          &InertialNav::SetClock,
          py::arg("cb"),
          py::arg("cd"),
          R"pbdoc(
          SetClock
          ========
          
          Set the clock states of the navigator
          
          Parameters
          ----------
          
          cb : double
          
              Clock bias [m]
          
          cd : double
          
              Clock drift [m/s]
          )pbdoc")
      .def(
          "SetPosition",
          &InertialNav::SetPosition,
          py::arg("lat"),
          py::arg("lon"),
          py::arg("alt"),
          R"pbdoc(
          SetPosition
          ===========
          
          Set the position of the navigator
          
          Parameters
          ----------
          
          lat : double
          
              Latitude [rad]
          
          lon : double
          
              Longitude [rad]
              
          alt : double
          
              Altitude/Height [m]
          )pbdoc")
      .def(
          "SetVelocity",
          &InertialNav::SetVelocity,
          py::arg("vn"),
          py::arg("ve"),
          py::arg("vd"),
          R"pbdoc(
          SetVelocity
          ===========
          
          Set the velocity of the navigator
          
          Parameters
          ----------
          
          vn : double
          
              north velocity [m/s]
          
          ve : double
          
              east velocity [m/s]
              
          vd : double
          
              down velocity [m/s]
          )pbdoc")
      .def(
          "SetAttitude",
          py::overload_cast<const Eigen::Ref<const Eigen::Matrix3d> &>(&InertialNav::SetAttitude),
          py::arg("C"),
          R"pbdoc(
          SetAttitude
          ===========
          
          Set the attitude of the navigator
          
          Parameters
          ----------
          
          C : np.ndarray
          
              body-to-ned rotation matrix
          )pbdoc")
      .def(
          "SetAttitude",
          py::overload_cast<const double &, const double &, const double &>(
              &InertialNav::SetAttitude),
          py::arg("r"),
          py::arg("p"),
          py::arg("y"),
          R"pbdoc(
          SetAttitude
          ===========
          
          Set the attitude of the navigator
          
          Parameters
          ----------
          
          r : double
          
              roll angle [rad]
          
          p : double
          
              pitch angle [rad]
              
          y : double
          
              yaw angle [rad]
          )pbdoc")
      .def(
          "Mechanize",
          &InertialNav::Mechanize,
          py::arg("w_ib_b"),
          py::arg("f_ib_b"),
          py::arg("dt"),
          R"pbdoc(
          Mechanize
          =========
          
          Integrate measured angular rates (delta thetas) and specific forces (delta velocities)
          
          Parameters
          ----------
          
          w_ib_b : np.ndarray
          
              Measured angular rates (delta thetas) in the body frame [rad/s]
          
          f_ib_b : np.ndarray
          
              Measured specific forces (delta velocities) in the body frame [m/s^2]
              
          dt : double
          
              Integration time [s]
          )pbdoc")
      .def(
          "Propagate",
          &InertialNav::Propagate,
          py::arg("w_ib_b"),
          py::arg("f_ib_b"),
          py::arg("dt"),
          R"pbdoc(
          Propagate
          =========
          
          Propagate the error state matrices
          
          Parameters
          ----------
          
          w_ib_b : np.ndarray
          
              Measured angular rates (delta thetas) in the body frame [rad/s]
          
          f_ib_b : np.ndarray
          
              Measured specific forces (delta velocities) in the body frame [m/s^2]
              
          dt : double
          
              Integration time [s]
          )pbdoc")
      .def(
          "GnssUpdate",
          &InertialNav::GnssUpdate,
          py::arg("sv_pos"),
          py::arg("sv_vel"),
          py::arg("psr"),
          py::arg("psrdot"),
          py::arg("psrvar"),
          py::arg("psrdotvar"),
          R"pbdoc(
          GnssUpdate
          ==========
          
          Correct state with GPS measurements
          
          Parameters
          ----------
          
          sv_pos : np.ndarray
          
              Satellite ECEF positions [m]
          
          sv_vel : np.ndarray
          
              Satellite ECEF velocities [m/s]
              
          psr : np.ndarray
          
              Pseudorange measurements [m]
                            
          psrdot : np.ndarray
          
              Pseudorange-rate measurements [m/s]
                            
          psrvar : np.ndarray
          
              Pseudorange measurement variance [m^2]
                            
          psrdotvar : np.ndarray
          
              Pseudorange-rate measurement variance [(m/s)^2]
          )pbdoc")
      .def(
          "PhasedArrayUpdate",
          &InertialNav::PhasedArrayUpdate,
          py::arg("sv_pos"),
          py::arg("sv_vel"),
          py::arg("psr"),
          py::arg("psrdot"),
          py::arg("phase"),
          py::arg("psrvar"),
          py::arg("psrdotvar"),
          py::arg("phasevar"),
          py::arg("ant_xyz"),
          py::arg("n_ant"),
          py::arg("lamb"),
          R"pbdoc(
          PhasedArrayUpdate
          =================
          
          Correct state with GPS measurements
          
          Parameters
          ----------
          
          sv_pos : np.ndarray
          
              Satellite ECEF positions [m]
          
          sv_vel : np.ndarray
          
              Satellite ECEF velocities [m/s]
              
          psr : np.ndarray
          
              Pseudorange measurements [m]
                            
          psrdot : np.ndarray
          
              Pseudorange-rate measurements [m/s]

          phase : np.ndarray
          
              Delta phase measurements [m/s]
                            
          psrvar : np.ndarray
          
              Pseudorange measurement variance [m^2]
                            
          psrdotvar : np.ndarray
          
              Pseudorange-rate measurement variance [(m/s)^2]

          phasevar : np.ndarray

              Variance of each phase measurement

          ant_xyz : np.ndarray

              Known antenna positions in the body frame

          n_ant : np.ndarray

              Known number of antennas in the array

          lamb : np.ndarray

              Wavelength for the signal of interest [m/rad]
          )pbdoc")
      .def_readwrite("phi_", &InertialNav::phi_)
      .def_readwrite("lam_", &InertialNav::lam_)
      .def_readwrite("h_", &InertialNav::h_)
      .def_readwrite("vn_", &InertialNav::vn_)
      .def_readwrite("ve_", &InertialNav::ve_)
      .def_readwrite("vd_", &InertialNav::vd_)
      .def_readwrite("q_b_l_", &InertialNav::q_b_l_)
      .def_readwrite("C_b_l_", &InertialNav::C_b_l_)
      .def_readwrite("bg_", &InertialNav::bg_)
      .def_readwrite("ba_", &InertialNav::ba_)
      .def_readwrite("cb_", &InertialNav::cb_)
      .def_readwrite("cd_", &InertialNav::cd_)
      .def_readwrite("P_", &InertialNav::P_)
      .doc() = R"pbdoc(
               InertialNav
               ===

               Inertial navigation Kalman Filter equations.
               )pbdoc";

  // KinematicNav
  py::class_<KinematicNav>(h, "KinematicNav")
      .def(py::init<>())
      .def(
          py::init<
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double>(),
          py::arg("lat"),
          py::arg("lon"),
          py::arg("alt"),
          py::arg("vn"),
          py::arg("ve"),
          py::arg("vd"),
          py::arg("cb"),
          py::arg("cd"))
      .def(
          py::init<
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double,
              const double>(),
          py::arg("lat"),
          py::arg("lon"),
          py::arg("alt"),
          py::arg("vn"),
          py::arg("ve"),
          py::arg("vd"),
          py::arg("r"),
          py::arg("p"),
          py::arg("y"),
          py::arg("cb"),
          py::arg("cd"))
      .def(
          "SetProcessNoise",
          &KinematicNav::SetProcessNoise,
          py::arg("Svel"),
          py::arg("Satt"),
          R"pbdoc(
          SetProcessNoise
          ===============
          
          Set the noise parameters of the KinematicNav (constant velocity) model
          
          Parameters
          ----------
          
          Svel : double
          
              PSD of expected acceleration white noise [(m/s^2)^2]
          
          Satt : double
          
              PSD of expected angular rate white noise [(rad/s)^2]
          )pbdoc")
      .def(
          "SetClockSpec",
          &KinematicNav::SetClockSpec,
          py::arg("h0"),
          py::arg("h1"),
          py::arg("h2"),
          R"pbdoc(
          SetClockSpec
          ============
          
          Set the noise parameters of the Clock
          
          Parameters
          ----------
          
          h0 : double
          
              white frequency modulation
          
          h1 : double
          
              flicker frequency modulation

          h2 : double
          
              random walk frequency modulation
          )pbdoc")
      .def(
          "SetClock",
          &KinematicNav::SetClock,
          py::arg("cb"),
          py::arg("cd"),
          R"pbdoc(
          SetClock
          ========
          
          Set the clock states of the navigator
          
          Parameters
          ----------
          
          cb : double
          
              Clock bias [m]
          
          cd : double
          
              Clock drift [m/s]
          )pbdoc")
      .def(
          "SetPosition",
          &KinematicNav::SetPosition,
          py::arg("lat"),
          py::arg("lon"),
          py::arg("alt"),
          R"pbdoc(
          SetPosition
          ===========
          
          Set the position of the navigator
          
          Parameters
          ----------
          
          lat : double
          
              Latitude [rad]
          
          lon : double
          
              Longitude [rad]
              
          alt : double
          
              Altitude/Height [m]
          )pbdoc")
      .def(
          "SetVelocity",
          &KinematicNav::SetVelocity,
          py::arg("vn"),
          py::arg("ve"),
          py::arg("vd"),
          R"pbdoc(
          SetVelocity
          ===========
          
          Set the velocity of the navigator
          
          Parameters
          ----------
          
          vn : double
          
              north velocity [m/s]
          
          ve : double
          
              east velocity [m/s]
              
          vd : double
          
              down velocity [m/s]
          )pbdoc")
      .def(
          "SetAttitude",
          py::overload_cast<const Eigen::Ref<const Eigen::Matrix3d> &>(&KinematicNav::SetAttitude),
          py::arg("C"),
          R"pbdoc(
          SetAttitude
          ===========
          
          Set the attitude of the navigator
          
          Parameters
          ----------
          
          C : np.ndarray
          
              body-to-ned rotation matrix
          )pbdoc")
      .def(
          "SetAttitude",
          py::overload_cast<const double &, const double &, const double &>(
              &KinematicNav::SetAttitude),
          py::arg("r"),
          py::arg("p"),
          py::arg("y"),
          R"pbdoc(
          SetAttitude
          ===========
          
          Set the attitude of the navigator
          
          Parameters
          ----------
          
          r : double
          
              roll angle [rad]
          
          p : double
          
              pitch angle [rad]
              
          y : double
          
              yaw angle [rad]
          )pbdoc")
      .def(
          "Propagate",
          &KinematicNav::Propagate,
          py::arg("dt"),
          R"pbdoc(
          Propagate
          =========
          
          Propagate the error state matrices
          
          Parameters
          ----------
              
          dt : double
          
              Integration time [s]
          )pbdoc")
      .def(
          "GnssUpdate",
          &KinematicNav::GnssUpdate,
          py::arg("sv_pos"),
          py::arg("sv_vel"),
          py::arg("psr"),
          py::arg("psrdot"),
          py::arg("psrvar"),
          py::arg("psrdotvar"),
          R"pbdoc(
          GnssUpdate
          ==========
          
          Correct state with GPS measurements
          
          Parameters
          ----------
          
          sv_pos : np.ndarray
          
              Satellite ECEF positions [m]
          
          sv_vel : np.ndarray
          
              Satellite ECEF velocities [m/s]
              
          psr : np.ndarray
          
              Pseudorange measurements [m]
                            
          psrdot : np.ndarray
          
              Pseudorange-rate measurements [m/s]
                            
          psrvar : np.ndarray
          
              Pseudorange measurement variance [m^2]
                            
          psrdotvar : np.ndarray
          
              Pseudorange-rate measurement variance [(m/s)^2]
          )pbdoc")
      .def(
          "PhasedArrayUpdate",
          &KinematicNav::PhasedArrayUpdate,
          py::arg("sv_pos"),
          py::arg("sv_vel"),
          py::arg("psr"),
          py::arg("psrdot"),
          py::arg("phase"),
          py::arg("psrvar"),
          py::arg("psrdotvar"),
          py::arg("phasevar"),
          py::arg("ant_xyz"),
          py::arg("n_ant"),
          py::arg("lamb"),
          R"pbdoc(
          PhasedArrayUpdate
          =================
          
          Correct state with GPS measurements
          
          Parameters
          ----------
          
          sv_pos : np.ndarray
          
              Satellite ECEF positions [m]
          
          sv_vel : np.ndarray
          
              Satellite ECEF velocities [m/s]
              
          psr : np.ndarray
          
              Pseudorange measurements [m]
                            
          psrdot : np.ndarray
          
              Pseudorange-rate measurements [m/s]

          phase : np.ndarray
          
              Delta phase measurements [m/s]
                            
          psrvar : np.ndarray
          
              Pseudorange measurement variance [m^2]
                            
          psrdotvar : np.ndarray
          
              Pseudorange-rate measurement variance [(m/s)^2]

          phasevar : np.ndarray

              Variance of each phase measurement

          ant_xyz : np.ndarray

              Known antenna positions in the body frame

          n_ant : np.ndarray

              Known number of antennas in the array

          lamb : np.ndarray

              Wavelength for the signal of interest [m/rad]
          )pbdoc")
      .def(
          "AttitudeUpdate",
          &KinematicNav::AttitudeUpdate,
          py::arg("C"),
          py::arg("R"),
          R"pbdoc(
            AttitudeUpdate
            ==============
            
            Update the navigator attitude with attitude measurement
            
            Parameters
            ----------
            
            C : np.ndarray
            
                body-to-ned rotation matrix

            R : np.ndarray
            
                DCM variance
            )pbdoc")
      .def_readwrite("phi_", &KinematicNav::phi_)
      .def_readwrite("lam_", &KinematicNav::lam_)
      .def_readwrite("h_", &KinematicNav::h_)
      .def_readwrite("vn_", &KinematicNav::vn_)
      .def_readwrite("ve_", &KinematicNav::ve_)
      .def_readwrite("vd_", &KinematicNav::vd_)
      .def_readwrite("q_b_l_", &KinematicNav::q_b_l_)
      .def_readwrite("C_b_l_", &KinematicNav::C_b_l_)
      .def_readwrite("cb_", &KinematicNav::cb_)
      .def_readwrite("cd_", &KinematicNav::cd_)
      .def_readwrite("P_", &KinematicNav::P_)
      .doc() = R"pbdoc(
               KinematicNav
               === 

               KinematicNav navigation Kalman Filter equations.
               )pbdoc";

  // Least Squares
  py::module_ ls = h.def_submodule("leastsquares", R"pbdoc(
      Least Squares
      =============

      Least squares algorithms for initializing the GNSS-InertialNav filter.)pbdoc");

  // RangeAndRate
  ls.def(
      "RangeAndRate",
      &RangeAndRate,
      py::arg("pos"),
      py::arg("vel"),
      py::arg("cb"),
      py::arg("cd"),
      py::arg("sv_pos"),
      py::arg("sv_vel"),
      py::arg("pred_u"),
      py::arg("pred_udot"),
      py::arg("pred_psr"),
      py::arg("pred_psrdot"),
      R"pbdoc(
      RangeAndRate
      ============

      Predicts a range and rate based on a satellite location and velocity

      Parameters
      ----------

      pos : np.ndarray

          3x1 User ECEF position [m]

      vel : np.ndarray

          3x1 User ECEF velocity [m/s]

      cb : double

          User clock bias [m]

      cd : double

          User clock drift [m/s]

      sv_pos : np.ndarray

          3x1 Satellite ECEF positions [m]

      sv_vel : np.ndarray

          3x1 Satellite ECEF velocities [m/s]

      pred_u : np.ndarray

          3x1 reference to unit vector to satellite

      pred_udot : np.ndarray

          3x1 reference to unit vector rate of change to satellite

      pred_psr : double

          Reference to pseudorange prediction

      pred_psrdot :double

          Reference to pseudorange-rate prediction
      )pbdoc");

  // GnssPvt
  ls.def(
      "GnssPVT",
      &GnssPVT,
      py::arg("x"),
      py::arg("P"),
      py::arg("sv_pos"),
      py::arg("sv_vel"),
      py::arg("psr"),
      py::arg("psrdot"),
      py::arg("psrvar"),
      py::arg("psrdotvar"),
      R"pbdoc(
      GnssPVT
      =======

      Least Squares solver for GNSS position, velocity, and timing terms

      Parameters
      ----------

      x : np.ndarray

          Initial state estimate

      P : np.ndarray

          Initial covariance estimate

      sv_pos : np.ndarray

          Satellite ECEF positions [m]

      sv_vel : np.ndarray

          Satellite ECEF velocities [m/s]

      psr : np.ndarray

          Pseudorange measurements [m]

      psrdot : np.ndarray

          Pseudorange-rate measurements [m/s]

      psrvar : np.ndarray

          Pseudorange measurement variance [m^2]

      psrdotvar : np.ndarray

          Pseudorange-rate measurement variance [(m/s)^2]

      Returns
      -------

      status : bool

          Convergence success
      )pbdoc");

  // PhasedArrayAttitude
  ls.def(
      "PhasedArrayAttitude",
      &PhasedArrayAttitude,
      py::arg("C_b_l"),
      py::arg("u_ned"),
      py::arg("phase"),
      py::arg("phasevar"),
      py::arg("ant_xyz"),
      py::arg("n_ant"),
      py::arg("lamb"),
      py::arg("thresh") = 1e-6,
      R"pbdoc(
      PhasedArrayAttitude
      ===================

      Iterative attitude estimate based on the known spatial phase of an antenna array

      Parameters
      ----------

      C_b_l : np.ndarray

          Initial estimate of the body to local-nav frame attitude dcm

      u_ned : np.ndarray

          3 x n_sv Ephemeris based unit vectors in the local-nav frame

      phase : np.ndarray

          n_ant x n_sv matrix of measured differential gnss phase values

      phasevar : np.ndarray

          Variance of each phase measurement

      ant_xyz : np.ndarray

          Known antenna positions in the body frame

      n_ant : np.ndarray

          Known number of antennas in the array

      lamb : np.ndarray

          Wavelength for the signal of interest [m/rad]

      thresh : np.ndarray

          Desired threshold of convergence

      Returns
      -------

      status : bool

          Convergence success
      )pbdoc");

  // Wahba
  ls.def(
      "Wahba",
      &Wahba,
      py::arg("C_b_l"),
      py::arg("u_body"),
      py::arg("u_ned"),
      py::arg("u_body_var"),
      R"pbdoc(
      Whaba
      =====

      Wahba's problem solver

      Parameters
      ----------

      C_l_b : np.ndarray

          Attitude DCM estimate (local-nav to body)

      u_body : np.ndarray

          Measured unit vectors in the body frame

      u_ned : np.ndarray

          Ephemeris based unit vectors in the local-nav frame

      u_body_var : np.ndarray

          Variance of the measured unit vectors
      )pbdoc");

  // MUSIC
  ls.def(
      "MUSIC",
      [](double &az_mean,
         double &el_mean,
         const Eigen::Ref<const Eigen::VectorXcd> &P,
         const Eigen::Ref<const Eigen::Matrix3Xd> &ant_xyz,
         const int &n_ant,
         const double &lambda,
         const double &thresh = 1e-4) {
        MUSIC(az_mean, el_mean, P, ant_xyz, n_ant, lambda, thresh);
        return std::pair<double, double>(az_mean, el_mean);
      },
      py::arg("az_mean"),
      py::arg("el_mean"),
      py::arg("P"),
      py::arg("ant_xyz"),
      py::arg("n_ant"),
      py::arg("lamb"),
      py::arg("thresh") = 1e-4,
      pybind11::return_value_policy::reference_internal,
      R"pbdoc(
      MUSIC
      ===================

      MUSIC estimator using Prompt correlators (I & Q)

      Parameters
      ----------

      az_mean : double

          Azimuth estimates [rad]

      el_mean : double

          Elevation estimates [rad]

      P : np.ndarray

          Measured prompt correlators

      ant_xyz : np.ndarray

          Known antenna positions in the body frame

      n_ant : np.ndarray

          Known number of antennas in the array

      lamb : np.ndarray

          Wavelength for the signal of interest [m/rad]

      thresh : np.ndarray

          Desired threshold of convergence

      Returns
      -------

      status : bool

          Convergence success
      )pbdoc");

  // Nav sensors
  py::module_ ns = h.def_submodule("navsense", R"pbdoc(
      Nav Sensors
      ===========

      Generators for navigation sensor error parameters.

      Currently Contains the following sensors:

      1. `NavigationIMU`
      2. `NavigationClock`
      )pbdoc");

  // NavigationIMU
  py::class_<NavigationIMU>(ns, "NavigationIMU")
      .def(py::init<>())
      .def_readwrite("Ba", &NavigationIMU::Ba)
      .def_readwrite("Ka", &NavigationIMU::Ka)
      .def_readwrite("Na", &NavigationIMU::Na)
      .def_readwrite("Ta", &NavigationIMU::Ta)
      .def_readwrite("Bg", &NavigationIMU::Bg)
      .def_readwrite("Kg", &NavigationIMU::Kg)
      .def_readwrite("Ng", &NavigationIMU::Ng)
      .def_readwrite("Tg", &NavigationIMU::Tg)
      .doc() = R"pbdoc(
               NavigationIMU
               =============

               Struct of navigation IMU Allan variance values
               )pbdoc";

  // NavigationClock
  py::class_<NavigationClock>(ns, "NavigationClock")
      .def(py::init<>())
      .def_readwrite("h0", &NavigationClock::h0)
      .def_readwrite("h1", &NavigationClock::h1)
      .def_readwrite("h2", &NavigationClock::h2)
      .doc() = R"pbdoc(
               NavigationClock
               =============== 

               Struct of navigation clock Allan variance values
               )pbdoc";

  ns.attr("LOW_QUALTIY_TCXO") = LOW_QUALITY_TCXO;
  ns.attr("HIGH_QUALITY_TCXO") = HIGH_QUALITY_TCXO;
  ns.attr("OCXO") = OCXO;
  ns.attr("RUBIDIUM") = RUBIDIUM;
  ns.attr("CESIUM") = CESIUM;
  ns.attr("TACTICAL") = TACTICAL;
  ns.attr("AUTOMOTIVE") = AUTOMOTIVE;
  ns.attr("CONSUMER") = CONSUMER;

  // GetNavImu
  ns.def(
      "GetNavImu",
      &GetNavImu,
      py::arg("imu_name"),
      R"pbdoc(
      GetNavImu
      =========

      Generator for the default navigation IMUs

      Parameters
      ----------

      imu_name : str

          Name of the desired default IMU

      Returns
      -------

      imu : NavigationIMU

          Allan variance parameters for the desired IMU
      )pbdoc");

  // GetNavClock
  ns.def(
      "GetNavClock",
      &GetNavClock,
      py::arg("clock_name"),
      R"pbdoc(
      GetNavClock
      ===========

      Generator for the default navigation clocks

      Parameters
      ----------

      clock_name : str

          Name of the desired default clock

      Returns
      -------

      imu : NavigationClock

          Allan variance parameters for the desired clock
      )pbdoc");

  // NavImuToSiUnits
  ns.def(
      "NavImuToSiUnits",
      &NavImuToSiUnits,
      py::arg("imu"),
      R"pbdoc(
      NavImuToSiUnits
      ===============

      Converts the NavigationIMU parameters to SI units of a navigation filter

      Parameters
      ----------

      imu : NavigationIMU

          Allan variance parameters for the IMU
      )pbdoc");
}