/**
 * *least-squares.hpp*
 *
 * =======  ========================================================================================
 * @file    sturdins/least-squares.hpp
 * @brief   Least squares algorithms for initializing the GNSS-INS filter.
 * @date    January 2025
 * @author  Daniel Sturdivant <sturdivant20@gmail.com>
 * @ref     1. "Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems", 2nd
 *              Edition, 2013 - Groves
 * =======  ========================================================================================
 */

// TODO: gradient descent, levenburg-marquardt
// TODO: maybe a particle filter?

#ifndef STURDINS_LEAST_SQUARES_HPP
#define STURDINS_LEAST_SQUARES_HPP

#include <Eigen/Dense>

namespace sturdins {

/**
 * *=== RangeAndRate ===*
 * @brief predicts a range and rate based on a satellite location and velocity
 * @param pos         3x1 User ECEF position [m]
 * @param vel         3x1 User ECEF velocity [m/s]
 * @param cb          User clock bias [m]
 * @param cd          User clock drift [m/s]
 * @param sv_pos      3xN Satellite ECEF positions [m]
 * @param sv_vel      3xN Satellite ECEF velocities [m/s]
 * @param pred_u      3x1 reference to unit vector to satellite
 * @param pred_udot   3x1 reference to unit vector rate of change to satellite
 * @param pred_psr    Reference to pseudorange prediction
 * @param pred_psrdot Reference to pseudorange-rate prediction
 */
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
    double &pred_psrdot);

/**
 * *=== GnssPVT ===*
 * @brief Least Squares solver for GNSS position, velocity, and timing terms
 * @param x           Initial state estimate
 * @param P           Initial covariance estimate
 * @param sv_pos      Satellite ECEF positions [m]
 * @param sv_vel      Satellite ECEF velocities [m/s]
 * @param psr         Pseudorange measurements [m]
 * @param psrdot      Pseudorange-rate measurements [m/s]
 * @param psr_var     Pseudorange measurement variance [m^2]
 * @param psrdot_var  Pseudorange-rate measurement variance [(m/s)^2]
 */
bool GnssPVT(
    Eigen::Ref<Eigen::VectorXd> x,
    Eigen::Ref<Eigen::MatrixXd> P,
    const Eigen::Ref<const Eigen::MatrixXd> &sv_pos,
    const Eigen::Ref<const Eigen::MatrixXd> &sv_vel,
    const Eigen::Ref<const Eigen::VectorXd> &psr,
    const Eigen::Ref<const Eigen::VectorXd> &psrdot,
    const Eigen::Ref<const Eigen::VectorXd> &psr_var,
    const Eigen::Ref<const Eigen::VectorXd> &psrdot_var);

/**
 * *=== Wahba ===*
 * @brief Solves Wahba's problem using least squares
 * @param C_b_l      Attitude DCM estimate (body to local-nav)
 * @param u_body     Measured unit vectors in the body frame
 * @param u_ned      Ephemeris based unit vectors in the local-nav frame
 * @param u_body_var Variance of the measured unit vectors
 */
void Wahba(
    Eigen::Ref<Eigen::Matrix3d> C_b_l,
    const Eigen::Ref<const Eigen::MatrixXd> &u_body,
    const Eigen::Ref<const Eigen::MatrixXd> &u_ned,
    const Eigen::Ref<const Eigen::VectorXd> &u_body_var);

/**
 * *=== MUSIC ===*
 * @brief MUSIC estimator using Prompt correlators (I & Q)
 * @param az        Azimuth estimates [rad]
 * @param el        Elevation estimates [rad]
 * @param P         Measured prompt correlators
 * @param thresh    Desired threshold of convergence
 */
void MUSIC(
    double &az_mean,
    double &el_mean,
    const Eigen::Ref<const Eigen::VectorXcd> &P,
    const Eigen::Ref<const Eigen::MatrixXd> &ant_xyz,
    const int &n_ant,
    const double &lambda,
    const double &thresh = 1e-4);

}  // namespace sturdins

#endif