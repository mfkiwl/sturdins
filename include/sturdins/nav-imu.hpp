/**
 * *nav-imu.hpp*
 *
 * =======  ========================================================================================
 * @file    sturdins/nav-imu.hpp
 * @brief   Implementation and definitions for calculating IMU parameters.
 * @date    January 2025
 * @author  Daniel Sturdivant <sturdivant20@gmail.com>
 * =======  ========================================================================================
 */

#ifndef STURDINS_NAV_IMU_HPP
#define STURDINS_NAV_IMU_HPP

#include <string>

namespace sturdins {

/**
 * @brief Struct of navigation IMU Allan variance values
 */
struct NavigationIMU {
  double Ba;  // accel bias instability [mg]
  double Ka;  // accel rate random walk [(m/s)/(hr*sqrt(hr)]
  double Na;  // accel random walk [m/s/sqrt(hr)]
  double Ta;  // accel correlation times [s]
  double Bg;  // gyro bias instability [deg/hr]
  double Kg;  // gyro rate random walk [deg/(hr*sqrt(hr))]
  double Ng;  // gyro random walk [deg/sqrt(hr)]
  double Tg;  // gyro correlation times [s]
};

inline constexpr NavigationIMU TACTICAL =
    NavigationIMU{0.5, 0.0, 0.2942, 60.0, 0.35, 0.0, 0.102, 100.0};
inline constexpr NavigationIMU AUTOMOTIVE =
    NavigationIMU{1.0, 0.0, 0.5884, 100.0, 180.0, 0.0, 3.0, 300.0};
inline constexpr NavigationIMU CONSUMER =
    NavigationIMU{2.4, 0.0, 0.5884, 100.0, 360.0, 0.0, 3.0, 300.0};

/**
 * *=== GetNavImu ===*
 * @brief Generator for the default navigation IMUs
 * @param imu_name  Name of the desired default IMU
 * @return Allan variance parameters for the desired IMU
 */
NavigationIMU GetNavImu(std::string imu_name);

/**
 * *=== NavImuToSiUnits ===*
 * @brief Converts the NavigationIMU parameters to SI units of a navigation filter
 * @param imu Allan variance parameters for the IMU
 */
void NavImuToSiUnits(NavigationIMU &imu);

}  // namespace sturdins

#endif