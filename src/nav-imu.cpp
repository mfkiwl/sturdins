/**
 * *nav-imu.cpp*
 *
 * =======  ========================================================================================
 * @file    sturdins/nav-imu.cpp
 * @brief   Implementation and definitions for calculating IMU parameters.
 * @date    January 2025
 * @author  Daniel Sturdivant <sturdivant20@gmail.com>
 * =======  ========================================================================================
 */

#include "sturdins/nav-imu.hpp"

#include <navtools/constants.hpp>

namespace sturdins {

// *=== GetNavImu ===*
NavigationIMU GetNavImu(std::string imu_name) {
  for (char &c : imu_name) {
    c = std::tolower(c);
  };
  if (imu_name == "tactical") {
    return TACTICAL;
  } else if (imu_name == "automotive") {
    return AUTOMOTIVE;
  } else if (imu_name == "consumer") {
    return CONSUMER;
  } else {
    // Default to automotive
    return AUTOMOTIVE;
  }
};

// *=== NavImuToSiUnits ===*
void NavImuToSiUnits(NavigationIMU &imu) {
  imu.Ba *= 9.80665 / 1000.0;                // [mg] -> [(m/s)/s]
  imu.Ka /= 216000.0;                        // [(m/s)/(hr*√hr] -> [(m/s)/(s*√s)]
  imu.Na /= 60.0;                            // [(m/s)/√hr] -> [(m/s)/√s]
  imu.Bg *= navtools::DEG2RAD<> / 3600.0;    // [deg/hr] ->  [rad/s]
  imu.Kg *= navtools::DEG2RAD<> / 216000.0;  // [deg/(hr*√hr)] ->  [rad/(s*√s)]
  imu.Ng *= navtools::DEG2RAD<> / 60.0;      // [deg/√hr]  ->  [rad/√s]
}

}  // namespace sturdins