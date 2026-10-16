#pragma once
#include <vector>
#include <string>
#include "radar/swerling.hpp"

namespace radar {

struct Target {
    double range_m = 10000.0;        // Target range in meters
    double velocity_mps = 0.0;      // Radial velocity in m/s
    double rcs_m2 = 1.0;            // Mean Radar Cross Section (sigma_0)
    SwerlingType swerling = SwerlingType::Swerling0;
    std::string label = "Target";
};

using TargetScene = std::vector<Target>;

} // namespace radar
