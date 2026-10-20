#pragma once
#include <algorithm>
#include <cmath>
#include "radar/channel.hpp"
#include "radar/types.hpp"

namespace radar {

// Offset (bins, within [-0.5, 0.5]) of the true peak from bin idx: parabola through the log-power of the three
// samples around idx (circular indexing). Returns 0 when the neighbourhood is not a clean peak.
inline double peak_offset(const RBuffer& power, std::size_t idx) {
    const std::size_t n = power.size();
    if (n < 3 || idx >= n) return 0.0;
    const double a = std::log(std::max(static_cast<double>(power[(idx + n - 1) % n]), 1e-300));
    const double b = std::log(std::max(static_cast<double>(power[idx]), 1e-300));
    const double c = std::log(std::max(static_cast<double>(power[(idx + 1) % n]), 1e-300));
    const double denom = a - 2.0 * b + c;
    if (!(denom < 0.0)) return 0.0;
    return std::clamp(0.5 * (a - c) / denom, -0.5, 0.5);
}

inline double range_from_bin(double bin, double fs) { return kSpeedOfLight * bin / (2.0 * fs); }

// Parabolic peak offset (bins, within [-0.5, 0.5]) from three LINEAR power values y(-1), y(0), y(+1): the parabola
// is fitted through their logarithms. 0 when the middle value is not a clean local maximum.
inline double parabolic_offset(double a, double b, double c) {
    a = std::log(std::max(a, 1e-300));
    b = std::log(std::max(b, 1e-300));
    c = std::log(std::max(c, 1e-300));
    const double denom = a - 2.0 * b + c;
    if (!(denom < 0.0)) return 0.0;
    return std::clamp(0.5 * (a - c) / denom, -0.5, 0.5);
}

}  // namespace radar
