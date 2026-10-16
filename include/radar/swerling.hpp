#pragma once
#include <random>
#include <cmath>

namespace radar {

enum class SwerlingType {
    Swerling0 = 0, // Constant RCS
    Swerling1 = 1, // Rayleigh amplitude, constant across CPI (slow)
    Swerling2 = 2, // Rayleigh amplitude, fluctuates pulse-to-pulse (fast)
    Swerling3 = 3, // Chi-square (4 DoF), constant across CPI (slow)
    Swerling4 = 4  // Chi-square (4 DoF), fluctuates pulse-to-pulse (fast)
};

class SwerlingGenerator {
public:
    explicit SwerlingGenerator(uint64_t seed = 1337) : rng_(seed) {}

    // Generates fluctuating RCS multiplier relative to mean RCS (sigma_0)
    double sample_rcs(SwerlingType type, double sigma_0) {
        if (type == SwerlingType::Swerling0 || sigma_0 <= 0.0) {
            return sigma_0;
        }

        if (type == SwerlingType::Swerling1 || type == SwerlingType::Swerling2) {
            // Exponential distribution with mean = sigma_0
            std::exponential_distribution<double> dist(1.0 / sigma_0);
            return dist(rng_);
        } else {
            // Chi-square 4 DoF: Gamma(shape=2, scale=sigma_0 / 2)
            std::gamma_distribution<double> dist(2.0, sigma_0 / 2.0);
            return dist(rng_);
        }
    }

private:
    std::mt19937_64 rng_;
};

} // namespace radar
