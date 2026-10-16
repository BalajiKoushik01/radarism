#include <gtest/gtest.h>
#include <cmath>
#include <numbers>
#include "radar/passband.hpp"

using namespace radar;
static constexpr double kTwoPi = 2.0 * std::numbers::pi;

TEST(Passband, ConstantIGivesCosineAndConstantQGivesMinusSine) {
    const double fs = 1e6, fc = 10e9, t0 = 500e-6, dt = 5e-12;
    CBuffer xi(1000, cf32{1.0f, 0.0f}), xq(1000, cf32{0.0f, 1.0f});
    RBuffer a, b;
    reconstruct_passband(xi, fs, fc, t0, dt, 2000, a);
    reconstruct_passband(xq, fs, fc, t0, dt, 2000, b);
    for (std::size_t k = 0; k < 2000; ++k) {
        const double t = t0 + k * dt;
        const double ph = kTwoPi * std::fmod(fc * t, 1.0);
        ASSERT_NEAR(a[k], std::cos(ph), 1e-4) << k;
        ASSERT_NEAR(b[k], -std::sin(ph), 1e-4) << k;
    }
}

TEST(Passband, BasebandToneAppearsAtCarrierPlusOffset) {
    const double fs = 1e6, fc = 1e9, f0 = 1e3;
    CBuffer x(2000);
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double ph = kTwoPi * f0 * static_cast<double>(i) / fs;
        x[i] = {static_cast<float>(std::cos(ph)), static_cast<float>(std::sin(ph))};
    }
    const double t0 = 1e-3, dt = 20e-12;
    RBuffer rf;
    reconstruct_passband(x, fs, fc, t0, dt, 500, rf);   // 10 ns = 10 carrier cycles
    for (std::size_t k = 0; k < 500; ++k) {
        const double t = t0 + k * dt;
        ASSERT_NEAR(rf[k], std::cos(kTwoPi * std::fmod((fc + f0) * t, 1.0)), 1e-4) << k;
    }
}

TEST(Passband, ZeroOutsideBufferAndAmplitudeEqualsMagnitude) {
    CBuffer x(100, cf32{1.0f, 0.0f});
    RBuffer out;
    reconstruct_passband(x, 1e6, 0.0, -5e-6, 1e-6, 5, out);   // fc = 0: RF = I
    EXPECT_EQ(out[0], 0.0f);                                  // t = -5 us: before the buffer
    reconstruct_passband(x, 1e6, 0.0, 50e-6, 1e-6, 5, out);
    EXPECT_NEAR(out[2], 1.0f, 1e-6f);
    reconstruct_passband(x, 1e6, 0.0, 200e-6, 1e-6, 5, out);  // beyond the buffer
    EXPECT_EQ(out[0], 0.0f);

    CBuffer y(100, cf32{0.6f, 0.8f});
    reconstruct_passband(y, 1e6, 1e9, 40e-6, 10e-12, 1000, out);   // 10 cycles, 100 points per cycle
    float peak = 0;
    for (float v : out) peak = std::max(peak, std::fabs(v));
    EXPECT_NEAR(peak, 1.0f, 5e-3f);
}

TEST(Passband, RejectsBadArguments) {
    CBuffer x(10);
    RBuffer out;
    EXPECT_THROW(reconstruct_passband(x, 0.0, 1e9, 0.0, 1e-12, 4, out), std::invalid_argument);
    EXPECT_THROW(reconstruct_passband(x, 1e6, 1e9, 0.0, 0.0, 4, out), std::invalid_argument);
    EXPECT_THROW(reconstruct_passband(CBuffer{}, 1e6, 1e9, 0.0, 1e-12, 4, out), std::invalid_argument);
}
