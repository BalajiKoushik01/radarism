#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
#include <numbers>
#include "radar/fft.hpp"

namespace { const int kEnv = (setenv("RADAR_FFT_THREADS", "4", 1), 0); }   // read once, before the first plan

TEST(FftThreads, LargeToneLandsInTheRightBinAndRoundTrips) {
    (void)kEnv;
    const std::size_t N = static_cast<std::size_t>(1) << 18, k = 12345;
    radar::CBuffer x(N), X(N), y(N);
    for (std::size_t n = 0; n < N; ++n) {
        const double ph = 2.0 * std::numbers::pi * static_cast<double>((k * n) % N) / static_cast<double>(N);
        x[n] = {static_cast<float>(std::cos(ph)), static_cast<float>(std::sin(ph))};
    }
    radar::Fft f(N, true), b(N, false);
    f.execute(x.data(), X.data());
    std::size_t pk = 0;
    for (std::size_t i = 1; i < N; ++i) if (std::norm(X[i]) > std::norm(X[pk])) pk = i;
    EXPECT_EQ(pk, k);
    EXPECT_NEAR(std::abs(X[k]), static_cast<double>(N), 0.01 * static_cast<double>(N));
    b.execute(X.data(), y.data());
    for (std::size_t i = 0; i < N; i += 997) {
        ASSERT_NEAR(y[i].real() / static_cast<float>(N), x[i].real(), 2e-3f) << i;
        ASSERT_NEAR(y[i].imag() / static_cast<float>(N), x[i].imag(), 2e-3f) << i;
    }
}
