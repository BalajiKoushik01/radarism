#include <gtest/gtest.h>
#include <fftw3.h>
#include <volk/volk.h>
#include <cmath>
#include <cstdint>
#include <numbers>
#include "radar/types.hpp"

using namespace radar;

TEST(Memory, BufferIs64ByteAligned) {
    CBuffer b(1000);
    EXPECT_EQ(reinterpret_cast<std::uintptr_t>(b.data()) % 64, 0u);
}

TEST(Volk, ComplexMultiplyMatchesScalar) {
    const unsigned n = 1023;  // odd length on purpose: tests the SIMD tail
    CBuffer a(n), b(n), out(n);
    for (unsigned i = 0; i < n; ++i) {
        a[i] = {0.01f * i, -0.02f * i};
        b[i] = {1.5f, 0.001f * i};
    }
    volk_32fc_x2_multiply_32fc(out.data(), a.data(), b.data(), n);
    for (unsigned i = 0; i < n; ++i) {
        cf32 ref = a[i] * b[i];
        EXPECT_NEAR(out[i].real(), ref.real(), 1e-2f * (1.0f + std::abs(ref)));
        EXPECT_NEAR(out[i].imag(), ref.imag(), 1e-2f * (1.0f + std::abs(ref)));
    }
}

TEST(Fftw, ToneLandsInCorrectBin) {
    const int N = 4096, k = 100;
    CBuffer x(N), X(N);
    for (int n = 0; n < N; ++n) {
        double ph = 2.0 * std::numbers::pi * k * n / N;
        x[n] = {static_cast<float>(std::cos(ph)), static_cast<float>(std::sin(ph))};
    }
    fftwf_plan p = fftwf_plan_dft_1d(N, reinterpret_cast<fftwf_complex*>(x.data()),
                                     reinterpret_cast<fftwf_complex*>(X.data()),
                                     FFTW_FORWARD, FFTW_ESTIMATE);
    ASSERT_NE(p, nullptr);
    fftwf_execute(p);
    fftwf_destroy_plan(p);
    int peak = 0;
    for (int i = 1; i < N; ++i)
        if (std::abs(X[i]) > std::abs(X[peak])) peak = i;
    EXPECT_EQ(peak, k);
    EXPECT_NEAR(std::abs(X[k]), static_cast<float>(N), 1.0f);
}

TEST(Fftw, ForwardInverseRoundTrip) {
    const int N = 8192;
    CBuffer x(N), X(N), y(N);
    for (int n = 0; n < N; ++n) x[n] = {std::sin(0.01f * n), std::cos(0.037f * n)};
    fftwf_plan f = fftwf_plan_dft_1d(N, reinterpret_cast<fftwf_complex*>(x.data()),
                                     reinterpret_cast<fftwf_complex*>(X.data()),
                                     FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_plan b = fftwf_plan_dft_1d(N, reinterpret_cast<fftwf_complex*>(X.data()),
                                     reinterpret_cast<fftwf_complex*>(y.data()),
                                     FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(f);
    fftwf_execute(b);
    fftwf_destroy_plan(f);
    fftwf_destroy_plan(b);
    for (int n = 0; n < N; ++n) {
        EXPECT_NEAR(y[n].real() / N, x[n].real(), 1e-4f);
        EXPECT_NEAR(y[n].imag() / N, x[n].imag(), 1e-4f);
    }
}
