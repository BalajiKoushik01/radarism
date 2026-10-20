#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <numbers>
#include <random>
#include "radar/fft.hpp"
#include "radar/simulation.hpp"

using namespace radar;

static CBuffer noise(std::size_t n, unsigned seed, double sigma) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<float> g(0.0f, static_cast<float>(sigma));
    CBuffer x(n);
    for (auto& v : x) v = {g(rng), g(rng)};
    return x;
}

TEST(Impairments, OffLeavesTheBufferUntouched) {
    CBuffer x = noise(1000, 1, 1.0), y = x;
    RxImpairments r;
    EXPECT_FALSE(r.active());
    EXPECT_FALSE(apply_rx_impairments(y, r, 1.0).applied);
    EXPECT_EQ(x, y);
}

TEST(Impairments, IqImbalanceImageMatchesTheory) {
    const std::size_t N = 4096, k = 300;
    CBuffer x(N), y, Y(N);
    for (std::size_t n = 0; n < N; ++n) {
        const double ph = 2.0 * std::numbers::pi * static_cast<double>(k * n) / static_cast<double>(N);
        x[n] = {static_cast<float>(std::cos(ph)), static_cast<float>(std::sin(ph))};
    }
    y = x;
    RxImpairments r;
    r.iq_gain_db = 1.0; r.iq_phase_deg = 5.0;
    const auto rep = apply_rx_impairments(y, r, 1.0);
    Fft f(N, true);
    f.execute(y.data(), Y.data());
    const double measured = 10.0 * std::log10(static_cast<double>(std::norm(Y[k])) / static_cast<double>(std::norm(Y[N - k])));
    const double g = std::pow(10.0, 1.0 / 20.0), phi = 5.0 * std::numbers::pi / 180.0;
    const std::complex<double> k1(0.5 * (1.0 + g * std::cos(phi)), 0.5 * g * std::sin(phi));
    const std::complex<double> k2(0.5 * (1.0 - g * std::cos(phi)), 0.5 * g * std::sin(phi));
    EXPECT_NEAR(measured, 20.0 * std::log10(std::abs(k1) / std::abs(k2)), 0.05);   // y = K1 x + K2 conj(x)
    EXPECT_NEAR(rep.image_rejection_db, measured, 0.05);
}

TEST(Impairments, DcOffsetShiftsTheMean) {
    const double sigma = 0.7;
    CBuffer x = noise(200000, 2, sigma), y = x;
    RxImpairments r;
    r.dc_i_sigma = 0.5; r.dc_q_sigma = -1.5;
    apply_rx_impairments(y, r, sigma);
    double di = 0.0, dq = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) { di += y[i].real() - x[i].real(); dq += y[i].imag() - x[i].imag(); }
    EXPECT_NEAR(di / x.size() / sigma, 0.5, 1e-4);
    EXPECT_NEAR(dq / x.size() / sigma, -1.5, 1e-4);
}

TEST(Impairments, AdcNoiseIsLsbSquaredOver12) {
    CBuffer x = noise(1'000'000, 3, 1.0), y = x;
    RxImpairments r;
    r.adc_bits = 8; r.adc_fs_sigma = 4.0;
    const auto rep = apply_rx_impairments(y, r, 1.0);
    const double lsb = 8.0 / 256.0;
    EXPECT_NEAR(rep.lsb, lsb, 1e-12);
    double acc = 0.0; std::size_t cnt = 0;
    for (std::size_t i = 0; i < x.size(); ++i)
        for (int c = 0; c < 2; ++c) {
            const double v = c ? x[i].imag() : x[i].real(), w = c ? y[i].imag() : y[i].real();
            if (std::fabs(v) < 3.9) { acc += (w - v) * (w - v); ++cnt; }
        }
    EXPECT_NEAR(acc / static_cast<double>(cnt) / (lsb * lsb / 12.0), 1.0, 0.03);
}

TEST(Impairments, ClippingFractionMatchesTheGaussianTail) {
    CBuffer x = noise(500000, 4, 1.0);
    RxImpairments r;
    r.adc_bits = 8; r.adc_fs_sigma = 1.0;
    const auto rep = apply_rx_impairments(x, r, 1.0);
    EXPECT_NEAR(static_cast<double>(rep.clipped) / (2.0 * 500000.0), 0.3173, 0.004);   // P(|N(0,1)| > 1)
}

TEST(Impairments, AutoFullScaleNeverClipsAndErrorIsHalfLsb) {
    CBuffer x = noise(100000, 5, 1.0), y = x;
    RxImpairments r;
    r.adc_bits = 8;
    const auto rep = apply_rx_impairments(y, r, 1.0);
    EXPECT_EQ(rep.clipped, 0u);
    for (std::size_t i = 0; i < x.size(); i += 7)
        ASSERT_LE(std::fabs(y[i].real() - x[i].real()), 0.5 * rep.lsb * 1.001);
}

TEST(Impairments, ValidationAndPipeline) {
    SimulationParams p;
    p.wf.type = WaveformType::LFM; p.wf.fs = 1e6; p.wf.pri = 10e-3; p.wf.tau = 500e-6; p.wf.bandwidth = 2e5;
    p.ch.snr_db = -10.0; p.ch.range_m = 150e3;
    const auto r0 = simulate(p);
    p.rx_imp.dc_i_sigma = 2.0;
    const auto r1 = simulate(p);
    const double sigma = std::sqrt(r0->report.pn_w / 2.0);
    EXPECT_FALSE(r0->imp.applied);
    EXPECT_TRUE(r1->imp.applied);
    for (std::size_t i : {10u, 4000u, 9000u}) {
        EXPECT_NEAR((r1->rx[i].real() - r0->rx[i].real()) / sigma, 2.0, 1e-3);
        EXPECT_NEAR((r1->rx[i].imag() - r0->rx[i].imag()) / sigma, 0.0, 1e-3);
    }
    p.rx_imp = RxImpairments{};
    p.rx_imp.adc_bits = 1;
    EXPECT_FALSE(validate_all(p).empty());
}
