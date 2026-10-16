#include <gtest/gtest.h>
#include <cmath>
#include <numbers>
#include <random>
#include "radar/fft.hpp"
#include "radar/phase_noise.hpp"
#include "radar/psd.hpp"
#include "radar/waveforms.hpp"

using namespace radar;

// Measured one-sided PSD of a phase sequence vs the target, averaged over a frequency band
static void band_compare(const RBuffer& phi, double fs, double f1, double f2,
                         const PhaseNoiseParams& p, double& measured, double& target) {
    const std::size_t n = phi.size();
    CBuffer x(n), X(n);
    for (std::size_t i = 0; i < n; ++i) x[i] = {phi[i], 0.0f};
    Fft fft(n, true);
    fft.execute(x.data(), X.data());
    const double df = fs / static_cast<double>(n);
    measured = 0.0; target = 0.0;
    int cnt = 0;
    for (std::size_t k = 1; k < n / 2; ++k) {
        const double f = k * df;
        if (f < f1 || f > f2) continue;
        measured += 2.0 * static_cast<double>(std::norm(X[k])) / (static_cast<double>(n) * fs);
        target += leeson_psd(p, f);
        ++cnt;
    }
    ASSERT_GT(cnt, 20);
    measured /= cnt; target /= cnt;
}

TEST(PhaseNoise, WhiteFloorVarianceMatchesTheory) {
    PhaseNoiseParams p; p.floor_dbc_hz = -100; p.leeson_corner_hz = 0; p.flicker_corner_hz = 0;
    const double fs = 1e6; const std::size_t n = 65536;
    RBuffer phi = generate_phase_noise(n, fs, p);
    double v = 0; for (float s : phi) v += static_cast<double>(s) * s;
    v /= n;
    const double expected = 2.0 * 1e-10 * (fs / 2.0);  // S_phi * Fs/2 = 1e-4 rad^2
    EXPECT_NEAR(v / expected, 1.0, 0.05);
}

TEST(PhaseNoise, SpectrumFollowsLeesonShape) {
    PhaseNoiseParams p; p.floor_dbc_hz = -130; p.leeson_corner_hz = 5e4; p.flicker_corner_hz = 0;
    const double fs = 1e6; const std::size_t n = 65536;
    RBuffer phi = generate_phase_noise(n, fs, p);
    double m, t;
    band_compare(phi, fs, 2e3, 4e3, p, m, t);      // steep 1/f^2 region
    EXPECT_NEAR(m / t, 1.0, 0.35);
    band_compare(phi, fs, 2e5, 3e5, p, m, t);      // flat floor region
    EXPECT_NEAR(m / t, 1.0, 0.10);
}

TEST(PhaseNoise, ApplyKeepsUnitModulusAndZeroPadding) {
    WaveformParams wp; wp.type = WaveformType::Pulse; wp.fs = 10e6; wp.pri = 1e-3; wp.tau = 100e-6;
    Waveform w = generate(wp);
    CBuffer orig = w.samples;
    PhaseNoiseParams p; p.floor_dbc_hz = -80; p.leeson_corner_hz = 1e4; p.flicker_corner_hz = 1e3;
    apply_phase_noise(w.samples, wp.fs, p);
    ASSERT_EQ(w.samples.size(), orig.size());
    for (std::size_t i = 0; i < w.n_pulse; ++i) ASSERT_NEAR(std::abs(w.samples[i]), 1.0f, 1e-4f);
    for (std::size_t i = w.n_pulse; i < w.samples.size(); ++i) ASSERT_EQ(w.samples[i], cf32(0.0f, 0.0f));
    EXPECT_GT(std::abs(w.samples[500] - orig[500]), 1e-6f);  // phase actually changed
}

TEST(Psd, ToneIntegratesToItsPowerAndPeaksAtRightFrequency) {
    const std::size_t n = 8192; const double fs = 1e6, f0 = 125e3, A = 0.5;  // f0 is exactly bin 1024
    CBuffer x(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double ph = 2.0 * std::numbers::pi * f0 * i / fs;
        x[i] = {static_cast<float>(A * std::cos(ph)), static_cast<float>(A * std::sin(ph))};
    }
    Psd s = compute_psd(x, fs, Window::Hann);
    double total = 0; std::size_t peak = 0;
    for (std::size_t i = 0; i < n; ++i) {
        total += static_cast<double>(s.density[i]) * s.df;
        if (s.density[i] > s.density[peak]) peak = i;
    }
    EXPECT_NEAR(total / (A * A), 1.0, 0.01);
    EXPECT_NEAR(s.f_start + peak * s.df, f0, s.df);
}

TEST(Psd, WhiteNoiseLevelIsSigmaSquaredOverFs) {
    const std::size_t n = 65536; const double fs = 1e6;
    std::mt19937_64 rng(7);
    std::normal_distribution<float> g(0.0f, std::sqrt(0.5f));  // complex variance = 1
    CBuffer x(n);
    for (auto& v : x) v = {g(rng), g(rng)};
    Psd s = compute_psd(x, fs, Window::Hann);
    double mean = 0; for (float d : s.density) mean += d;
    mean /= n;
    EXPECT_NEAR(mean * fs, 1.0, 0.03);
}

TEST(PulsePsd, EnergyConservedAndLfmPowerInsideBandwidth) {
    WaveformParams wp; wp.type = WaveformType::LFM; wp.fs = 1e6; wp.pri = 10e-3; wp.tau = 1e-3; wp.bandwidth = 2e5;
    Waveform w = generate(wp);   // 1000 active samples, time-bandwidth product 200
    for (Window win : {Window::Rectangular, Window::Hann}) {
        Psd s = compute_pulse_psd(w.samples, w.n_pulse, wp.fs, win);
        double total = 0, inband = 0;
        for (std::size_t i = 0; i < s.density.size(); ++i) {
            const double f = s.f_start + static_cast<double>(i) * s.df;
            const double p = static_cast<double>(s.density[i]) * s.df;
            total += p;
            if (std::fabs(f) <= 0.55 * wp.bandwidth) inband += p;
        }
        EXPECT_NEAR(total, 1.0, 1e-3);
        EXPECT_GT(inband / total, 0.9);
    }
}

TEST(PulsePsd, HannLowersFarSidelobesVersusRectangular) {
    WaveformParams wp; wp.type = WaveformType::Pulse; wp.fs = 1e6; wp.pri = 10e-3; wp.tau = 500e-6;
    Waveform w = generate(wp);
    auto mean_db = [&](Window win) {
        Psd s = compute_pulse_psd(w.samples, w.n_pulse, wp.fs, win);
        const double peak = *std::max_element(s.density.begin(), s.density.end());
        double acc = 0; int cnt = 0;
        for (std::size_t i = 0; i < s.density.size(); ++i) {
            const double f = std::fabs(s.f_start + static_cast<double>(i) * s.df);
            if (f > 3e5 && f < 4.5e5) { acc += 10.0 * std::log10(std::max(static_cast<double>(s.density[i]) / peak, 1e-18)); ++cnt; }
        }
        return acc / cnt;
    };
    EXPECT_LT(mean_db(Window::Hann), mean_db(Window::Rectangular) - 20.0);
}
