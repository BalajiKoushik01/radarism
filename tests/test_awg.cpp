#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <random>
#include <sstream>
#include "radar/awg_profile.hpp"
#include "radar/waveforms.hpp"

using namespace radar;
static constexpr double kTwoPi = 2.0 * std::numbers::pi;

TEST(Awg, FftResampleOfBandlimitedToneIsExact) {
    const std::size_t N = 64, up = 4;
    CBuffer x(N);
    for (std::size_t n = 0; n < N; ++n) {
        const double ph = kTwoPi * 5.0 * static_cast<double>(n) / N;
        x[n] = {static_cast<float>(std::cos(ph)), static_cast<float>(std::sin(ph))};
    }
    const CBuffer y = resample_fft(x, up);
    ASSERT_EQ(y.size(), N * up);
    for (std::size_t m = 0; m < y.size(); ++m) {
        const double ph = kTwoPi * 5.0 * static_cast<double>(m) / (N * up);
        ASSERT_NEAR(y[m].real(), std::cos(ph), 1e-4) << m;
        ASSERT_NEAR(y[m].imag(), std::sin(ph), 1e-4) << m;
    }
}

TEST(Awg, ResamplePreservesOriginalSamplesForEvenAndOddLengths) {
    std::mt19937 rng(3);
    std::normal_distribution<float> g(0.0f, 1.0f);
    for (std::size_t N : {100u, 99u}) {
        CBuffer x(N);
        for (auto& v : x) v = {g(rng), g(rng)};
        const CBuffer y = resample_fft(x, 5);
        for (std::size_t m = 0; m < N; ++m) {
            ASSERT_NEAR(y[m * 5].real(), x[m].real(), 1e-4) << N << " " << m;
            ASSERT_NEAR(y[m * 5].imag(), x[m].imag(), 1e-4) << N << " " << m;
        }
    }
}

TEST(Awg, RealIfIsLoopCoherentAndSnapsToWholeCycles) {
    CBuffer x(1000, cf32{1.0f, 0.0f});   // 1 ms PRI at 1 MHz
    AwgParams p;
    p.fs_awg = 100e6; p.mode = AwgMode::RealIF; p.if_hz = 10.0003e6; p.snap_if = true;
    p.int16_out = false; p.headroom = 1.0;
    const auto w = make_awg(x, 1e6, p);
    EXPECT_NEAR(w.if_used_hz, 10e6, 1e-3);
    ASSERT_EQ(w.samples, 100000u);
    ASSERT_EQ(w.channels, 1u);
    for (std::size_t m = 0; m < w.samples; m += 997)
        ASSERT_NEAR(w.f32[m], std::cos(kTwoPi * 10e6 * static_cast<double>(m) / 100e6), 1e-4) << m;
    const double cycles = w.if_used_hz * w.loop_period_s;
    EXPECT_NEAR(cycles, std::round(cycles), 1e-9);
}

TEST(Awg, ValidationMessages) {
    CBuffer x(1000, cf32{1.0f, 0.0f});
    AwgParams p;
    p.fs_awg = 1.5e6;                                   // not an integer multiple of 1 MHz
    EXPECT_THROW(make_awg(x, 1e6, p), std::invalid_argument);
    p.fs_awg = 100e6; p.mode = AwgMode::RealIF; p.if_hz = 0.2e6;
    EXPECT_THROW(make_awg(x, 1e6, p), std::invalid_argument);   // IF below half the baseband rate
    p.if_hz = 49e6;
    EXPECT_THROW(make_awg(x, 1e6, p), std::invalid_argument);   // too close to the AWG Nyquist
    p.mode = AwgMode::IQ; p.dac_bits = 4;
    EXPECT_THROW(make_awg(x, 1e6, p), std::invalid_argument);
    p.dac_bits = 16; p.headroom = 1.5;
    EXPECT_THROW(make_awg(x, 1e6, p), std::invalid_argument);
}

TEST(Awg, GranularityAndMinimumLengthPadWithZeros) {
    CBuffer x(1000, cf32{1.0f, 0.0f});
    AwgParams p;
    p.fs_awg = 100e6; p.granularity = 64; p.int16_out = false;
    auto w = make_awg(x, 1e6, p);
    EXPECT_EQ(w.samples, 100032u);
    EXPECT_EQ(w.padded, 32u);
    EXPECT_EQ(w.channels, 2u);
    ASSERT_EQ(w.f32.size(), 2u * 100032u);
    EXPECT_EQ(w.f32[2 * (w.samples - 1)], 0.0f);
    EXPECT_NEAR(w.loop_period_s, 100032.0 / 100e6, 1e-15);
    p.min_samples = 200000;
    w = make_awg(x, 1e6, p);
    EXPECT_EQ(w.samples, 200000u);
    EXPECT_EQ(w.padded, 100000u);
}

TEST(Awg, DacQuantisationAndJustification) {
    CBuffer x(100);
    for (std::size_t n = 0; n < x.size(); ++n)
        x[n] = {static_cast<float>(std::sin(0.3 * n)), static_cast<float>(std::cos(0.2 * n))};
    AwgParams p;
    p.fs_awg = 1e6; p.dac_bits = 12; p.headroom = 0.9;
    const auto w = make_awg(x, 1e6, p);
    int mx = 0;
    for (auto v : w.i16) mx = std::max(mx, std::abs(static_cast<int>(v)));
    EXPECT_NEAR(mx, 0.9 * 2047.0, 1.0);
    for (std::size_t n = 0; n < x.size(); ++n)
        ASSERT_NEAR(w.i16[2 * n] / w.counts_per_unit, x[n].real(), 0.51 / w.counts_per_unit) << n;
    p.left_justify = true;
    const auto wl = make_awg(x, 1e6, p);
    int mxl = 0;
    for (auto v : wl.i16) {
        EXPECT_EQ(static_cast<int>(v) % 16, 0);
        mxl = std::max(mxl, std::abs(static_cast<int>(v)));
    }
    EXPECT_NEAR(mxl, 0.9 * 2047.0 * 16.0, 16.0);
    p.int16_out = false; p.headroom = 0.7;
    const auto wf = make_awg(x, 1e6, p);
    float mf = 0.0f;
    for (float v : wf.f32) mf = std::max(mf, std::fabs(v));
    EXPECT_NEAR(mf, 0.7f, 1e-5f);
}

TEST(Awg, LfmWaveformEndToEnd) {
    WaveformParams wp;
    wp.type = WaveformType::LFM; wp.fs = 1e6; wp.pri = 1e-3; wp.tau = 100e-6; wp.bandwidth = 2e5;
    const Waveform tx = generate(wp);
    AwgParams p;
    p.fs_awg = 4e6; p.headroom = 0.95;
    const auto w = make_awg(tx.samples, wp.fs, p);
    ASSERT_EQ(w.samples, 4000u);
    int mx = 0;
    for (auto v : w.i16) mx = std::max(mx, std::abs(static_cast<int>(v)));
    EXPECT_NEAR(mx, 0.95 * 32767.0, 1.0);
    // far from the pulse the interpolated listening window stays small (band-limited ringing only)
    for (std::size_t m = 600; m < 3800; ++m) {
        ASSERT_LT(std::abs(static_cast<int>(w.i16[2 * m])), 0.03 * 32767.0) << m;
        ASSERT_LT(std::abs(static_cast<int>(w.i16[2 * m + 1])), 0.03 * 32767.0) << m;
    }
}

TEST(Awg, WritesBinAndSidecar) {
    CBuffer x(200, cf32{0.5f, -0.25f});
    AwgParams p;
    p.fs_awg = 4e6;
    const auto w = make_awg(x, 1e6, p);
    const auto path = (std::filesystem::temp_directory_path() / "radarsim_awg_test.bin").string();
    write_awg(path, w, "tx_clean_for_awg");
    EXPECT_EQ(std::filesystem::file_size(path), w.samples * w.channels * 2);
    std::ifstream f(path + ".json");
    std::stringstream ss;
    ss << f.rdbuf();
    const std::string js = ss.str();
    EXPECT_NE(js.find("\"sample_rate_hz\""), std::string::npos);
    EXPECT_NE(js.find("\"loop_period_s\""), std::string::npos);
    EXPECT_NE(js.find("\"format\": \"int16\""), std::string::npos);
}
