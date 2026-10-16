#include <gtest/gtest.h>
#include <cmath>
#include <numbers>
#include "radar/channel.hpp"
#include "radar/waveforms.hpp"

using namespace radar;

static WaveformParams wpar(WaveformType t = WaveformType::Pulse) {
    WaveformParams p; p.type = t; p.fs = 1e6; p.pri = 10e-3; p.tau = 2e-3;  // N = 10000, pulse = 2000
    return p;
}
static double range_for_delay_samples(double d, double fs) { return kSpeedOfLight * (d / fs) / 2.0; }

TEST(Channel, IntegerDelayPlacesPulseExactlyWithRightAmplitude) {
    auto wp = wpar();
    ChannelParams c; c.snr_db = 20.0; c.range_m = range_for_delay_samples(3000.0, wp.fs);
    auto tx = generate(wp);
    auto rep = analyze(wp, c);
    auto echo = make_echo(tx, wp, c, rep);
    const double amp = std::sqrt(rep.pr_w);
    for (std::size_t n = 0; n < 2990; ++n) ASSERT_LT(std::abs(echo[n]), 1e-3 * amp) << n;
    for (std::size_t n = 3005; n < 4995; ++n) ASSERT_NEAR(std::abs(echo[n]), amp, 1e-3 * amp) << n;
    for (std::size_t n = 5010; n < echo.size(); ++n) ASSERT_LT(std::abs(echo[n]), 1e-3 * amp) << n;
}

TEST(Channel, FractionalDelayMatchesAnalyticChirp) {
    auto wp = wpar(WaveformType::LFM); wp.bandwidth = 1e5;
    ChannelParams c; c.fc = 10e9; c.snr_db = 20.0;
    c.range_m = range_for_delay_samples(3000.37, wp.fs);
    auto tx = generate(wp);
    auto rep = analyze(wp, c);
    auto echo = make_echo(tx, wp, c, rep);
    const double amp = std::sqrt(rep.pr_w);
    const double k = wp.bandwidth / wp.tau;
    const double carrier = -2.0 * std::numbers::pi * std::fmod(c.fc * rep.td_s, 1.0);
    const std::size_t lo = static_cast<std::size_t>(std::ceil(rep.delay_samples)) + 20;
    const std::size_t hi = static_cast<std::size_t>(std::floor(rep.delay_samples)) + tx.n_pulse - 20;
    for (std::size_t n = lo; n < hi; ++n) {
        const double t = (static_cast<double>(n) - rep.delay_samples) / wp.fs - 0.5 * wp.tau;
        const double ph = std::numbers::pi * k * t * t + carrier;
        const cf32 expected{static_cast<float>(amp * std::cos(ph)), static_cast<float>(amp * std::sin(ph))};
        ASSERT_LE(std::abs(echo[n] - expected), 2e-3 * amp) << "n=" << n;
    }
}

TEST(Channel, DopplerShiftsEchoFrequency) {
    auto wp = wpar();
    ChannelParams c; c.fc = 10e9; c.snr_db = 20.0;
    c.range_m = range_for_delay_samples(3000.0, wp.fs);
    c.velocity_mps = 50e3 * kSpeedOfLight / (2.0 * c.fc);   // fd = +50 kHz
    auto tx = generate(wp);
    auto rep = analyze(wp, c);
    EXPECT_NEAR(rep.doppler_hz, 50e3, 1e-6);
    auto echo = make_echo(tx, wp, c, rep);
    for (std::size_t n : {3500u, 4000u, 4500u}) {
        const double f = std::arg(echo[n + 1] * std::conj(echo[n])) * wp.fs / (2.0 * std::numbers::pi);
        EXPECT_NEAR(f, 50e3, 20.0) << n;
    }
}

TEST(Noise, PowerMatchesKT0FsF) {
    const double fs = 1e6;
    const double pn = thermal_noise_power(fs, 3.0);
    EXPECT_NEAR(pn / 7.98879e-15, 1.0, 1e-3);
    CBuffer x(1'000'000, cf32{0.0f, 0.0f});
    add_awgn(x, pn, 99);
    double sum_p = 0, sum_i2 = 0, sum_q2 = 0, sum_i = 0;
    for (const auto& s : x) {
        const double re = s.real(), im = s.imag();
        sum_p += re * re + im * im;
        sum_i2 += re * re;
        sum_q2 += im * im;
        sum_i += re;
    }
    const double n = static_cast<double>(x.size());
    EXPECT_NEAR(sum_p / n / pn, 1.0, 0.01);
    EXPECT_NEAR(sum_i2 / n / (pn / 2.0), 1.0, 0.015);
    EXPECT_NEAR(sum_q2 / n / (pn / 2.0), 1.0, 0.015);
    EXPECT_LT(std::fabs(sum_i / n), 0.01 * std::sqrt(pn / 2.0));
}

TEST(Radar, EquationNumericAndRangeFourthPower) {
    ChannelParams c; c.fc = 10e9; c.range_m = 10e3; c.snr_mode = SnrMode::Physical;
    c.pt_w = 1e3; c.gt_db = 30; c.gr_db = 30; c.rcs_m2 = 1.0; c.loss_db = 0.0;
    const double p1 = radar_equation_power(c);
    EXPECT_NEAR(p1 / 4.5291e-14, 1.0, 1e-3);
    c.range_m = 20e3;
    const double p2 = radar_equation_power(c);
    EXPECT_NEAR(p1 / p2, 16.0, 1e-9);
}

TEST(Radar, SnrModesAreConsistent) {
    auto wp = wpar();
    ChannelParams c; c.snr_mode = SnrMode::Manual; c.snr_db = 13.0;
    auto r = analyze(wp, c);
    EXPECT_NEAR(r.snr_in_db, 13.0, 1e-9);
    EXPECT_NEAR(r.pr_w / r.pn_w, std::pow(10.0, 1.3), 1e-9 * std::pow(10.0, 1.3));
    // Physical: doubling range costs 40*log10(2) dB of SNR
    c.snr_mode = SnrMode::Physical; c.range_m = 10e3;
    const double s1 = analyze(wp, c).snr_in_db;
    c.range_m = 20e3;
    const double s2 = analyze(wp, c).snr_in_db;
    EXPECT_NEAR(s1 - s2, 40.0 * std::log10(2.0), 1e-9);
}

TEST(Geometry, EclipsingBlindAndSecondTimeAroundFlags) {
    WaveformParams wp; wp.type = WaveformType::Pulse; wp.fs = 1e6; wp.pri = 1e-3; wp.tau = 1e-4;
    ChannelParams c;
    c.range_m = 100e3;  auto a = analyze(wp, c);
    EXPECT_FALSE(a.in_blind_zone); EXPECT_FALSE(a.straddles_pri); EXPECT_FALSE(a.second_time_around);
    c.range_m = 140e3;  a = analyze(wp, c);
    EXPECT_TRUE(a.straddles_pri); EXPECT_FALSE(a.second_time_around);
    c.range_m = 200e3;  a = analyze(wp, c);
    EXPECT_TRUE(a.second_time_around); EXPECT_FALSE(a.in_blind_zone); EXPECT_FALSE(a.straddles_pri);
    c.range_m = 10.0;   a = analyze(wp, c);
    EXPECT_TRUE(a.in_blind_zone);
    EXPECT_NEAR(a.r_unambiguous_m, kSpeedOfLight * 1e-3 / 2.0, 1e-6);
    EXPECT_NEAR(a.r_blind_m, kSpeedOfLight * 1e-4 / 2.0, 1e-6);
}

TEST(Channel, SecondTimeAroundEchoWrapsToCorrectPosition) {
    WaveformParams wp; wp.type = WaveformType::Pulse; wp.fs = 1e6; wp.pri = 1e-3; wp.tau = 1e-4;  // N=1000, pulse=100
    ChannelParams c; c.snr_db = 20.0; c.range_m = 200e3;  // delay ~1334.2 samples -> appears at ~334.2
    auto tx = generate(wp);
    auto rep = analyze(wp, c);
    ASSERT_TRUE(rep.second_time_around);
    auto echo = make_echo(tx, wp, c, rep);
    const double amp = std::sqrt(rep.pr_w);
    for (std::size_t n = 0; n < 318; ++n) ASSERT_LT(std::abs(echo[n]), 1e-3 * amp) << n;
    for (std::size_t n = 350; n < 420; ++n) ASSERT_NEAR(std::abs(echo[n]), amp, 2e-3 * amp) << n;
    for (std::size_t n = 452; n < echo.size(); ++n) ASSERT_LT(std::abs(echo[n]), 1e-3 * amp) << n;
}

TEST(Channel, EndToEndPowersMatchSnrDefinition) {
    auto wp = wpar();
    ChannelParams c; c.snr_db = 20.0; c.range_m = range_for_delay_samples(3000.5, wp.fs);
    auto tx = generate(wp);
    auto res = simulate_rx(tx, wp, c);
    ASSERT_EQ(res.rx.size(), tx.samples.size());
    auto mean_pow = [&](std::size_t a, std::size_t b) {
        double s = 0; for (std::size_t n = a; n < b; ++n) s += std::norm(res.rx[n]);
        return s / static_cast<double>(b - a);
    };
    const double pn = res.report.pn_w;
    EXPECT_NEAR(mean_pow(3050, 4950) / pn, 101.0, 101.0 * 0.02);   // Pr + Pn
    EXPECT_NEAR(mean_pow(6000, 10000) / pn, 1.0, 0.06);            // noise only
    EXPECT_NEAR(res.report.mf_gain_db, 10.0 * std::log10(2000.0), 1e-9);
}

TEST(Channel, ValidationCatchesBadInputs) {
    auto wp = wpar();
    ChannelParams c;
    EXPECT_TRUE(validate_channel(wp, c).empty());
    c.range_m = -1.0;
    EXPECT_FALSE(validate_channel(wp, c).empty());
    c.range_m = 10e3; c.velocity_mps = 1e5;   // fd ~ 6.7 MHz > 0.45 Fs
    EXPECT_FALSE(validate_channel(wp, c).empty());
    EXPECT_THROW(simulate_rx(generate(wp), wp, c), std::invalid_argument);
}
