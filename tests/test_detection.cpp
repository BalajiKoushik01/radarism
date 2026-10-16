#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include "radar/cfar.hpp"
#include "radar/channel.hpp"
#include "radar/matched_filter.hpp"
#include "radar/waveforms.hpp"

using namespace radar;

static WaveformParams lfm() {
    WaveformParams p; p.type = WaveformType::LFM;
    p.fs = 1e6; p.pri = 10e-3; p.tau = 2e-3; p.bandwidth = 2e5;   // N = 10000, pulse = 2000, TBP = 400
    return p;
}
static double range_for_delay_samples(double d, double fs) { return kSpeedOfLight * (d / fs) / 2.0; }

static std::size_t argmax_abs(const CBuffer& y, std::size_t a, std::size_t b) {
    std::size_t best = a;
    for (std::size_t i = a; i < b; ++i) if (std::norm(y[i]) > std::norm(y[best])) best = i;
    return best;
}

TEST(MatchedFilter, NoiseFreeIntegerDelayGivesExactPeakAndAmplitude) {
    auto wp = lfm();
    ChannelParams c; c.snr_db = 20.0; c.range_m = range_for_delay_samples(3000.0, wp.fs);
    auto tx = generate(wp);
    auto rep = analyze(wp, c);
    auto echo = make_echo(tx, wp, c, rep);
    auto y = matched_filter(echo, tx.samples);
    EXPECT_EQ(argmax_abs(y, 0, y.size()), 3000u);
    EXPECT_NEAR(std::abs(y[3000]), std::sqrt(rep.pr_w), 1e-3 * std::sqrt(rep.pr_w));
}

TEST(MatchedFilter, FractionalDelayPeakWithinOneSampleAndLowLoss) {
    auto wp = lfm();
    ChannelParams c; c.snr_db = 20.0; c.range_m = range_for_delay_samples(3000.4, wp.fs);
    auto tx = generate(wp);
    auto rep = analyze(wp, c);
    auto echo = make_echo(tx, wp, c, rep);
    auto y = matched_filter(echo, tx.samples);
    const std::size_t pk = argmax_abs(y, 0, y.size());
    EXPECT_LE(std::llabs(static_cast<long long>(pk) - 3000), 1);
    EXPECT_NEAR(std::abs(y[pk]), std::sqrt(rep.pr_w), 0.05 * std::sqrt(rep.pr_w));
}

TEST(MatchedFilter, NoiseFloorIsPnOverEnergyAndGainIsPulseLength) {
    auto wp = lfm();
    ChannelParams c; c.snr_db = 0.0; c.range_m = range_for_delay_samples(3000.0, wp.fs);
    auto tx = generate(wp);
    auto res = simulate_rx(tx, wp, c);
    auto y = matched_filter(res.rx, tx.samples);
    auto p = magnitude_squared(y);
    const double E = static_cast<double>(tx.n_pulse);              // unit-amplitude replica
    double floor_mean = 0; std::size_t cnt = 0;
    for (std::size_t i = 5500; i < 10000; ++i) { floor_mean += p[i]; ++cnt; }   // > pulse length from the peak
    floor_mean /= static_cast<double>(cnt);
    EXPECT_NEAR(floor_mean / (res.report.pn_w / E), 1.0, 0.15);
    const std::size_t pk = argmax_abs(y, 2990, 3010);
    const double gain_db = 10.0 * std::log10(static_cast<double>(p[pk]) / floor_mean);
    EXPECT_NEAR(gain_db, res.report.snr_in_db + res.report.mf_gain_db, 1.0);
}

TEST(Cfar, AlphaFormulaIsExactForTwoCells) {
    EXPECT_NEAR(cfar_alpha(2, 0.01), 18.0, 1e-9);   // 2 * (100^0.5 - 1)
}

TEST(Cfar, FalseAlarmRateMatchesDesignOnIidExponentialPower) {
    const std::size_t n = 200000;
    std::mt19937_64 rng(5);
    std::exponential_distribution<float> ex(1.0f);
    RBuffer pw(n);
    for (auto& v : pw) v = ex(rng);
    CfarParams p; p.guard_cells = 2; p.train_cells = 16; p.pfa = 1e-3;
    auto r = ca_cfar(pw, p);
    EXPECT_GT(r.cells_over_threshold, 130u);   // expect ~200
    EXPECT_LT(r.cells_over_threshold, 280u);
}

TEST(Cfar, DetectsTargetAtCorrectRangeEndToEnd) {
    auto wp = lfm();
    ChannelParams c; c.snr_db = -17.0; c.range_m = range_for_delay_samples(3000.0, wp.fs);  // below the noise per sample
    auto tx = generate(wp);
    auto res = simulate_rx(tx, wp, c);
    auto p = magnitude_squared(matched_filter(res.rx, tx.samples));
    CfarParams cp; cp.guard_cells = 20; cp.train_cells = 64; cp.pfa = 1e-6;
    auto r = ca_cfar(p, cp);
    bool found = false;
    for (auto d : r.detections)
        if (std::llabs(static_cast<long long>(d) - 3000) <= 2) found = true;
    EXPECT_TRUE(found);
    EXPECT_LE(r.detections.size(), 3u);
    if (found) {
        const double rng_m = bin_to_range_m(3000, wp.fs);
        EXPECT_NEAR(rng_m, c.range_m, kSpeedOfLight / (2.0 * wp.fs));
    }
}

TEST(MatchedFilter, RejectsLengthMismatch) {
    CBuffer a(100), b(50);
    EXPECT_THROW(matched_filter(a, b), std::invalid_argument);
}
