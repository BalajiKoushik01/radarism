#include <gtest/gtest.h>
#include <cmath>
#include <initializer_list>
#include <random>
#include <utility>
#include "radar/peak_detect.hpp"
#include "radar/simulation.hpp"

using namespace radar;

static RBuffer profile(std::size_t n, unsigned seed, std::initializer_list<std::pair<double, double>> peaks) {
    std::mt19937 rng(seed);
    std::exponential_distribution<float> ex(1.0f);
    RBuffer p(n);
    for (auto& v : p) v = ex(rng);
    for (const auto& pr : peaks)
        for (std::size_t i = 0; i < n; ++i) {
            const double d = static_cast<double>(i) - pr.first;
            p[i] += static_cast<float>(pr.second * std::exp(-d * d / (2.0 * 1.5 * 1.5)));
        }
    return p;
}

TEST(PeakDetect, FindsPeaksStrongestFirstWithSubBinRange) {
    const double fs = 1e6;
    const RBuffer p = profile(10000, 1, {{1000.3, 200.0}, {4000.0, 100.0}, {7000.0, 60.0}});
    const auto pk = find_peaks(p, 1.0, 15.0, 5, 5, fs);
    ASSERT_GE(pk.size(), 3u);
    EXPECT_EQ(pk[0].bin, 1000u);
    EXPECT_EQ(pk[1].bin, 4000u);
    EXPECT_EQ(pk[2].bin, 7000u);
    EXPECT_NEAR(pk[0].bin_refined, 1000.3, 0.15);
    EXPECT_NEAR(pk[0].range_m, kSpeedOfLight * pk[0].bin_refined / (2.0 * fs), 1e-6);
    EXPECT_TRUE(pk[0].above_threshold);
    EXPECT_GT(pk[0].snr_db, pk[1].snr_db);
    EXPECT_GT(pk[1].snr_db, pk[2].snr_db);
}

TEST(PeakDetect, NonMaximumSuppressionHonoursTheSeparation) {
    const RBuffer p = profile(10000, 2, {{2000.0, 100.0}, {2006.0, 40.0}});
    EXPECT_EQ(find_peaks(p, 1.0, 15.0, 8, 5, 1e6).size(), 1u);
    EXPECT_EQ(find_peaks(p, 1.0, 15.0, 3, 5, 1e6).size(), 2u);
}

TEST(PeakDetect, TheStrongestPeakIsReportedEvenBelowTheThreshold) {
    const RBuffer p = profile(10000, 3, {});
    const auto pk = find_peaks(p, 1.0, 40.0, 5, 5, 1e6);
    ASSERT_EQ(pk.size(), 1u);
    EXPECT_FALSE(pk[0].above_threshold);
}

TEST(PeakDetect, SeparationFollowsTheCompressedPulseWidth) {
    WaveformParams w;
    w.fs = 1e6; w.tau = 130e-6; w.bandwidth = 2e5; w.code_length = 13;
    w.type = WaveformType::LFM;    EXPECT_EQ(suggest_peak_separation(w), 10u);
    w.type = WaveformType::Pulse;  EXPECT_EQ(suggest_peak_separation(w), 130u);
    w.type = WaveformType::Barker; EXPECT_EQ(suggest_peak_separation(w), 20u);
    w.type = WaveformType::CW;     EXPECT_EQ(suggest_peak_separation(w), 3u);
}

TEST(PeakDetect, EndToEndRangeOfASimulatedTarget) {
    SimulationParams p;
    p.wf.type = WaveformType::LFM; p.wf.fs = 1e6; p.wf.pri = 10e-3; p.wf.tau = 500e-6; p.wf.bandwidth = 2e5;
    p.ch.snr_db = -5.0;
    p.ch.range_m = kSpeedOfLight * (3000.4 / p.wf.fs) / 2.0;
    const auto r = simulate(p);
    const auto pk = find_peaks(r->mf_power, r->mf_noise_floor_w, 15.0, suggest_peak_separation(p.wf), 5, p.wf.fs);
    ASSERT_FALSE(pk.empty());
    EXPECT_TRUE(pk[0].above_threshold);
    EXPECT_LT(std::fabs(pk[0].range_m - p.ch.range_m), 0.5 * kSpeedOfLight / (2.0 * p.wf.fs));
}
