#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include "radar/channel.hpp"
#include "radar/matched_filter.hpp"
#include "radar/peak_interp.hpp"
#include "radar/simulation.hpp"
#include "radar/waveforms.hpp"

using namespace radar;

static WaveformParams lfm() {
    WaveformParams p; p.type = WaveformType::LFM;
    p.fs = 1e6; p.pri = 10e-3; p.tau = 2e-3; p.bandwidth = 2e5;   // N = 10000, pulse = 2000
    return p;
}
static double range_for_delay(double d, double fs) { return kSpeedOfLight * (d / fs) / 2.0; }
static std::size_t argmax(const RBuffer& p) {
    return static_cast<std::size_t>(std::max_element(p.begin(), p.end()) - p.begin());
}

TEST(PeakOffset, RecoversFractionalDelayOnNoiseFreeLfm) {
    auto wp = lfm();
    for (double d : {3000.0, 3000.3, 3000.4, 2999.7, 3000.49}) {
        ChannelParams c; c.snr_db = 20.0; c.range_m = range_for_delay(d, wp.fs);
        auto tx = generate(wp);
        auto rep = analyze(wp, c);
        auto echo = make_echo(tx, wp, c, rep);
        auto p = magnitude_squared(matched_filter(echo, tx.samples));
        const std::size_t pk = argmax(p);
        EXPECT_NEAR(static_cast<double>(pk) + peak_offset(p, pk), d, 0.1) << "delay " << d;
    }
}

TEST(PeakOffset, ExactOnGaussianPeakAndWrapsAroundTheEnd) {
    const std::size_t n = 100;
    RBuffer p(n);
    for (std::size_t i = 0; i < n; ++i) {
        double d = std::fabs(static_cast<double>(i) - 0.3);
        d = std::min(d, static_cast<double>(n) - d);
        p[i] = static_cast<float>(std::exp(-d * d / (2.0 * 1.5 * 1.5)));
    }
    EXPECT_NEAR(peak_offset(p, 0), 0.3, 1e-4);   // neighbours of bin 0 are bins 99 and 1
    RBuffer flat(n, 1.0f);
    EXPECT_EQ(peak_offset(flat, 10), 0.0);
    RBuffer valley(n, 1.0f);
    valley[10] = 0.1f;
    EXPECT_EQ(peak_offset(valley, 10), 0.0);
}

TEST(Weighting, LossMatchesClosedForm) {
    EXPECT_EQ(weighting_loss_db(MfWeighting::None, 2000), 0.0);
    EXPECT_NEAR(weighting_loss_db(MfWeighting::Hamming, 2000), -1.344, 0.05);
    EXPECT_NEAR(weighting_loss_db(MfWeighting::Hann, 2000), -1.761, 0.05);
}

TEST(Weighting, HammingKeepsUnityPeakAndLowersSidelobes) {
    auto wp = lfm();
    ChannelParams c; c.snr_db = 20.0; c.range_m = range_for_delay(3000.0, wp.fs);
    auto tx = generate(wp);
    auto rep = analyze(wp, c);
    auto echo = make_echo(tx, wp, c, rep);
    const long long N = static_cast<long long>(tx.samples.size());
    auto run = [&](MfWeighting w, double& peak_amp) {
        const CBuffer h = weighted_replica(tx.samples, tx.n_pulse, w);
        double coh = 0.0;
        for (std::size_t i = 0; i < tx.n_pulse; ++i) coh += std::real(std::conj(h[i]) * tx.samples[i]);
        auto y = matched_filter(echo, h, coh);
        auto p = magnitude_squared(y);
        peak_amp = std::abs(y[3000]);
        double worst = 0.0;
        for (long long i = 0; i < N; ++i) {
            long long dist = std::llabs(i - 3000);
            dist = std::min(dist, N - dist);
            if (dist > 20) worst = std::max(worst, static_cast<double>(p[static_cast<std::size_t>(i)]));
        }
        return 10.0 * std::log10(worst / static_cast<double>(p[3000]));
    };
    double a0 = 0, a1 = 0;
    const double s0 = run(MfWeighting::None, a0);
    const double s1 = run(MfWeighting::Hamming, a1);
    EXPECT_NEAR(a0, std::sqrt(rep.pr_w), 1e-3 * std::sqrt(rep.pr_w));
    EXPECT_NEAR(a1, std::sqrt(rep.pr_w), 1e-3 * std::sqrt(rep.pr_w));
    EXPECT_LT(s1, s0 - 15.0) << "unweighted " << s0 << " dB, Hamming " << s1 << " dB";
}

TEST(Pipeline, HammingReportsLossAndRefinedRangeIsNearTruth) {
    SimulationParams p;
    p.wf = lfm();
    p.ch.snr_db = -14.0;
    p.ch.range_m = range_for_delay(3000.4, p.wf.fs);
    p.mf_weighting = MfWeighting::Hamming;
    auto r = simulate(p);
    EXPECT_NEAR(r->report.mf_gain_db, 10.0 * std::log10(2000.0) + weighting_loss_db(MfWeighting::Hamming, 2000), 1e-3);
    EXPECT_GT(r->mf_noise_floor_w, r->report.pn_w / 2000.0);   // tapered replica: higher floor per unit peak
    ASSERT_EQ(r->detection_bin.size(), r->cfar.detections.size());
    bool found = false;
    for (double b : r->detection_bin) if (std::fabs(b - 3000.4) < 1.0) found = true;
    EXPECT_TRUE(found);
}
