#include <gtest/gtest.h>
#include <cmath>
#include "radar/cpi.hpp"
#include "radar/simulation.hpp"

using namespace radar;

static double range_for_bin(double b, double fs) { return kSpeedOfLight * (b / fs) / 2.0; }

static SimulationParams single() {
    SimulationParams p;
    p.wf.type = WaveformType::LFM; p.wf.fs = 1e6; p.wf.pri = 10e-3; p.wf.tau = 500e-6; p.wf.bandwidth = 2e5;
    p.ch.range_m = range_for_bin(3000.0, p.wf.fs);
    return p;
}

TEST(Targets, TwoTargetsAreBothDetected) {
    auto p = single();
    p.ch.snr_db = -5.0;
    p.extra_targets.push_back({range_for_bin(5000.0, p.wf.fs), 0.0, -3.0});
    p.cfar.pfa = 1e-6; p.cfar_correlation_aware = true;
    const auto r = simulate(p);
    bool a = false, b = false;
    for (double bin : r->detection_bin) {
        if (std::fabs(bin - 3000.0) < 2.0) a = true;
        if (std::fabs(bin - 5000.0) < 2.0) b = true;
    }
    EXPECT_TRUE(a);
    EXPECT_TRUE(b);
}

TEST(Targets, RelativeLevelIsRespected) {
    auto p = single();
    p.ch.snr_db = 13.0;
    p.extra_targets.push_back({range_for_bin(5000.0, p.wf.fs), 0.0, -3.0});
    const auto r = simulate(p);
    const double db = 10.0 * std::log10(static_cast<double>(r->mf_power[5000]) / static_cast<double>(r->mf_power[3000]));
    EXPECT_NEAR(db, -3.0, 0.3);
}

TEST(Targets, ValidationMessagesNameTheTarget) {
    auto p = single();
    p.extra_targets.push_back({-5.0, 0.0, 0.0});
    EXPECT_NE(validate_all(p).find("Target 2"), std::string::npos);
    EXPECT_THROW(simulate(p), std::invalid_argument);
    p.extra_targets.assign(9, TargetSpec{40e3, 0.0, 0.0});
    EXPECT_FALSE(validate_all(p).empty());
}

TEST(Targets, CpiShowsBothTargetsAtTheirDopplersAndRanges) {
    SimulationParams p;
    p.wf.type = WaveformType::LFM; p.wf.fs = 1e6; p.wf.pri = 1e-3; p.wf.tau = 100e-6; p.wf.bandwidth = 2e5;
    p.ch.fc = 10e9; p.ch.range_m = 30e3; p.ch.snr_db = -10.0;
    p.ch.velocity_mps = 250.0 * kSpeedOfLight / (2.0 * p.ch.fc);
    p.extra_targets.push_back({33e3, -125.0 * kSpeedOfLight / (2.0 * p.ch.fc), -3.0});
    CpiParams c;
    c.pulses = 64; c.gate_len = 256;
    const auto r = simulate_cpi(p, c);
    bool a = false, b = false;
    for (const auto& d : r->detections) {
        if (std::fabs(d.doppler_hz - 250.0) < 20.0 && std::fabs(d.range_m - 30e3) < 400.0) a = true;
        if (std::fabs(d.doppler_hz + 125.0) < 20.0 && std::fabs(d.range_m - 33e3) < 400.0) b = true;
    }
    EXPECT_TRUE(a);
    EXPECT_TRUE(b);
}
