#include <gtest/gtest.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <random>
#include "radar/cfar2d.hpp"
#include "radar/cpi.hpp"

using namespace radar;

static SimulationParams base() {
    SimulationParams p;
    p.wf.type = WaveformType::LFM; p.wf.fs = 1e6; p.wf.pri = 1e-3; p.wf.tau = 100e-6; p.wf.bandwidth = 2e5;  // N=1000, pulse=100
    p.ch.fc = 10e9; p.ch.range_m = 30e3; p.ch.snr_db = 10.0;
    return p;
}
static double v_for_doppler(double fd, double fc) { return fd * kSpeedOfLight / (2.0 * fc); }
static std::size_t argmax(const RBuffer& p) {
    return static_cast<std::size_t>(std::max_element(p.begin(), p.end()) - p.begin());
}

TEST(Cpi, DopplerPeakAtExpectedRowRangeBinAndUnityPower) {
    auto p = base();
    p.ch.velocity_mps = v_for_doppler(250.0, p.ch.fc);   // +250 Hz = exactly 16 bins at PRF 1 kHz, L = 64
    CpiParams c; c.pulses = 64; c.gate_len = 256;
    auto r = simulate_cpi(p, c);
    ASSERT_EQ(r->doppler_rows, 64u);
    const std::size_t pk = argmax(r->rd_power);
    EXPECT_EQ(pk / r->gate_len, 32u + 16u);
    const double r_mid = p.ch.range_m - p.ch.velocity_mps * 31.5 * p.wf.pri;
    const double d_mid = 2.0 * r_mid / kSpeedOfLight * p.wf.fs;
    const double abs_bin = static_cast<double>((r->gate_start + pk % r->gate_len) % r->n);
    EXPECT_NEAR(abs_bin, d_mid, 1.5);
    EXPECT_NEAR(10.0 * std::log10(static_cast<double>(r->rd_power[pk]) / r->first_pulse.pr_w), 0.0, 1.5);
    EXPECT_NEAR(r->report.doppler_apparent_hz, 250.0, 1e-6);
    EXPECT_FALSE(r->report.doppler_ambiguous);
}

TEST(Cpi, StationaryTargetLandsOnZeroDopplerRow) {
    auto p = base();
    CpiParams c; c.pulses = 64; c.gate_len = 256;
    auto r = simulate_cpi(p, c);
    EXPECT_EQ(argmax(r->rd_power) / r->gate_len, 32u);
}

TEST(Cpi, FastTargetFoldsToApparentDopplerAndIsFlagged) {
    auto p = base();
    p.ch.velocity_mps = v_for_doppler(1250.0, p.ch.fc);   // PRF/2 = 500 Hz, so 1250 Hz folds to +250 Hz
    CpiParams c; c.pulses = 64; c.gate_len = 256;
    auto r = simulate_cpi(p, c);
    EXPECT_TRUE(r->report.doppler_ambiguous);
    EXPECT_NEAR(r->report.doppler_true_hz, 1250.0, 1e-6);
    EXPECT_NEAR(r->report.doppler_apparent_hz, 250.0, 1e-6);
    EXPECT_EQ(argmax(r->rd_power) / r->gate_len, 48u);
    EXPECT_NEAR(r->report.v_unambiguous_mps, 0.0299792458 * 1000.0 / 4.0, 1e-6);
}

TEST(Cpi, NoiseFloorMatchesReportedFloor) {
    auto p = base();
    p.ch.snr_db = -60.0;   // echo is negligible: the map is noise
    CpiParams c; c.pulses = 64; c.gate_len = 256;
    auto r = simulate_cpi(p, c);
    double sum = 0.0;
    for (float v : r->rd_power) sum += static_cast<double>(v);
    EXPECT_NEAR(sum / static_cast<double>(r->rd_power.size()) / r->rd_noise_floor_w, 1.0, 0.10);
}

TEST(Cpi, DetectorReportsTargetDopplerVelocityAndRange) {
    auto p = base();
    p.ch.velocity_mps = v_for_doppler(250.0, p.ch.fc);
    CpiParams c; c.pulses = 64; c.gate_len = 256;
    auto r = simulate_cpi(p, c);
    ASSERT_FALSE(r->detections.empty());
    const auto& d = r->detections.front();   // strongest first
    EXPECT_EQ(d.doppler_bin, 48u);
    EXPECT_NEAR(d.doppler_hz, 250.0, 1e-6);
    EXPECT_NEAR(d.velocity_mps, p.ch.velocity_mps, 1e-6);
    EXPECT_NEAR(d.range_m, p.ch.range_m, 200.0);
}

TEST(Cpi, ValidationAndCancellation) {
    EXPECT_TRUE(validate_cpi(base(), CpiParams()).empty());
    CpiParams c;
    auto p = base(); p.phase_noise_enabled = true;
    EXPECT_FALSE(validate_cpi(p, c).empty());
    p = base(); p.wf.type = WaveformType::CW;
    EXPECT_FALSE(validate_cpi(p, c).empty());
    p = base(); c.pulses = 1;
    EXPECT_FALSE(validate_cpi(p, c).empty());
    c = CpiParams();
    p = base(); p.ch.range_m = 50.0; p.ch.velocity_mps = 1000.0;   // range goes negative within the CPI
    EXPECT_FALSE(validate_cpi(p, c).empty());
    EXPECT_THROW(simulate_cpi(p, c), std::invalid_argument);
    std::atomic<bool> stop{true};
    EXPECT_THROW(simulate_cpi(base(), CpiParams(), &stop), Cancelled);
}

TEST(Cfar2d, FalseAlarmCountMatchesDesignOnIidExponentialMap) {
    const std::size_t rows = 128, cols = 512;
    std::mt19937_64 rng(11);
    std::exponential_distribution<float> ex(1.0f);
    RBuffer m(rows * cols);
    for (auto& v : m) v = ex(rng);
    Cfar2dParams cp; cp.guard_range = 2; cp.guard_doppler = 1; cp.train_range = 8; cp.train_doppler = 3; cp.pfa = 1e-3;
    auto r = ca_cfar_2d(m, rows, cols, cp);
    EXPECT_GT(r.cells_over_threshold, 30u);    // expect about 65
    EXPECT_LT(r.cells_over_threshold, 130u);
}

TEST(Cfar2d, FindsAnInjectedSpikeAndNothingElseOnAFlatFloor) {
    const std::size_t rows = 20, cols = 40;
    RBuffer m(rows * cols, 1.0f);
    m[10 * cols + 20] = 1000.0f;
    Cfar2dParams cp; cp.guard_range = 2; cp.guard_doppler = 1; cp.train_range = 6; cp.train_doppler = 2; cp.pfa = 1e-3;
    auto r = ca_cfar_2d(m, rows, cols, cp);
    ASSERT_EQ(r.detections.size(), 1u);
    EXPECT_EQ(r.detections[0].row, 10u);
    EXPECT_EQ(r.detections[0].col, 20u);
    EXPECT_THROW(ca_cfar_2d(m, rows, cols + 1, cp), std::invalid_argument);
}
