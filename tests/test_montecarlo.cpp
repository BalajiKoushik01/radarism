#include <gtest/gtest.h>
#include <atomic>
#include <cmath>
#include "radar/montecarlo.hpp"

using namespace radar;

static McParams base_params() {
    McParams mp;
    mp.sim.wf.type = WaveformType::LFM;
    mp.sim.wf.fs = 1e6; mp.sim.wf.pri = 10e-3; mp.sim.wf.tau = 500e-6; mp.sim.wf.bandwidth = 2e5;   // N=10000, pulse=500
    mp.pfa = 1e-3; mp.trials = 300; mp.noise_trials = 100; mp.threads = 4;
    return mp;
}

TEST(Theory, MarcumQ1KnownValuesAndLimits) {
    EXPECT_NEAR(marcum_q1(1.0, 1.0), 0.7329, 2e-4);
    EXPECT_NEAR(marcum_q1(0.0, 3.0), std::exp(-4.5), 1e-12);
    EXPECT_EQ(marcum_q1(2.0, 0.0), 1.0);
    EXPECT_NEAR(marcum_q1(12.0, 3.0), 1.0, 1e-6);
    EXPECT_GT(marcum_q1(2.0, 2.0), marcum_q1(1.5, 2.0));   // increasing in a
    EXPECT_LT(marcum_q1(2.0, 3.0), marcum_q1(2.0, 2.0));   // decreasing in b
}

TEST(Theory, PdAtZeroSnrIsPfaAndAlbersheimAgreesWithMarcum) {
    EXPECT_NEAR(pd_theory(-80.0, 1e-3), 1e-3, 1e-5);
    EXPECT_NEAR(albersheim_snr_db(0.9, 1e-6, 1), 13.114, 0.05);
    for (double pfa : {1e-4, 1e-6})
        EXPECT_NEAR(pd_theory(albersheim_snr_db(0.9, pfa), pfa), 0.9, 0.025) << "pfa " << pfa;
    const double s = theory_snr_for_pd(0.9, 1e-6);
    EXPECT_NEAR(pd_theory(s, 1e-6), 0.9, 1e-6);
    EXPECT_NEAR(s, albersheim_snr_db(0.9, 1e-6), 0.3);
}

TEST(MonteCarlo, NoiseOnlyFalseAlarmRateMatchesDesign) {
    auto mp = base_params();
    mp.pfa = 1e-2; mp.noise_trials = 60;
    const auto r = run_monte_carlo(mp, {});
    EXPECT_NEAR(r.pfa.pfa_fixed / 1e-2, 1.0, 0.12);   // validates the noise scaling and the MF floor
    EXPECT_GT(r.pfa.pfa_cfar / 1e-2, 0.7);            // CA-CFAR on correlated noise: near design, typically a little above
    EXPECT_LT(r.pfa.pfa_cfar / 1e-2, 2.5);
}

static void check_curve(const McParams& mp, const std::vector<double>& snrs) {
    const auto r = run_monte_carlo(mp, snrs);
    ASSERT_EQ(r.curve.size(), snrs.size());
    for (const auto& pt : r.curve) {
        const double tol = 4.0 * std::sqrt(pt.pd_theory * (1.0 - pt.pd_theory) / pt.trials) + 0.01;
        EXPECT_NEAR(pt.pd_fixed, pt.pd_theory, tol) << "snr_out " << pt.snr_out_db;
        EXPECT_NEAR(pt.pd_cfar, pt.pd_theory, 0.07) << "snr_out " << pt.snr_out_db;
    }
}

TEST(MonteCarlo, PdCurveMatchesMarcumUnweighted) {
    check_curve(base_params(), {8.0, 10.7, 13.0});
}

TEST(MonteCarlo, PdCurveMatchesMarcumWithHammingWeighting) {   // checks the weighting-loss bookkeeping end to end
    auto mp = base_params();
    mp.sim.mf_weighting = MfWeighting::Hamming;
    check_curve(mp, {9.0, 11.0});
}

TEST(MonteCarlo, ValidationAndCancellation) {
    auto mp = base_params();
    mp.sim.wf.type = WaveformType::CW;
    EXPECT_THROW(run_monte_carlo(mp, {10.0}), std::invalid_argument);
    mp = base_params(); mp.pfa = 0.7;
    EXPECT_THROW(run_monte_carlo(mp, {10.0}), std::invalid_argument);
    mp = base_params(); mp.target_bin = 100;   // inside the 500-sample pulse
    EXPECT_THROW(run_monte_carlo(mp, {10.0}), std::invalid_argument);
    mp = base_params();
    std::atomic<bool> stop{true};
    EXPECT_THROW(run_monte_carlo(mp, {10.0}, &stop), Cancelled);
}
