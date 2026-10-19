#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
#include "radar/cfar.hpp"
#include "radar/cfar_neff.hpp"
#include "radar/montecarlo.hpp"

using namespace radar;

TEST(CfarNeff, ImpulseReplicaMeansIndependentCells) {
    CBuffer h(8, cf32{0.0f, 0.0f});
    h[0] = {1.0f, 0.0f};
    const auto c = noise_correlation_power(h, 1, 2 * (3 + 5));
    ASSERT_EQ(c.size(), 17u);
    EXPECT_NEAR(c[0], 1.0, 1e-6);
    for (std::size_t k = 1; k < c.size(); ++k) EXPECT_NEAR(c[k], 0.0, 1e-6);
    EXPECT_NEAR(cfar_effective_cells(c, 3, 5), 10.0, 1e-6);
}

TEST(CfarNeff, RectPulseMatchesClosedForm) {
    const std::size_t n = 50;
    CBuffer h(n, cf32{1.0f, 0.0f});
    const int g = 2, t = 4;
    const std::size_t K = 2 * (g + t);
    const auto c = noise_correlation_power(h, n, K);
    for (std::size_t k = 0; k <= K; ++k) {
        const double rho = 1.0 - static_cast<double>(k) / static_cast<double>(n);
        EXPECT_NEAR(c[k], rho * rho, 1e-4) << k;
    }
    std::vector<int> pos;
    for (int i = 1; i <= t; ++i) { pos.push_back(-(g + i)); pos.push_back(g + i); }
    double s = 0.0;
    for (int a : pos)
        for (int b : pos) {
            const double rho = 1.0 - std::abs(a - b) / static_cast<double>(n);
            s += rho * rho;
        }
    const double expected = static_cast<double>(pos.size() * pos.size()) / s;
    EXPECT_NEAR(cfar_effective_cells(c, g, t), expected, 1e-2);
    EXPECT_LT(expected, 3.0);   // strongly correlated: far fewer than the 8 training cells
}

TEST(CfarNeff, AlphaGrowsAsEffectiveCellsShrink) {
    EXPECT_GT(cfar_alpha(30.0, 1e-6), cfar_alpha(128.0, 1e-6));
    EXPECT_NEAR(cfar_alpha(2.0, 0.01), 18.0, 1e-9);
    EXPECT_NEAR(cfar_alpha(128.5, 1e-3), 128.5 * (std::pow(1e-3, -1.0 / 128.5) - 1.0), 1e-9);
}

static McParams params(double pfa, bool corr, int noise_runs) {
    McParams mp;
    mp.sim.wf.type = WaveformType::LFM;
    mp.sim.wf.fs = 1e6; mp.sim.wf.pri = 10e-3; mp.sim.wf.tau = 500e-6; mp.sim.wf.bandwidth = 2e5;
    mp.sim.cfar_correlation_aware = corr;
    mp.pfa = pfa; mp.trials = 10; mp.noise_trials = noise_runs; mp.threads = 4;
    return mp;
}

TEST(CfarNeff, CorrelationAwareThresholdRestoresDesignFalseAlarmRate) {
    for (double pfa : {1e-3, 1e-4}) {
        const auto off = run_monte_carlo(params(pfa, false, 120), {10.0});
        const auto on = run_monte_carlo(params(pfa, true, 120), {10.0});
        const double ro = off.pfa.pfa_cfar / pfa, rn = on.pfa.pfa_cfar / pfa;
        EXPECT_GT(ro, 1.3) << "pfa " << pfa;     // the iid assumption is optimistic on oversampled noise
        EXPECT_GT(rn, 0.5) << "pfa " << pfa;
        EXPECT_LT(rn, 1.7) << "pfa " << pfa;
        EXPECT_LT(rn, ro) << "pfa " << pfa;
        EXPECT_LT(on.pfa.effective_cells, 2.0 * 64);
    }
}
