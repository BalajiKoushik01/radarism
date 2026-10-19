#pragma once
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>
#include "radar/simulation.hpp"

namespace radar {

// ---------------- theory ----------------

// Marcum Q1(a, b) = P(noncentral chi-square, 2 dof, lambda = a^2 > b^2): a Poisson mixture of regularized
// upper incomplete gamma functions. For a square-law detector with a fixed threshold T = -ln(Pfa) on a
// non-fluctuating target: Pd = Q1(sqrt(2*SNR), sqrt(2*T)).
inline double marcum_q1(double a, double b) {
    if (!(b > 0.0)) return 1.0;
    if (!(a > 0.0)) return std::exp(-0.5 * b * b);
    const double lam = 0.5 * a * a, x = 0.5 * b * b;
    const double logx = std::log(x), logl = std::log(lam);
    const int kmax = static_cast<int>(lam + 12.0 * std::sqrt(lam) + 40.0);
    double sum = 0.0, q = 0.0;   // q = Q(k+1, x) accumulated term by term
    for (int k = 0; k <= kmax; ++k) {
        const double lg = std::lgamma(static_cast<double>(k) + 1.0);
        q += std::exp(-x + k * logx - lg);
        sum += std::exp(k * logl - lam - lg) * std::min(q, 1.0);
    }
    return std::clamp(sum, 0.0, 1.0);
}

inline double pd_theory(double snr_out_db, double pfa) {
    const double snr = std::pow(10.0, snr_out_db / 10.0);
    return marcum_q1(std::sqrt(2.0 * snr), std::sqrt(-2.0 * std::log(pfa)));
}

// Albersheim's approximation: single-sample SNR (dB) needed for (Pd, Pfa) with n non-coherently integrated
// pulses. Quoted accuracy about 0.2 dB for 1e-7 < Pfa < 1e-3 and 0.1 < Pd < 0.9.
inline double albersheim_snr_db(double pd, double pfa, int n = 1) {
    const double A = std::log(0.62 / pfa), B = std::log(pd / (1.0 - pd));
    return -5.0 * std::log10(static_cast<double>(n)) +
           (6.2 + 4.54 / std::sqrt(static_cast<double>(n) + 0.44)) * std::log10(A + 0.12 * A * B + 1.7 * B);
}

inline double theory_snr_for_pd(double pd, double pfa) {
    double lo = -30.0, hi = 60.0;
    for (int i = 0; i < 80; ++i) {
        const double mid = 0.5 * (lo + hi);
        (pd_theory(mid, pfa) < pd ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

inline void wilson_interval(int k, int n, double& lo, double& hi) {
    const double z = 1.96, p = static_cast<double>(k) / n;
    const double den = 1.0 + z * z / n;
    const double c = (p + z * z / (2.0 * n)) / den;
    const double h = z * std::sqrt(p * (1.0 - p) / n + z * z / (4.0 * static_cast<double>(n) * n)) / den;
    lo = c - h;
    hi = c + h;
}

// Matched-filter SNR gain (dB) of the (optionally tapered) replica, exactly as simulate() computes it.
inline double matched_filter_gain_db(const SimulationParams& p) {
    const Waveform tx = generate(p.wf);
    const CBuffer h = weighted_replica(tx.samples, tx.n_pulse, p.mf_weighting);
    double coh = 0.0, eh = 0.0;
    for (std::size_t i = 0; i < tx.n_pulse; ++i) {
        coh += static_cast<double>(std::real(std::conj(h[i]) * tx.samples[i]));
        eh += static_cast<double>(std::norm(h[i]));
    }
    return lin_to_db(coh * coh / eh);
}

// ---------------- Monte Carlo ----------------

template <class F>
void parallel_for(std::size_t n, unsigned threads, F fn) {
    if (n == 0) return;
    const std::size_t nt = std::max<std::size_t>(1, std::min<std::size_t>(threads, n));
    std::atomic<std::size_t> next{0};
    std::exception_ptr err;
    std::mutex m;
    auto worker = [&] {
        try {
            for (;;) {
                const std::size_t i = next.fetch_add(1);
                if (i >= n) break;
                fn(i);
            }
        } catch (...) {
            std::lock_guard<std::mutex> lk(m);
            if (!err) err = std::current_exception();
            next.store(n);
        }
    };
    std::vector<std::thread> pool;
    for (std::size_t t = 1; t < nt; ++t) pool.emplace_back(worker);
    worker();
    for (auto& th : pool) th.join();
    if (err) std::rethrow_exception(err);
}

struct McParams {
    SimulationParams sim;          // waveform, channel (SNR is overridden), MF weighting, 1D CFAR window
    double pfa = 1e-3;             // design Pfa for both the fixed threshold and the CA-CFAR
    std::size_t target_bin = 0;    // integer delay in samples; 0 = 30 % of the PRI
    int trials = 400;              // per SNR point
    int noise_trials = 200;        // noise-only runs (each contributes N cells to the Pfa estimate)
    unsigned threads = 0;          // 0 = hardware concurrency
    double match_tol_bins = 2.0;   // a CFAR detection within this many bins of the target counts as a hit
};

struct McPoint {
    double snr_out_db = 0.0, snr_in_db = 0.0;
    int trials = 0, hits_fixed = 0, hits_cfar = 0;
    double pd_theory = 0.0, pd_fixed = 0.0, pd_cfar = 0.0;
    double lo_fixed = 0.0, hi_fixed = 0.0, lo_cfar = 0.0, hi_cfar = 0.0;   // 95 % Wilson intervals
};

struct McPfa {
    int runs = 0;
    std::size_t cells = 0, over_fixed = 0, over_cfar = 0, clusters = 0;
    double pfa_design = 0.0, pfa_fixed = 0.0, pfa_cfar = 0.0, clusters_per_run = 0.0;
    double alpha = 0.0, expected_cfar_loss_db = 0.0, effective_cells = 0.0;
};

struct McResult {
    std::vector<McPoint> curve;
    McPfa pfa;
    double mf_gain_db = 0.0;
    std::size_t target_bin = 0, n = 0, n_pulse = 0;
    unsigned threads = 1;
    double theory_snr90 = 0.0, albersheim_snr90 = 0.0;
    double measured_snr90_fixed = std::numeric_limits<double>::quiet_NaN();
    double measured_snr90_cfar = std::numeric_limits<double>::quiet_NaN();
};

inline double mc_crossing_snr(const std::vector<McPoint>& c, bool cfar, double level = 0.9) {
    for (std::size_t i = 1; i < c.size(); ++i) {
        const double a = cfar ? c[i - 1].pd_cfar : c[i - 1].pd_fixed;
        const double b = cfar ? c[i].pd_cfar : c[i].pd_fixed;
        if (a < level && b >= level) {
            const double f = (level - a) / (b - a);
            return c[i - 1].snr_out_db + f * (c[i].snr_out_db - c[i - 1].snr_out_db);
        }
    }
    return std::numeric_limits<double>::quiet_NaN();
}

// snr_out_db: matched-filter OUTPUT SNR points (the axis the theory uses); input SNR = output - MF gain.
inline McResult run_monte_carlo(const McParams& mp, std::vector<double> snr_out_db,
                                const std::atomic<bool>* cancel = nullptr) {
    SimulationParams base = mp.sim;
    if (base.wf.type == WaveformType::CW) throw std::invalid_argument("Monte Carlo needs a pulsed waveform.");
    if (!(mp.pfa > 0.0 && mp.pfa < 0.5)) throw std::invalid_argument("Design Pfa must be in (0, 0.5).");
    if (mp.trials < 10 || mp.noise_trials < 1) throw std::invalid_argument("Need >= 10 trials per point and >= 1 noise run.");
    if (auto e = validate(base.wf); !e.empty()) throw std::invalid_argument(e);

    base.cfar.pfa = mp.pfa;
    base.run_detection = true;
    base.phase_noise_enabled = false;
    base.ch.velocity_mps = 0.0;
    base.ch.snr_mode = SnrMode::Manual;
    base.ch.snr_db = 0.0;

    const std::size_t N = total_samples(base.wf);
    const std::size_t np = static_cast<std::size_t>(std::llround(base.wf.tau * base.wf.fs));
    const std::size_t d = mp.target_bin ? mp.target_bin : static_cast<std::size_t>(0.3 * static_cast<double>(N));
    if (d < np || d + np >= N) throw std::invalid_argument("Target bin must satisfy pulse length <= bin < N - pulse length.");
    base.ch.range_m = kSpeedOfLight * (static_cast<double>(d) / base.wf.fs) / 2.0;   // integer delay: no straddle loss
    if (auto e = validate_all(base); !e.empty()) throw std::invalid_argument(e);

    std::sort(snr_out_db.begin(), snr_out_db.end());
    McResult res;
    res.mf_gain_db = matched_filter_gain_db(base);
    res.target_bin = d;
    res.n = N;
    res.n_pulse = np;
    res.threads = mp.threads ? mp.threads : std::max(1u, std::thread::hardware_concurrency());
    const double T = -std::log(mp.pfa);   // fixed threshold in units of the matched-filter noise floor

    const ResultPtr warm = simulate(base);   // warm-up on one thread (initializes VOLK dispatch) and CFAR readout
    const double warm_alpha = warm->cfar.alpha, warm_neff = warm->cfar_effective_cells;

    const std::size_t P = snr_out_db.size(), trials = static_cast<std::size_t>(mp.trials);
    std::vector<std::atomic<int>> hf(P), hc(P);
    for (auto& a : hf) a.store(0);
    for (auto& a : hc) a.store(0);
    parallel_for(P * trials, res.threads, [&](std::size_t k) {
        if (cancel && cancel->load(std::memory_order_relaxed)) throw Cancelled();
        const std::size_t i = k / trials, j = k % trials;
        SimulationParams p = base;
        p.ch.snr_db = snr_out_db[i] - res.mf_gain_db;
        p.ch.seed = base.ch.seed + (static_cast<std::uint64_t>(i) + 1) * 1000003ULL + j;
        const ResultPtr r = simulate(p);
        if (static_cast<double>(r->mf_power[d]) > T * r->mf_noise_floor_w) hf[i].fetch_add(1);
        for (double b : r->detection_bin)
            if (std::fabs(b - static_cast<double>(d)) <= mp.match_tol_bins) { hc[i].fetch_add(1); break; }
    });

    std::atomic<std::size_t> over_fixed{0}, over_cfar{0}, clusters{0};
    parallel_for(static_cast<std::size_t>(mp.noise_trials), res.threads, [&](std::size_t j) {
        if (cancel && cancel->load(std::memory_order_relaxed)) throw Cancelled();
        SimulationParams p = base;
        p.ch.snr_db = -200.0;   // echo negligible: noise only
        p.ch.seed = base.ch.seed + 7777777ULL + j;
        const ResultPtr r = simulate(p);
        const double thr = T * r->mf_noise_floor_w;
        std::size_t c = 0;
        for (float v : r->mf_power) if (static_cast<double>(v) > thr) ++c;
        over_fixed.fetch_add(c);
        over_cfar.fetch_add(r->cfar.cells_over_threshold);
        clusters.fetch_add(r->cfar.detections.size());
    });

    for (std::size_t i = 0; i < P; ++i) {
        McPoint pt;
        pt.snr_out_db = snr_out_db[i];
        pt.snr_in_db = snr_out_db[i] - res.mf_gain_db;
        pt.trials = mp.trials;
        pt.hits_fixed = hf[i].load();
        pt.hits_cfar = hc[i].load();
        pt.pd_theory = pd_theory(pt.snr_out_db, mp.pfa);
        pt.pd_fixed = static_cast<double>(pt.hits_fixed) / mp.trials;
        pt.pd_cfar = static_cast<double>(pt.hits_cfar) / mp.trials;
        wilson_interval(pt.hits_fixed, mp.trials, pt.lo_fixed, pt.hi_fixed);
        wilson_interval(pt.hits_cfar, mp.trials, pt.lo_cfar, pt.hi_cfar);
        res.curve.push_back(pt);
    }
    McPfa& f = res.pfa;
    f.runs = mp.noise_trials;
    f.cells = N * static_cast<std::size_t>(mp.noise_trials);
    f.over_fixed = over_fixed.load();
    f.over_cfar = over_cfar.load();
    f.clusters = clusters.load();
    f.pfa_design = mp.pfa;
    f.pfa_fixed = static_cast<double>(f.over_fixed) / static_cast<double>(f.cells);
    f.pfa_cfar = static_cast<double>(f.over_cfar) / static_cast<double>(f.cells);
    f.clusters_per_run = static_cast<double>(f.clusters) / mp.noise_trials;
    f.alpha = warm_alpha;
    f.effective_cells = warm_neff;
    f.expected_cfar_loss_db = 10.0 * std::log10(f.alpha / T);

    res.theory_snr90 = theory_snr_for_pd(0.9, mp.pfa);
    res.albersheim_snr90 = albersheim_snr_db(0.9, mp.pfa, 1);
    res.measured_snr90_fixed = mc_crossing_snr(res.curve, false);
    res.measured_snr90_cfar = mc_crossing_snr(res.curve, true);
    return res;
}

}  // namespace radar
