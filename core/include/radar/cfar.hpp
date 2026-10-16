#pragma once
#include <cmath>
#include <string>
#include <vector>
#include "radar/channel.hpp"   // kSpeedOfLight
#include "radar/types.hpp"

namespace radar {

struct CfarParams {
    int guard_cells = 20;   // each side; use >= compressed mainlobe width (in samples)
    int train_cells = 64;   // each side
    double pfa = 1e-6;
};

struct CfarResult {
    RBuffer threshold;                     // per-cell threshold (same units as input power)
    std::vector<std::size_t> detections;   // one index (the local maximum) per contiguous run
    std::size_t cells_over_threshold = 0;  // raw count, before clustering
    double alpha = 0.0;
};

inline std::string validate_cfar(const CfarParams& p, std::size_t n) {
    if (p.guard_cells < 0) return "CFAR guard cells must be >= 0.";
    if (p.train_cells < 4) return "CFAR training cells must be >= 4 per side.";
    if (!(p.pfa > 0.0 && p.pfa < 1.0)) return "CFAR Pfa must be in (0, 1).";
    if (n < static_cast<std::size_t>(2 * (p.guard_cells + p.train_cells) + 1))
        return "CFAR window is larger than the data.";
    return "";
}

// Square-law detector, exponential (Gaussian-noise) power statistics.
inline double cfar_alpha(int n_train_total, double pfa) {
    const double n = static_cast<double>(n_train_total);
    return n * (std::pow(pfa, -1.0 / n) - 1.0);
}

inline double bin_to_range_m(std::size_t bin, double fs) {
    return kSpeedOfLight * static_cast<double>(bin) / (2.0 * fs);
}

// Cell-averaging CFAR on a circular power array (the range profile is periodic with the PRI).
inline CfarResult ca_cfar(const RBuffer& power, const CfarParams& p) {
    const long long n = static_cast<long long>(power.size());
    if (auto e = validate_cfar(p, power.size()); !e.empty()) throw std::invalid_argument(e);

    std::vector<double> pre(static_cast<std::size_t>(n) + 1, 0.0);
    for (long long i = 0; i < n; ++i)
        pre[static_cast<std::size_t>(i) + 1] = pre[static_cast<std::size_t>(i)] + static_cast<double>(power[static_cast<std::size_t>(i)]);

    auto window_sum = [&](long long start, long long len) {
        const long long s = ((start % n) + n) % n;
        const long long e = s + len;
        if (e <= n) return pre[static_cast<std::size_t>(e)] - pre[static_cast<std::size_t>(s)];
        return pre[static_cast<std::size_t>(n)] - pre[static_cast<std::size_t>(s)] + pre[static_cast<std::size_t>(e - n)];
    };

    const long long g = p.guard_cells, t = p.train_cells;
    CfarResult r;
    r.alpha = cfar_alpha(static_cast<int>(2 * t), p.pfa);
    r.threshold.resize(static_cast<std::size_t>(n));

    std::vector<char> hit(static_cast<std::size_t>(n), 0);
    for (long long i = 0; i < n; ++i) {
        const double lead = window_sum(i - g - t, t);
        const double lag = window_sum(i + g + 1, t);
        const double mean = (lead + lag) / static_cast<double>(2 * t);
        const double thr = r.alpha * mean;
        r.threshold[static_cast<std::size_t>(i)] = static_cast<float>(thr);
        if (static_cast<double>(power[static_cast<std::size_t>(i)]) > thr) {
            hit[static_cast<std::size_t>(i)] = 1;
            ++r.cells_over_threshold;
        }
    }

    // One detection per contiguous run: the strongest cell of the run
    long long i = 0;
    while (i < n) {
        if (!hit[static_cast<std::size_t>(i)]) { ++i; continue; }
        long long best = i;
        while (i < n && hit[static_cast<std::size_t>(i)]) {
            if (power[static_cast<std::size_t>(i)] > power[static_cast<std::size_t>(best)]) best = i;
            ++i;
        }
        r.detections.push_back(static_cast<std::size_t>(best));
    }
    return r;
}

}  // namespace radar
