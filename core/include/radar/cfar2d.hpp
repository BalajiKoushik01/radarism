#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include "radar/cfar.hpp"
#include "radar/types.hpp"

namespace radar {

struct Cfar2dParams {
    int guard_range = 10;     // cells each side
    int guard_doppler = 2;
    int train_range = 20;
    int train_doppler = 4;
    double pfa = 1e-6;
};

struct Cfar2dDetection {
    std::size_t row = 0, col = 0;
    float power = 0.0f, threshold = 0.0f;
};

struct Cfar2dResult {
    std::vector<Cfar2dDetection> detections;   // local maxima of over-threshold cells, strongest first
    std::size_t cells_over_threshold = 0;
};

inline std::string validate_cfar2d(const Cfar2dParams& p) {
    if (p.guard_range < 0 || p.guard_doppler < 0) return "CFAR guard cells must be >= 0.";
    if (p.train_range < 1 || p.train_doppler < 0) return "CFAR training cells: range >= 1, Doppler >= 0.";
    if (!(p.pfa > 0.0 && p.pfa < 1.0)) return "CFAR Pfa must be in (0, 1).";
    return "";
}

// Cell-averaging CFAR on a rows x cols power map (row-major). Windows are clamped at the map edges, and the
// per-cell scale alpha adapts to the actual number of training cells, so Pfa holds at the edges too.
inline Cfar2dResult ca_cfar_2d(const RBuffer& map, std::size_t rows, std::size_t cols, const Cfar2dParams& p) {
    if (auto e = validate_cfar2d(p); !e.empty()) throw std::invalid_argument(e);
    if (rows == 0 || cols == 0 || map.size() != rows * cols) throw std::invalid_argument("ca_cfar_2d: bad map size");

    const std::size_t W = cols + 1;
    std::vector<double> S((rows + 1) * W, 0.0);   // summed-area table
    for (std::size_t i = 0; i < rows; ++i) {
        double run = 0.0;
        for (std::size_t j = 0; j < cols; ++j) {
            run += static_cast<double>(map[i * cols + j]);
            S[(i + 1) * W + (j + 1)] = S[i * W + (j + 1)] + run;
        }
    }
    const long long R = static_cast<long long>(rows), C = static_cast<long long>(cols);
    auto at = [&](long long r, long long c) { return S[static_cast<std::size_t>(r) * W + static_cast<std::size_t>(c)]; };
    auto rect = [&](long long r0, long long r1, long long c0, long long c1, long long& count) {   // inclusive, clamped
        r0 = std::max(r0, 0LL); c0 = std::max(c0, 0LL);
        r1 = std::min(r1, R - 1); c1 = std::min(c1, C - 1);
        if (r0 > r1 || c0 > c1) { count = 0; return 0.0; }
        count = (r1 - r0 + 1) * (c1 - c0 + 1);
        return at(r1 + 1, c1 + 1) - at(r0, c1 + 1) - at(r1 + 1, c0) + at(r0, c0);
    };

    const long long gr = p.guard_range, gd = p.guard_doppler, tr = p.train_range, td = p.train_doppler;
    std::vector<float> thr(map.size());
    std::vector<char> over(map.size(), 0);
    Cfar2dResult res;
    long long last_nt = -1;
    double alpha = 0.0;
    for (long long i = 0; i < R; ++i) {
        for (long long j = 0; j < C; ++j) {
            long long n_out = 0, n_in = 0;
            const double s_out = rect(i - gd - td, i + gd + td, j - gr - tr, j + gr + tr, n_out);
            const double s_in = rect(i - gd, i + gd, j - gr, j + gr, n_in);
            const long long nt = n_out - n_in;
            float t = std::numeric_limits<float>::infinity();
            if (nt >= 4) {
                if (nt != last_nt) { alpha = cfar_alpha(static_cast<int>(nt), p.pfa); last_nt = nt; }
                t = static_cast<float>(alpha * std::max(s_out - s_in, 0.0) / static_cast<double>(nt));
            }
            const std::size_t idx = static_cast<std::size_t>(i) * cols + static_cast<std::size_t>(j);
            thr[idx] = t;
            if (map[idx] > t) { over[idx] = 1; ++res.cells_over_threshold; }
        }
    }
    for (long long i = 0; i < R; ++i) {
        for (long long j = 0; j < C; ++j) {
            const std::size_t idx = static_cast<std::size_t>(i) * cols + static_cast<std::size_t>(j);
            if (!over[idx]) continue;
            bool is_max = true;
            for (long long di = -1; di <= 1 && is_max; ++di)
                for (long long dj = -1; dj <= 1; ++dj) {
                    if (di == 0 && dj == 0) continue;
                    const long long ni = i + di, nj = j + dj;
                    if (ni < 0 || nj < 0 || ni >= R || nj >= C) continue;
                    if (map[static_cast<std::size_t>(ni) * cols + static_cast<std::size_t>(nj)] > map[idx]) { is_max = false; break; }
                }
            if (is_max) res.detections.push_back({static_cast<std::size_t>(i), static_cast<std::size_t>(j), map[idx], thr[idx]});
        }
    }
    std::sort(res.detections.begin(), res.detections.end(),
              [](const Cfar2dDetection& a, const Cfar2dDetection& b) { return a.power > b.power; });
    return res;
}

}  // namespace radar
