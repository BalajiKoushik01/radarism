#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>
#include "radar/peak_interp.hpp"
#include "radar/types.hpp"
#include "radar/waveforms.hpp"

namespace radar {

struct Peak {
    std::size_t bin = 0;          // integer bin of the local maximum
    double bin_refined = 0.0;     // sub-bin position (log-parabola through the three bins around the maximum)
    double range_m = 0.0;         // range of the refined position
    double power = 0.0;           // linear power at the bin
    double snr_db = 0.0;          // dB above the matched-filter noise floor
    bool above_threshold = false;
};

// Width of the compressed pulse in bins: peaks closer than this belong to the same return.
inline std::size_t suggest_peak_separation(const WaveformParams& w) {
    double s = 3.0;
    switch (w.type) {
        case WaveformType::LFM: if (w.bandwidth > 0.0) s = 2.0 * w.fs / w.bandwidth; break;
        case WaveformType::Pulse: s = w.fs * w.tau; break;
        case WaveformType::Barker:
        case WaveformType::Frank: s = 2.0 * w.fs * w.tau / static_cast<double>(std::max(1, w.code_length)); break;
        default: break;
    }
    return std::max<std::size_t>(3, static_cast<std::size_t>(std::llround(s)));
}

// Peak detection on a matched-filter POWER profile (circular):
//  1. candidates = local maxima at least (threshold_db - 6 dB) above the noise floor (plus the global maximum),
//  2. ranked by power, accepted greedily with non-maximum suppression (min_sep bins),
//  3. each accepted peak refined to sub-bin accuracy,
//  4. at most max_peaks returned, strongest first; above_threshold marks peaks >= threshold_db above the floor.
// The strongest peak is always returned, even when it is below the threshold.
inline std::vector<Peak> find_peaks(const RBuffer& p, double floor_w, double threshold_db, std::size_t min_sep,
                                    std::size_t max_peaks, double fs) {
    std::vector<Peak> out;
    const std::size_t n = p.size();
    if (n < 3 || !(floor_w > 0.0) || max_peaks == 0 || !(fs > 0.0)) return out;
    float gmax = 0.0f;
    for (float v : p) gmax = std::max(gmax, v);
    const double cut = floor_w * std::pow(10.0, (threshold_db - 6.0) / 10.0);
    std::vector<std::size_t> cand;
    for (std::size_t i = 0; i < n; ++i) {
        const float a = p[(i + n - 1) % n], b = p[i], c = p[(i + 1) % n];
        if (b > a && b >= c && (static_cast<double>(b) >= cut || b == gmax)) cand.push_back(i);
    }
    std::sort(cand.begin(), cand.end(), [&](std::size_t x, std::size_t y) { return p[x] > p[y]; });
    const std::size_t sep = std::max<std::size_t>(min_sep, 1);
    for (std::size_t idx : cand) {
        bool ok = true;
        for (const Peak& q : out) {
            const std::size_t d = idx > q.bin ? idx - q.bin : q.bin - idx;
            if (std::min(d, n - d) < sep) { ok = false; break; }
        }
        if (!ok) continue;
        Peak pk;
        pk.bin = idx;
        pk.power = static_cast<double>(p[idx]);
        pk.snr_db = 10.0 * std::log10(std::max(pk.power / floor_w, 1e-30));
        pk.above_threshold = (pk.snr_db >= threshold_db);
        if (!pk.above_threshold && !out.empty()) break;   // sorted by power: everything after is weaker
        const double b = std::fmod(static_cast<double>(idx) + peak_offset(p, idx) + static_cast<double>(n), static_cast<double>(n));
        pk.bin_refined = b;
        pk.range_m = range_from_bin(b, fs);
        out.push_back(pk);
        if (out.size() >= max_peaks) break;
    }
    return out;
}

}  // namespace radar
