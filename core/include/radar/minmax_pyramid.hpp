#pragma once
#include <algorithm>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>
#include "radar/types.hpp"

namespace radar {

// Min/max pyramid over a real-valued trace. Level j stores min/max per block of 16^(j+1) samples.
// query() costs O(16 * levels) per output column regardless of window size, so pan/zoom stays
// cheap at 10M+ samples. Every spike survives decimation (that is the point of min/max).
//
// Two modes:
//   - owning:   MinMaxPyramid(RBuffer)               takes the trace
//   - borrowed: MinMaxPyramid(ptr, n, stride)        views external data (e.g. the I or Q part of an
//               interleaved complex buffer); the data must outlive the pyramid.
class MinMaxPyramid {
public:
    static constexpr std::size_t kFactor = 16;

    MinMaxPyramid() = default;
    MinMaxPyramid(const MinMaxPyramid&) = delete;
    MinMaxPyramid& operator=(const MinMaxPyramid&) = delete;
    MinMaxPyramid(MinMaxPyramid&&) noexcept = default;
    MinMaxPyramid& operator=(MinMaxPyramid&&) noexcept = default;

    explicit MinMaxPyramid(RBuffer base) : base_(std::move(base)), p_(base_.data()), n_(base_.size()) { build(); }
    MinMaxPyramid(const float* data, std::size_t n, std::size_t stride = 1) : p_(data), n_(n), stride_(stride) { build(); }

    std::size_t size() const { return n_; }

    // Window [i0, i1) split into min(columns, i1-i0) columns; column c covers
    // [i0 + span*c/cols, i0 + span*(c+1)/cols). Returns the number of columns written.
    std::size_t query(std::size_t i0, std::size_t i1, std::size_t columns,
                      std::vector<float>& mn, std::vector<float>& mx) const {
        i1 = std::min(i1, n_);
        if (i0 >= i1 || columns == 0) { mn.clear(); mx.clear(); return 0; }
        const std::size_t span = i1 - i0;
        const std::size_t cols = std::min(columns, span);
        mn.resize(cols);
        mx.resize(cols);
        for (std::size_t c = 0; c < cols; ++c) {
            const std::size_t a = i0 + span * c / cols;
            const std::size_t b = i0 + span * (c + 1) / cols;   // b > a because cols <= span
            float lo = std::numeric_limits<float>::infinity();
            float hi = -std::numeric_limits<float>::infinity();
            int top = -1;
            for (std::size_t j = 0; j < levels_.size() && levels_[j].block <= b - a; ++j) top = static_cast<int>(j);
            scan(a, b, top, lo, hi);
            mn[c] = lo;
            mx[c] = hi;
        }
        return cols;
    }

private:
    struct Level {
        std::size_t block = 0;
        std::vector<float> mn, mx;
    };

    float at(std::size_t i) const { return p_[i * stride_]; }

    void build() {
        const std::size_t n = n_;
        std::size_t block = kFactor;
        while (n >= block) {
            Level lv;
            lv.block = block;
            const std::size_t count = (n + block - 1) / block;
            lv.mn.resize(count);
            lv.mx.resize(count);
            if (levels_.empty()) {
                for (std::size_t b = 0; b < count; ++b) {
                    const std::size_t s = b * block, e = std::min(n, s + block);
                    float lo = at(s), hi = lo;
                    for (std::size_t i = s + 1; i < e; ++i) { const float v = at(i); lo = std::min(lo, v); hi = std::max(hi, v); }
                    lv.mn[b] = lo; lv.mx[b] = hi;
                }
            } else {
                const Level& pl = levels_.back();
                const std::size_t pcount = pl.mn.size();
                for (std::size_t b = 0; b < count; ++b) {
                    const std::size_t s = b * kFactor, e = std::min(pcount, s + kFactor);
                    float lo = pl.mn[s], hi = pl.mx[s];
                    for (std::size_t i = s + 1; i < e; ++i) { lo = std::min(lo, pl.mn[i]); hi = std::max(hi, pl.mx[i]); }
                    lv.mn[b] = lo; lv.mx[b] = hi;
                }
            }
            levels_.push_back(std::move(lv));
            block *= kFactor;
        }
    }

    void scan(std::size_t a, std::size_t b, int level, float& lo, float& hi) const {
        if (a >= b) return;
        if (level < 0) {
            for (std::size_t i = a; i < b; ++i) { const float v = at(i); lo = std::min(lo, v); hi = std::max(hi, v); }
            return;
        }
        const Level& L = levels_[static_cast<std::size_t>(level)];
        const std::size_t first = (a + L.block - 1) / L.block;
        const std::size_t last = b / L.block;
        if (first >= last) { scan(a, b, level - 1, lo, hi); return; }
        scan(a, first * L.block, level - 1, lo, hi);
        for (std::size_t k = first; k < last; ++k) { lo = std::min(lo, L.mn[k]); hi = std::max(hi, L.mx[k]); }
        scan(last * L.block, b, level - 1, lo, hi);
    }

    RBuffer base_;                 // owning mode only
    const float* p_ = nullptr;     // points into base_ or into external data
    std::size_t n_ = 0;
    std::size_t stride_ = 1;
    std::vector<Level> levels_;
};

}  // namespace radar
