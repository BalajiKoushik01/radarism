#pragma once
#include <cmath>
#include <stdexcept>
#include <volk/volk.h>
#include "radar/fft.hpp"
#include "radar/types.hpp"

namespace radar {

// Reusable matched filter: FFT plans and the replica spectrum are built once, then compress() runs per pulse.
//   y[n] = (1/normalization) * sum_m rx[m+n] * conj(h[m])     (circular)
class RangeCompressor {
public:
    RangeCompressor(const CBuffer& replica, double normalization)
        : n_(replica.size()), fwd_(n_, true), inv_(n_, false), TX_(n_), RX_(n_), prod_(n_) {
        if (n_ < 2) throw std::invalid_argument("RangeCompressor: replica too short");
        if (!(normalization > 0.0)) throw std::invalid_argument("RangeCompressor: normalization must be > 0");
        fwd_.execute(replica.data(), TX_.data());
        scale_ = static_cast<float>(1.0 / (normalization * static_cast<double>(n_)));
    }
    std::size_t size() const { return n_; }

    void compress(const CBuffer& rx, CBuffer& y) {
        if (rx.size() != n_) throw std::invalid_argument("RangeCompressor: length mismatch");
        y.resize(n_);
        fwd_.execute(rx.data(), RX_.data());
        volk_32fc_x2_multiply_conjugate_32fc(prod_.data(), RX_.data(), TX_.data(), static_cast<unsigned int>(n_));
        inv_.execute(prod_.data(), y.data());
        for (auto& s : y) s *= scale_;
    }

private:
    std::size_t n_;
    Fft fwd_, inv_;
    CBuffer TX_, RX_, prod_;
    float scale_ = 1.0f;
};

}  // namespace radar
