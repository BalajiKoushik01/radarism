#pragma once
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <volk/volk.h>
#include "radar/fft.hpp"
#include "radar/types.hpp"

namespace radar {

// Taper applied to the receive replica over the active pulse. Lowers range sidelobes at the cost of a small
// matched-filter SNR loss (a mismatched filter).
enum class MfWeighting { None, Hamming, Hann };

inline double taper_value(MfWeighting w, std::size_t i, std::size_t n) {
    if (w == MfWeighting::None || n < 2) return 1.0;
    const double c = std::cos(2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(n - 1));
    return (w == MfWeighting::Hamming) ? 0.54 - 0.46 * c : 0.5 - 0.5 * c;
}

inline CBuffer weighted_replica(const CBuffer& tx, std::size_t n_pulse, MfWeighting w) {
    CBuffer h = tx;
    if (w == MfWeighting::None) return h;
    for (std::size_t i = 0; i < n_pulse && i < h.size(); ++i)
        h[i] *= static_cast<float>(taper_value(w, i, n_pulse));
    return h;
}

// SNR change of the weighted filter versus the matched one for a constant-modulus pulse (dB, <= 0).
inline double weighting_loss_db(MfWeighting w, std::size_t n_pulse) {
    if (w == MfWeighting::None || n_pulse < 2) return 0.0;
    double s1 = 0.0, s2 = 0.0;
    for (std::size_t i = 0; i < n_pulse; ++i) {
        const double v = taper_value(w, i, n_pulse);
        s1 += v;
        s2 += v * v;
    }
    return 10.0 * std::log10(s1 * s1 / (static_cast<double>(n_pulse) * s2));
}

// Circular cross-correlation via FFT:  y[n] = (1/E) * sum_m rx[m+n] * conj(h[m]).
// Default E = sum |h|^2 (classic matched filter: |y[d]| = A for an echo of amplitude A).
// Pass normalization = sum Re(conj(h) * tx) to keep unity gain on the signal for a tapered replica.
inline CBuffer matched_filter(const CBuffer& rx, const CBuffer& replica, double normalization = 0.0) {
    const std::size_t n = rx.size();
    if (n < 2 || replica.size() != n)
        throw std::invalid_argument("matched_filter: rx and replica must have the same length (>= 2)");
    double e = 0.0;
    for (const auto& s : replica) e += static_cast<double>(std::norm(s));
    if (!(e > 0.0)) throw std::invalid_argument("matched_filter: replica has zero energy");
    if (normalization > 0.0) e = normalization;

    CBuffer RX(n), TX(n), prod(n), y(n);
    Fft fwd(n, true), inv(n, false);
    fwd.execute(rx.data(), RX.data());
    fwd.execute(replica.data(), TX.data());
    volk_32fc_x2_multiply_conjugate_32fc(prod.data(), RX.data(), TX.data(), static_cast<unsigned int>(n));
    inv.execute(prod.data(), y.data());   // unnormalized inverse: gives n * correlation

    const float scale = static_cast<float>(1.0 / (e * static_cast<double>(n)));
    for (auto& s : y) s *= scale;
    return y;
}

inline RBuffer magnitude_squared(const CBuffer& x) {
    RBuffer p(x.size());
    volk_32fc_magnitude_squared_32f(p.data(), x.data(), static_cast<unsigned int>(x.size()));
    return p;
}

}  // namespace radar
