#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include "radar/types.hpp"

namespace radar {

// Lazy passband reconstruction  RF(t) = I(t) cos(2 pi fc t) - Q(t) sin(2 pi fc t)
// for the window t_k = t_start + k*dt, k = 0..n-1 (time measured from the first buffer sample).
// I(t), Q(t) are linearly interpolated between baseband samples (fine for views much shorter than 1/Fs)
// and are zero outside the buffer. The carrier phase is reduced modulo one cycle in double precision,
// so fc*t up to ~1e8 cycles stays accurate. Only the n visible points are ever allocated.
inline void reconstruct_passband(const CBuffer& x, double fs, double fc, double t_start, double dt,
                                 std::size_t n, RBuffer& out) {
    if (x.empty() || !(fs > 0.0) || !(fc >= 0.0) || !(dt > 0.0))
        throw std::invalid_argument("reconstruct_passband: bad arguments");
    out.assign(n, 0.0f);
    const double two_pi = 2.0 * std::numbers::pi;
    const long long N = static_cast<long long>(x.size());
    for (std::size_t k = 0; k < n; ++k) {
        const double t = t_start + static_cast<double>(k) * dt;
        const double u = t * fs;
        const double fl = std::floor(u);
        const double fr = u - fl;
        const long long i0 = static_cast<long long>(fl);
        double I = 0.0, Q = 0.0;
        if (i0 >= 0 && i0 < N) {
            I += (1.0 - fr) * static_cast<double>(x[static_cast<std::size_t>(i0)].real());
            Q += (1.0 - fr) * static_cast<double>(x[static_cast<std::size_t>(i0)].imag());
        }
        if (i0 + 1 >= 0 && i0 + 1 < N) {
            I += fr * static_cast<double>(x[static_cast<std::size_t>(i0 + 1)].real());
            Q += fr * static_cast<double>(x[static_cast<std::size_t>(i0 + 1)].imag());
        }
        const double cycles = fc * t;
        const double ph = two_pi * (cycles - std::floor(cycles));
        out[k] = static_cast<float>(I * std::cos(ph) - Q * std::sin(ph));
    }
}

}  // namespace radar
