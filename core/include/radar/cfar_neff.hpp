#pragma once
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include "radar/fft.hpp"
#include "radar/types.hpp"

namespace radar {

// Matched-filter output noise is correlated between neighbouring range bins. For white complex Gaussian input,
// y[n] has complex correlation rho(k) = r(k)/r(0), r(k) = sum_m conj(h[m]) h[m+k] (h = receive replica), and the
// POWER samples |y|^2 are correlated with coefficient |rho(k)|^2. Returns c[k] = |rho(k)|^2 for k = 0..max_lag.
inline std::vector<double> noise_correlation_power(const CBuffer& h, std::size_t n_active, std::size_t max_lag) {
    if (n_active == 0 || n_active > h.size()) throw std::invalid_argument("noise_correlation_power: bad length");
    std::size_t M = 1;
    while (M < n_active + max_lag + 1) M <<= 1;   // zero padding: lags <= max_lag do not alias
    CBuffer a(M, cf32{0.0f, 0.0f}), A(M), P(M), r(M);
    std::copy(h.begin(), h.begin() + static_cast<std::ptrdiff_t>(n_active), a.begin());
    Fft fwd(M, true), inv(M, false);
    fwd.execute(a.data(), A.data());
    for (std::size_t i = 0; i < M; ++i) P[i] = cf32{std::norm(A[i]), 0.0f};
    inv.execute(P.data(), r.data());
    std::vector<double> c(max_lag + 1, 0.0);
    const double r0 = static_cast<double>(std::abs(r[0]));
    if (!(r0 > 0.0)) return c;
    for (std::size_t k = 0; k <= max_lag; ++k) {
        const double rho = static_cast<double>(std::abs(r[k])) / r0;
        c[k] = rho * rho;
    }
    return c;
}

// Independent-cell equivalent of the two CA-CFAR training blocks (positions -g-t..-g-1 and g+1..g+t around the
// cell under test): N_eff = N^2 / sum_ij c(|i-j|). The training mean is then ~Gamma(N_eff), so
// Pfa = (1 + alpha/N_eff)^(-N_eff). Equals N for independent cells. Needs c up to lag 2*(g+t).
inline double cfar_effective_cells(const std::vector<double>& c, int guard, int train) {
    std::vector<long long> pos;
    for (int i = 1; i <= train; ++i) {
        pos.push_back(-(static_cast<long long>(guard) + i));
        pos.push_back(static_cast<long long>(guard) + i);
    }
    double s = 0.0;
    for (long long a : pos)
        for (long long b : pos) {
            const std::size_t d = static_cast<std::size_t>(std::llabs(a - b));
            if (d < c.size()) s += c[d];
        }
    const double n = static_cast<double>(pos.size());
    return (s > 0.0) ? std::max(1.0, n * n / s) : n;
}

}  // namespace radar
