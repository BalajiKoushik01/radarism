#pragma once
#include <cmath>
#include <cstdint>
#include <numbers>
#include <random>
#include <stdexcept>
#include <volk/volk.h>
#include "radar/fft.hpp"
#include "radar/types.hpp"

namespace radar {

// Leeson-shaped SSB phase noise:
//   L(f) = floor * (1 + (f_leeson/f)^2) * (1 + f_flicker/f)       [linear, 1/Hz]
// One-sided phase PSD: S_phi(f) = 2 * L(f)  [rad^2/Hz].
// Set a corner to 0 to disable that term.
struct PhaseNoiseParams {
    double floor_dbc_hz = -140.0;
    double leeson_corner_hz = 1e3;    // f0 / (2 Q_loaded)
    double flicker_corner_hz = 100.0;
    std::uint64_t seed = 1;           // fixed seed => reproducible runs and deterministic tests
};

inline double leeson_psd(const PhaseNoiseParams& p, double f) {  // f > 0, returns rad^2/Hz
    f = std::fabs(f);
    double shape = 1.0;
    if (p.leeson_corner_hz > 0.0) { const double r = p.leeson_corner_hz / f; shape *= 1.0 + r * r; }
    if (p.flicker_corner_hz > 0.0) shape *= 1.0 + p.flicker_corner_hz / f;
    return 2.0 * std::pow(10.0, p.floor_dbc_hz / 10.0) * shape;
}

// Frequency-domain synthesis: shape complex Gaussian noise by sqrt(PSD), inverse FFT, take the real part.
// Bins k > N/2 are mapped to NEGATIVE frequencies (k - N). Offsets below Fs/N cannot be represented
// (DC bin is zeroed, so the mean phase is removed).
inline RBuffer generate_phase_noise(std::size_t n, double fs, const PhaseNoiseParams& p) {
    if (n < 2 || !(fs > 0.0)) throw std::invalid_argument("generate_phase_noise: bad arguments");
    CBuffer X(n), z(n);
    std::mt19937_64 rng(p.seed);
    std::normal_distribution<float> g(0.0f, 1.0f);
    const double df = fs / static_cast<double>(n);
    X[0] = {0.0f, 0.0f};
    for (std::size_t k = 1; k < n; ++k) {
        const long long m = (k <= n / 2) ? static_cast<long long>(k) : static_cast<long long>(k) - static_cast<long long>(n);
        const double c = std::sqrt(leeson_psd(p, static_cast<double>(m) * df) * df / 2.0);
        const float a = g(rng), b = g(rng);
        X[k] = {static_cast<float>(c) * a, static_cast<float>(c) * b};
    }
    Fft ifft(n, false);  // unnormalized backward FFT: exactly the scaling derived for c_k
    ifft.execute(X.data(), z.data());
    RBuffer phi(n);
    for (std::size_t i = 0; i < n; ++i) phi[i] = z[i].real();
    return phi;
}

// x[n] *= exp(j*phi[n])  (multiplicative, as the spec requires)
inline void apply_phase_noise(CBuffer& x, double fs, const PhaseNoiseParams& p) {
    const std::size_t n = x.size();
    RBuffer phi = generate_phase_noise(n, fs, p);
    CBuffer ph(n), out(n);
    for (std::size_t i = 0; i < n; ++i) ph[i] = std::polar(1.0f, phi[i]);
    volk_32fc_x2_multiply_32fc(out.data(), x.data(), ph.data(), static_cast<unsigned int>(n));
    x.swap(out);
}

}  // namespace radar
