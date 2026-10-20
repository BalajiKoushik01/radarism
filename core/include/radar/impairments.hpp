#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include "radar/types.hpp"

namespace radar {

// Receiver impairments applied to the noisy Rx after AWGN, in this order: I/Q imbalance, DC offset, ADC.
// Everything defaults to off; then the pipeline output is unchanged.
struct RxImpairments {
    int adc_bits = 0;            // 0 = ideal converter; otherwise 2..16
    double adc_fs_sigma = 0.0;   // ADC full scale in noise sigmas per I/Q component; 0 = auto (no clipping)
    double iq_gain_db = 0.0;     // Q-branch gain error relative to I
    double iq_phase_deg = 0.0;   // Q-branch phase error
    double dc_i_sigma = 0.0;     // DC offsets in noise sigmas per component
    double dc_q_sigma = 0.0;
    bool active() const {
        return adc_bits != 0 || iq_gain_db != 0.0 || iq_phase_deg != 0.0 || dc_i_sigma != 0.0 || dc_q_sigma != 0.0;
    }
};

struct ImpairmentReport {
    bool applied = false;
    std::size_t clipped = 0;              // I or Q values clipped by the ADC
    double lsb = 0.0, full_scale = 0.0;   // sqrt(W)
    double image_rejection_db = 0.0;      // from the I/Q imbalance (0 when none)
};

inline std::string validate_impairments(const RxImpairments& r) {
    if (r.adc_bits != 0 && (r.adc_bits < 2 || r.adc_bits > 16)) return "ADC resolution must be off (0) or 2..16 bits.";
    if (!(r.adc_fs_sigma >= 0.0 && r.adc_fs_sigma <= 1000.0)) return "ADC full scale must be 0 (auto) up to 1000 sigma.";
    if (!(std::fabs(r.iq_gain_db) <= 20.0)) return "I/Q gain imbalance must be within +/-20 dB.";
    if (!(std::fabs(r.iq_phase_deg) <= 45.0)) return "I/Q phase error must be within +/-45 degrees.";
    if (!(std::fabs(r.dc_i_sigma) <= 100.0 && std::fabs(r.dc_q_sigma) <= 100.0)) return "DC offsets must be within +/-100 sigma.";
    return "";
}

// y = K1 x + K2 conj(x),  K1 = (1 + g e^{j phi})/2 form: see derivation in the tests.  Returns 10 log10(|K1|^2/|K2|^2).
inline double iq_image_rejection_db(double gain_db, double phase_deg) {
    const double g = std::pow(10.0, gain_db / 20.0), phi = phase_deg * std::numbers::pi / 180.0;
    const double k1r = 0.5 * (1.0 + g * std::cos(phi)), k1i = 0.5 * g * std::sin(phi);
    const double k2r = 0.5 * (1.0 - g * std::cos(phi)), k2i = k1i;
    const double p1 = k1r * k1r + k1i * k1i, p2 = k2r * k2r + k2i * k2i;
    return (p2 > 1e-30) ? 10.0 * std::log10(p1 / p2) : 300.0;
}

// sigma = noise standard deviation per I/Q component. fs_override > 0 forces the ADC full scale (used to keep it
// identical across the pulses of a CPI).
inline ImpairmentReport apply_rx_impairments(CBuffer& x, const RxImpairments& r, double sigma, double fs_override = 0.0) {
    ImpairmentReport rep;
    if (!r.active() || x.empty()) return rep;
    rep.applied = true;
    const bool iq = (r.iq_gain_db != 0.0 || r.iq_phase_deg != 0.0);
    const double g = std::pow(10.0, r.iq_gain_db / 20.0), phi = r.iq_phase_deg * std::numbers::pi / 180.0;
    const double cg = g * std::cos(phi), sg = g * std::sin(phi);
    if (iq) rep.image_rejection_db = iq_image_rejection_db(r.iq_gain_db, r.iq_phase_deg);
    const double dci = r.dc_i_sigma * sigma, dcq = r.dc_q_sigma * sigma;
    if (iq || dci != 0.0 || dcq != 0.0) {
        for (auto& s : x) {
            const double i = static_cast<double>(s.real()), q = static_cast<double>(s.imag());
            const double q2 = iq ? cg * q + sg * i : q;   // Q' = g (Q cos(phi) + I sin(phi)); I unchanged
            s = cf32{static_cast<float>(i + dci), static_cast<float>(q2 + dcq)};
        }
    }
    if (r.adc_bits != 0) {
        double fs = fs_override;
        if (!(fs > 0.0)) {
            if (r.adc_fs_sigma > 0.0) {
                fs = r.adc_fs_sigma * sigma;
            } else {
                double pk = 0.0;
                for (const auto& s : x) pk = std::max({pk, static_cast<double>(std::fabs(s.real())), static_cast<double>(std::fabs(s.imag()))});
                fs = std::max(1.05 * pk, 6.0 * sigma);
            }
        }
        if (!(fs > 0.0)) fs = 1.0;
        const double levels = std::ldexp(1.0, r.adc_bits), lsb = 2.0 * fs / levels;
        const double lo = -levels / 2.0, hi = levels / 2.0 - 1.0;
        auto q = [&](double v) {
            double c = std::nearbyint(v / lsb);
            if (c < lo) { c = lo; ++rep.clipped; } else if (c > hi) { c = hi; ++rep.clipped; }
            return static_cast<float>(c * lsb);
        };
        for (auto& s : x) s = cf32{q(static_cast<double>(s.real())), q(static_cast<double>(s.imag()))};
        rep.lsb = lsb;
        rep.full_scale = fs;
    }
    return rep;
}

}  // namespace radar
