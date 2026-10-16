#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include "radar/types.hpp"
#include "radar/waveforms.hpp"

namespace radar {

inline constexpr double kSpeedOfLight = 299792458.0;
inline constexpr double kBoltzmann = 1.380649e-23;
inline constexpr double kRefTemp = 290.0;

enum class SnrMode { Manual, Physical };

struct ChannelParams {
    double fc = 10e9;             // carrier frequency [Hz]
    double range_m = 10e3;        // target range [m]
    double velocity_mps = 0.0;    // radial velocity, positive = closing
    SnrMode snr_mode = SnrMode::Manual;
    double snr_db = 13.0;         // manual mode: per-sample input SNR
    double pt_w = 1e3;            // physical mode parameters
    double gt_db = 30.0;
    double gr_db = 30.0;
    double rcs_m2 = 1.0;
    double loss_db = 0.0;
    double noise_figure_db = 3.0; // used by both modes (sets the noise floor)
    std::uint64_t seed = 12345;
    int delay_half_taps = 16;     // fractional-delay filter has 2*this taps
};

struct ChannelReport {
    double td_s = 0.0;
    double delay_samples = 0.0;
    double doppler_hz = 0.0;
    double r_unambiguous_m = 0.0;
    double r_blind_m = 0.0;
    double pr_w = 0.0;            // received pulse power
    double pn_w = 0.0;            // noise power per sample
    double snr_in_db = 0.0;       // per-sample input SNR
    double mf_gain_db = 0.0;      // 10*log10(n_pulse)
    bool in_blind_zone = false;   // echo overlaps the transmit pulse
    bool straddles_pri = false;   // echo tail runs past the PRI boundary (eclipsed)
    bool second_time_around = false;
};

inline double db_to_lin(double db) { return std::pow(10.0, db / 10.0); }
inline double lin_to_db(double x) { return 10.0 * std::log10(x); }

inline double thermal_noise_power(double fs, double nf_db) {
    return kBoltzmann * kRefTemp * fs * db_to_lin(nf_db);
}

inline double radar_equation_power(const ChannelParams& c) {
    const double lambda = kSpeedOfLight / c.fc;
    const double four_pi_cubed = std::pow(4.0 * std::numbers::pi, 3.0);
    return c.pt_w * db_to_lin(c.gt_db) * db_to_lin(c.gr_db) * lambda * lambda * c.rcs_m2 /
           (four_pi_cubed * std::pow(c.range_m, 4.0) * db_to_lin(c.loss_db));
}

// "" if valid, else a UI-ready message
inline std::string validate_channel(const WaveformParams& wp, const ChannelParams& c) {
    if (!(c.fc > 0.0)) return "Carrier frequency must be > 0.";
    if (!(c.range_m > 0.0) || c.range_m > 1e9) return "Target range must be in (0, 1e9] m.";
    if (c.delay_half_taps < 4 || c.delay_half_taps > 64) return "Delay filter half-taps must be 4..64.";
    if (c.snr_mode == SnrMode::Physical) {
        if (!(c.pt_w > 0.0) || !(c.rcs_m2 > 0.0)) return "Pt and RCS must be > 0.";
    }
    const double fd = 2.0 * c.velocity_mps * c.fc / kSpeedOfLight;
    const double occupied = (wp.type == WaveformType::LFM) ? 0.5 * wp.bandwidth : 0.0;
    if (std::fabs(fd) + occupied > 0.45 * wp.fs)
        return "Doppler shift plus signal bandwidth exceeds 0.45 x Fs; raise Fs or lower velocity.";
    return "";
}

// Pure bookkeeping (no signals): the UI can call this live while the user edits fields.
inline ChannelReport analyze(const WaveformParams& wp, const ChannelParams& c) {
    ChannelReport r;
    const bool cw = (wp.type == WaveformType::CW);
    r.td_s = 2.0 * c.range_m / kSpeedOfLight;
    r.delay_samples = r.td_s * wp.fs;
    r.doppler_hz = 2.0 * c.velocity_mps * c.fc / kSpeedOfLight;
    r.r_unambiguous_m = kSpeedOfLight * wp.pri / 2.0;
    r.r_blind_m = cw ? 0.0 : kSpeedOfLight * wp.tau / 2.0;
    r.pn_w = thermal_noise_power(wp.fs, c.noise_figure_db);
    r.pr_w = (c.snr_mode == SnrMode::Physical) ? radar_equation_power(c)
                                               : r.pn_w * db_to_lin(c.snr_db);
    r.snr_in_db = lin_to_db(r.pr_w / r.pn_w);
    const double n_pulse = cw ? static_cast<double>(total_samples(wp))
                              : static_cast<double>(std::llround(wp.tau * wp.fs));
    r.mf_gain_db = lin_to_db(n_pulse);
    if (!cw) {
        const double td_eff = std::fmod(r.td_s, wp.pri);
        r.second_time_around = r.td_s >= wp.pri;
        r.in_blind_zone = td_eff < wp.tau;
        r.straddles_pri = td_eff > wp.pri - wp.tau;
    }
    return r;
}

// Noise-free echo: amplitude sqrt(Pr), carrier phase exp(-j2*pi*fc*td), Doppler ramp over the pulse,
// then a windowed-sinc fractional delay, accumulated modulo N (periodic pulse train).
// Apply Tx phase noise to tx.samples BEFORE calling this.
inline CBuffer make_echo(const Waveform& tx, const WaveformParams& wp, const ChannelParams& c,
                         const ChannelReport& rep) {
    const std::size_t N = tx.samples.size();
    const std::size_t np = tx.n_pulse;
    const double two_pi = 2.0 * std::numbers::pi;
    CBuffer echo(N, cf32{0.0f, 0.0f});

    const double amp = std::sqrt(rep.pr_w);
    const double carrier = -two_pi * std::fmod(c.fc * rep.td_s, 1.0);
    const double dop_w = two_pi * rep.doppler_hz / wp.fs;
    CBuffer pulse(np);
    for (std::size_t m = 0; m < np; ++m) {
        const double ph = carrier + dop_w * static_cast<double>(m);
        const cf32 g{static_cast<float>(amp * std::cos(ph)), static_cast<float>(amp * std::sin(ph))};
        pulse[m] = tx.samples[m] * g;
    }

    const int M = c.delay_half_taps;
    const double fl = std::floor(rep.delay_samples);
    const double frac = rep.delay_samples - fl;
    const long long n_int = static_cast<long long>(fl);

    std::vector<double> h(static_cast<std::size_t>(2 * M));
    double sum = 0.0;
    for (int j = 0; j < 2 * M; ++j) {
        const int k = j - M + 1;                    // tap offsets -M+1 .. M
        const double x = static_cast<double>(k) - frac;
        const double sinc = (std::fabs(x) < 1e-12) ? 1.0 : std::sin(std::numbers::pi * x) / (std::numbers::pi * x);
        const double u = x / M;
        const double w = 0.42 + 0.5 * std::cos(std::numbers::pi * u) + 0.08 * std::cos(two_pi * u);
        h[static_cast<std::size_t>(j)] = sinc * w;
        sum += h[static_cast<std::size_t>(j)];
    }
    const long long Nll = static_cast<long long>(N);
    for (int j = 0; j < 2 * M; ++j) {
        const int k = j - M + 1;
        const float tap = static_cast<float>(h[static_cast<std::size_t>(j)] / sum);
        long long b = (n_int + k) % Nll;
        if (b < 0) b += Nll;
        const std::size_t base = static_cast<std::size_t>(b);
        const std::size_t first = std::min(np, N - base);
        for (std::size_t m = 0; m < first; ++m) echo[base + m] += pulse[m] * tap;
        for (std::size_t m = first; m < np; ++m) echo[m - first] += pulse[m] * tap;
    }
    return echo;
}

// Complex AWGN: each of I and Q has variance noise_power/2.
inline void add_awgn(CBuffer& x, double noise_power_w, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<float> g(0.0f, 1.0f);
    const float sigma = static_cast<float>(std::sqrt(noise_power_w / 2.0));
    for (auto& s : x) {
        const float a = g(rng);
        const float b = g(rng);
        s += cf32{a * sigma, b * sigma};
    }
}

struct RxResult {
    CBuffer echo;   // noise-free (kept for now; we will drop it in the memory-budget step)
    CBuffer rx;     // echo + AWGN over the full PRI
    ChannelReport report;
};

inline RxResult simulate_rx(const Waveform& tx, const WaveformParams& wp, const ChannelParams& c) {
    if (auto e = validate_channel(wp, c); !e.empty()) throw std::invalid_argument(e);
    RxResult r;
    r.report = analyze(wp, c);
    r.report.mf_gain_db = lin_to_db(static_cast<double>(tx.n_pulse));
    r.echo = make_echo(tx, wp, c, r.report);
    r.rx = r.echo;
    add_awgn(r.rx, r.report.pn_w, c.seed);
    return r;
}

}  // namespace radar
