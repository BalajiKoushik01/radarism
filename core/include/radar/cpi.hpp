#pragma once
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include "radar/cfar2d.hpp"
#include "radar/fft.hpp"
#include "radar/range_compressor.hpp"
#include "radar/simulation.hpp"

namespace radar {

using SlowWindow = MfWeighting;   // None = rectangular, Hamming, Hann (symmetric, across the M pulses)

struct CpiParams {
    int pulses = 64;
    int gate_len = 512;               // range bins kept per pulse (clamped to N)
    long long gate_start = -1;        // first absolute range bin; -1 = auto (centred on the target at mid-CPI)
    SlowWindow window = SlowWindow::Hann;
    int doppler_fft = 0;              // 0 = smallest power of two >= pulses
    Cfar2dParams cfar;
    bool run_detection = true;
};

struct CpiReport {
    double prf_hz = 0.0, cpi_s = 0.0;
    double doppler_true_hz = 0.0, doppler_apparent_hz = 0.0, velocity_apparent_mps = 0.0;
    double v_unambiguous_mps = 0.0;                                   // +/- this value
    double doppler_resolution_hz = 0.0, velocity_resolution_mps = 0.0; // PRF/M and lambda*PRF/(2M)
    double migration_bins = 0.0;                                      // range migration over the CPI
    double integration_gain_db = 0.0;                                 // coherent slow-time gain incl. window loss
    bool doppler_ambiguous = false;
    bool target_leaves_gate = false;
};

struct CpiDetection {
    std::size_t doppler_bin = 0, gate_bin = 0, range_bin = 0;   // range_bin = (gate_start + gate_bin) % N
    double doppler_hz = 0.0, velocity_mps = 0.0, range_m = 0.0;
    double range_refined_m = 0.0, doppler_refined_hz = 0.0, velocity_refined_mps = 0.0;   // sub-cell estimates
    double power_w = 0.0, snr_db = 0.0;                         // power relative to the RD noise floor
};

struct CpiResult {
    SimulationParams sim;
    CpiParams cpi;
    CpiReport report;
    ChannelReport first_pulse;
    std::size_t n = 0, gate_start = 0, gate_len = 0, doppler_rows = 0;
    RBuffer rd_power;              // [doppler_row * gate_len + gate_bin], fftshifted: row L/2 = zero Doppler
    double rd_noise_floor_w = 0.0; // expected mean noise power per RD cell
    Cfar2dResult cfar;
    ImpairmentReport imp;
    std::vector<CpiDetection> detections;   // strongest first
};
using CpiResultPtr = std::shared_ptr<const CpiResult>;

inline std::size_t cpi_doppler_rows(const CpiParams& c) {
    if (c.doppler_fft > 0) return static_cast<std::size_t>(c.doppler_fft);
    std::size_t L = 1;
    while (L < static_cast<std::size_t>(std::max(c.pulses, 1))) L <<= 1;
    return L;
}

inline CpiReport analyze_cpi(const SimulationParams& p, const CpiParams& c) {
    CpiReport r;
    const double lambda = kSpeedOfLight / p.ch.fc;
    const std::size_t N = total_samples(p.wf);
    r.prf_hz = 1.0 / p.wf.pri;
    r.cpi_s = static_cast<double>(c.pulses) * p.wf.pri;
    r.doppler_true_hz = 2.0 * p.ch.velocity_mps / lambda;
    double f = std::fmod(r.doppler_true_hz + 0.5 * r.prf_hz, r.prf_hz);
    if (f < 0.0) f += r.prf_hz;
    r.doppler_apparent_hz = f - 0.5 * r.prf_hz;
    r.velocity_apparent_mps = r.doppler_apparent_hz * lambda / 2.0;
    r.v_unambiguous_mps = lambda * r.prf_hz / 4.0;
    r.doppler_resolution_hz = r.prf_hz / static_cast<double>(c.pulses);
    r.velocity_resolution_mps = lambda * r.doppler_resolution_hz / 2.0;
    r.doppler_ambiguous = std::fabs(r.doppler_true_hz) >= 0.5 * r.prf_hz;
    r.migration_bins = std::fabs(p.ch.velocity_mps) * static_cast<double>(c.pulses - 1) * p.wf.pri * 2.0 * p.wf.fs / kSpeedOfLight;
    const double G = static_cast<double>(std::min<std::size_t>(static_cast<std::size_t>(std::max(c.gate_len, 0)), N));
    r.target_leaves_gate = (G < static_cast<double>(N)) && (r.migration_bins > 0.8 * G);
    double s1 = 0.0, s2 = 0.0;
    for (int k = 0; k < c.pulses; ++k) {
        const double w = taper_value(c.window, static_cast<std::size_t>(k), static_cast<std::size_t>(c.pulses));
        s1 += w;
        s2 += w * w;
    }
    if (s2 > 0.0) r.integration_gain_db = 10.0 * std::log10(s1 * s1 / s2);
    return r;
}

// "" if valid, else a UI-ready message
inline std::string validate_cpi(const SimulationParams& p, const CpiParams& c) {
    if (auto e = validate(p.wf); !e.empty()) return e;
    if (p.wf.type == WaveformType::CW) return "CPI mode needs a pulsed waveform (CW has no range gate).";
    if (p.phase_noise_enabled) return "Tx phase noise is not supported in CPI mode yet (pulse-to-pulse correlation is not modelled).";
    if (c.pulses < 2 || c.pulses > 4096) return "Pulses per CPI must be 2..4096.";
    if (c.gate_len < 8) return "Range gate must be at least 8 bins.";
    if (c.doppler_fft != 0 && (c.doppler_fft < c.pulses || c.doppler_fft > 65536))
        return "Doppler FFT size must be >= the pulse count and <= 65536.";
    const std::size_t N = total_samples(p.wf);
    const std::size_t G = std::min<std::size_t>(static_cast<std::size_t>(c.gate_len), N);
    const std::size_t L = cpi_doppler_rows(c);
    if (G * static_cast<std::size_t>(c.pulses) > 64'000'000 || G * L > 64'000'000)
        return "CPI too large: pulses x gate (or Doppler FFT x gate) exceeds 64M cells.";
    if (auto e = validate_channel(p.wf, p.ch); !e.empty()) return e;
    if (auto e = validate_impairments(p.rx_imp); !e.empty()) return e;
    if (auto e = validate_extra_targets(p, static_cast<double>(c.pulses - 1) * p.wf.pri); !e.empty()) return e;
    const double r_last = p.ch.range_m - p.ch.velocity_mps * static_cast<double>(c.pulses - 1) * p.wf.pri;
    if (!(r_last > 0.0) || r_last > 1e9) return "Target range leaves (0, 1e9] m during the CPI; lower the velocity or the pulse count.";
    if (c.run_detection)
        if (auto e = validate_cfar2d(c.cfar); !e.empty()) return e;
    return "";
}

inline CpiResultPtr simulate_cpi(const SimulationParams& p, const CpiParams& cp,
                                 const std::atomic<bool>* cancel = nullptr,
                                 const std::function<void(int, int)>& progress = {}) {
    if (auto e = validate_cpi(p, cp); !e.empty()) throw std::invalid_argument(e);
    auto check = [&] { if (cancel && cancel->load(std::memory_order_relaxed)) throw Cancelled(); };

    auto res = std::make_shared<CpiResult>();
    res->sim = p;
    res->cpi = cp;
    const std::size_t N = total_samples(p.wf);
    const std::size_t M = static_cast<std::size_t>(cp.pulses);
    const std::size_t G = std::min<std::size_t>(static_cast<std::size_t>(cp.gate_len), N);
    const std::size_t L = cpi_doppler_rows(cp);
    const double fs = p.wf.fs;
    res->n = N;
    res->gate_len = G;
    res->doppler_rows = L;
    res->report = analyze_cpi(p, cp);
    res->first_pulse = analyze(p.wf, p.ch);

    // Clean, optionally tapered replica; unity gain on the signal
    Waveform tx = generate(p.wf);
    const CBuffer h = weighted_replica(tx.samples, tx.n_pulse, p.mf_weighting);
    double coherent = 0.0, eh = 0.0;
    for (std::size_t i = 0; i < tx.n_pulse; ++i) {
        coherent += static_cast<double>(std::real(std::conj(h[i]) * tx.samples[i]));
        eh += static_cast<double>(std::norm(h[i]));
    }
    RangeCompressor rc(h, coherent);
    const double mf_floor = res->first_pulse.pn_w * eh / (coherent * coherent);

    // Range gate, centred on the target at mid-CPI unless the caller fixed it
    const double r_last = p.ch.range_m - p.ch.velocity_mps * static_cast<double>(M - 1) * p.wf.pri;
    long long gs = cp.gate_start;
    if (gs < 0) {
        const double r_mid = 0.5 * (p.ch.range_m + r_last);
        gs = static_cast<long long>(std::llround(2.0 * r_mid / kSpeedOfLight * fs)) - static_cast<long long>(G / 2);
    }
    const long long Nll = static_cast<long long>(N);
    gs %= Nll;
    if (gs < 0) gs += Nll;
    res->gate_start = static_cast<std::size_t>(gs);
    if (G < N) {
        auto rel = [&](double range) {
            double x = std::fmod(2.0 * range / kSpeedOfLight * fs - static_cast<double>(gs), static_cast<double>(N));
            if (x < 0.0) x += static_cast<double>(N);
            return x;
        };
        const double a = rel(p.ch.range_m), b = rel(r_last), g = static_cast<double>(G);
        res->report.target_leaves_gate = (a < 2.0 || a > g - 3.0 || b < 2.0 || b > g - 3.0);
    }

    // Fast time: one pulse at a time, keep only the gate
    CBuffer slow(G * M);   // [gate_bin * M + pulse]
    CBuffer y;
    const std::uint64_t mix = 0x9E3779B97F4A7C15ULL;
    for (std::size_t k = 0; k < M; ++k) {
        check();
        ChannelParams ck = p.ch;
        ck.range_m = p.ch.range_m - p.ch.velocity_mps * static_cast<double>(k) * p.wf.pri;
        const ChannelReport rk = analyze(p.wf, ck);
        CBuffer rx = make_echo(tx, p.wf, ck, rk);
        for (const TargetSpec& t : p.extra_targets) {
            ChannelParams ce = p.ch;
            ce.range_m = t.range_m - t.velocity_mps * static_cast<double>(k) * p.wf.pri;
            ce.velocity_mps = t.velocity_mps;
            ChannelReport re = analyze(p.wf, ce);
            re.pr_w *= db_to_lin(t.rel_db);
            const CBuffer e = make_echo(tx, p.wf, ce, re);
            for (std::size_t i = 0; i < N; ++i) rx[i] += e[i];
        }
        add_awgn(rx, rk.pn_w, p.ch.seed + mix * (static_cast<std::uint64_t>(k) + 1));
        {
            const ImpairmentReport ir = apply_rx_impairments(rx, p.rx_imp, std::sqrt(rk.pn_w / 2.0), k == 0 ? 0.0 : res->imp.full_scale);
            if (k == 0) res->imp = ir; else res->imp.clipped += ir.clipped;
        }
        rc.compress(rx, y);
        for (std::size_t g = 0; g < G; ++g) slow[g * M + k] = y[(res->gate_start + g) % N];
        if (progress) progress(static_cast<int>(k + 1), static_cast<int>(M));
    }

    // Slow time: window, zero-padded FFT per range bin, fftshift. Power is normalized by (sum w)^2 so that an
    // on-bin target's peak power equals the single-pulse compressed peak power.
    std::vector<double> w(M);
    double s1 = 0.0, s2 = 0.0;
    for (std::size_t k = 0; k < M; ++k) {
        w[k] = taper_value(cp.window, k, M);
        s1 += w[k];
        s2 += w[k] * w[k];
    }
    const double norm = 1.0 / (s1 * s1);
    res->rd_noise_floor_w = mf_floor * s2 / (s1 * s1);
    res->report.integration_gain_db = 10.0 * std::log10(s1 * s1 / s2);
    res->rd_power.assign(L * G, 0.0f);
    Fft sfft(L, true);
    CBuffer row(L), spec(L);
    for (std::size_t g = 0; g < G; ++g) {
        check();
        std::fill(row.begin(), row.end(), cf32{0.0f, 0.0f});
        for (std::size_t k = 0; k < M; ++k) row[k] = slow[g * M + k] * static_cast<float>(w[k]);
        sfft.execute(row.data(), spec.data());
        for (std::size_t f = 0; f < L; ++f) {
            const std::size_t d = (f + L / 2) % L;
            res->rd_power[d * G + g] = static_cast<float>(static_cast<double>(std::norm(spec[f])) * norm);
        }
    }

    if (cp.run_detection) {
        res->cfar = ca_cfar_2d(res->rd_power, L, G, cp.cfar);
        const double lambda = kSpeedOfLight / p.ch.fc;
        for (const auto& d : res->cfar.detections) {
            CpiDetection cd;
            cd.doppler_bin = d.row;
            cd.gate_bin = d.col;
            cd.range_bin = (res->gate_start + d.col) % N;
            cd.doppler_hz = (static_cast<double>(d.row) - static_cast<double>(L / 2)) * res->report.prf_hz / static_cast<double>(L);
            cd.velocity_mps = cd.doppler_hz * lambda / 2.0;
            cd.range_m = range_from_bin(static_cast<double>(cd.range_bin), fs);
            cd.power_w = static_cast<double>(d.power);
            cd.snr_db = 10.0 * std::log10(cd.power_w / res->rd_noise_floor_w);
            {
                auto P = [&](std::size_t rr, std::size_t cc) { return static_cast<double>(res->rd_power[rr * G + cc]); };
                double dr = 0.0;
                if (d.col > 0 && d.col + 1 < G) dr = parabolic_offset(P(d.row, d.col - 1), P(d.row, d.col), P(d.row, d.col + 1));
                const double dd = parabolic_offset(P((d.row + L - 1) % L, d.col), P(d.row, d.col), P((d.row + 1) % L, d.col));
                cd.range_refined_m = range_from_bin(static_cast<double>(cd.range_bin) + dr, fs);
                cd.doppler_refined_hz = (static_cast<double>(d.row) + dd - static_cast<double>(L / 2)) * res->report.prf_hz / static_cast<double>(L);
                cd.velocity_refined_mps = cd.doppler_refined_hz * lambda / 2.0;
            }
            res->detections.push_back(cd);
        }
    }
    return res;
}

}  // namespace radar
