#pragma once
#include <atomic>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include "radar/cfar.hpp"
#include "radar/cfar_neff.hpp"
#include "radar/channel.hpp"
#include "radar/impairments.hpp"
#include "radar/matched_filter.hpp"
#include "radar/peak_interp.hpp"
#include "radar/phase_noise.hpp"
#include "radar/targets.hpp"
#include "radar/waveforms.hpp"

namespace radar {

struct SimulationParams {
    WaveformParams wf;
    ChannelParams ch;
    PhaseNoiseParams pn;
    bool phase_noise_enabled = false;
    CfarParams cfar;
    bool run_detection = true;
    MfWeighting mf_weighting = MfWeighting::None;   // taper on the receive replica
    RxImpairments rx_imp;                           // receiver impairments (default: ideal)
    std::vector<TargetSpec> extra_targets;          // additional point targets (up to 8)
    bool cfar_correlation_aware = false;            // size the CFAR threshold for correlated (oversampled) noise
};

struct SimulationResult {
    SimulationParams params;
    ChannelReport report;
    std::size_t n_pulse = 0;
    CBuffer tx;                        // clean, ideal Tx (this is what an AWG should play)
    CBuffer tx_impaired;               // Tx after LO phase noise (empty unless phase noise is enabled)
    CBuffer rx;                        // echo (from phase-noisy Tx if enabled) + AWGN
    bool has_range_profile = false;    // false for CW
    RBuffer mf_power;                  // |matched filter output|^2 (peak = Pr, floor = Pn/E)
    CfarResult cfar;
    std::vector<double> detection_range_m;   // from the refined (sub-bin) peak position
    std::vector<double> detection_bin;       // refined peak position in bins
    double cfar_effective_cells = 0.0;       // independent-cell equivalent of the CFAR training window (0 = assumed iid)
    ImpairmentReport imp;                    // what the receiver impairments did
    double mf_noise_floor_w = 0.0;           // noise floor of mf_power; the echo peak equals Pr
};

struct Cancelled : std::runtime_error {
    Cancelled() : std::runtime_error("simulation cancelled") {}
};

using ResultPtr = std::shared_ptr<const SimulationResult>;

inline const char* waveform_name(WaveformType t) {
    switch (t) {
        case WaveformType::CW: return "CW";
        case WaveformType::Pulse: return "Pulse";
        case WaveformType::LFM: return "LFM";
        case WaveformType::Barker: return "Barker";
        case WaveformType::Frank: return "Frank";
    }
    return "?";
}

// "" if everything is valid, else a UI-ready message
// Checks the additional targets: each must satisfy the channel limits for as long as it keeps moving (span_s).
inline std::string validate_extra_targets(const SimulationParams& p, double span_s = 0.0) {
    if (p.extra_targets.size() > 8) return "At most 8 additional targets are supported.";
    for (std::size_t i = 0; i < p.extra_targets.size(); ++i) {
        const TargetSpec& t = p.extra_targets[i];
        const std::string who = "Target " + std::to_string(i + 2) + ": ";
        if (!(std::fabs(t.rel_db) <= 100.0)) return who + "relative level must be within +/-100 dB.";
        ChannelParams ce = p.ch;
        ce.range_m = t.range_m;
        ce.velocity_mps = t.velocity_mps;
        if (auto e = validate_channel(p.wf, ce); !e.empty()) return who + e;
        const double r_end = t.range_m - t.velocity_mps * span_s;
        if (!(r_end > 0.0) || r_end > 1e9) return who + "range leaves (0, 1e9] m during the CPI.";
    }
    return "";
}

inline std::string validate_all(const SimulationParams& p) {
    if (auto e = validate(p.wf); !e.empty()) return e;
    if (auto e = validate_channel(p.wf, p.ch); !e.empty()) return e;
    if (auto e = validate_impairments(p.rx_imp); !e.empty()) return e;
    if (auto e = validate_extra_targets(p); !e.empty()) return e;
    if (p.run_detection && p.wf.type != WaveformType::CW)
        if (auto e = validate_cfar(p.cfar, total_samples(p.wf)); !e.empty()) return e;
    return "";
}

inline void add_extra_targets(CBuffer& rx, const Waveform& tx, const SimulationParams& p) {
    for (const TargetSpec& t : p.extra_targets) {
        ChannelParams ce = p.ch;
        ce.range_m = t.range_m;
        ce.velocity_mps = t.velocity_mps;
        ChannelReport re = analyze(p.wf, ce);
        re.pr_w *= db_to_lin(t.rel_db);
        const CBuffer e = make_echo(tx, p.wf, ce, re);
        for (std::size_t i = 0; i < rx.size(); ++i) rx[i] += e[i];
    }
}

inline ResultPtr simulate(const SimulationParams& p, const std::atomic<bool>* cancel = nullptr) {
    if (auto e = validate_all(p); !e.empty()) throw std::invalid_argument(e);
    auto check = [&] { if (cancel && cancel->load(std::memory_order_relaxed)) throw Cancelled(); };

    auto res = std::make_shared<SimulationResult>();
    res->params = p;

    Waveform tx = generate(p.wf);
    res->n_pulse = tx.n_pulse;
    res->tx = tx.samples;                       // clean copy
    check();

    if (p.phase_noise_enabled) {
        apply_phase_noise(tx.samples, p.wf.fs, p.pn);   // impairs the working copy only
        res->tx_impaired = tx.samples;                   // kept for the PSD / analysis views
    }
    check();

    res->report = analyze(p.wf, p.ch);
    res->report.mf_gain_db = lin_to_db(static_cast<double>(tx.n_pulse));
    CBuffer rx = make_echo(tx, p.wf, p.ch, res->report);
    add_extra_targets(rx, tx, p);
    tx.samples = CBuffer();                     // release the impaired copy
    check();

    add_awgn(rx, res->report.pn_w, p.ch.seed);
    res->imp = apply_rx_impairments(rx, p.rx_imp, std::sqrt(res->report.pn_w / 2.0));
    res->rx = std::move(rx);
    check();

    if (p.wf.type != WaveformType::CW) {
        double neff = 0.0;
        {
            const CBuffer h = weighted_replica(res->tx, res->n_pulse, p.mf_weighting);   // clean (tapered) replica
            double coherent = 0.0, eh = 0.0;
            for (std::size_t i = 0; i < res->n_pulse; ++i) {
                coherent += static_cast<double>(std::real(std::conj(h[i]) * res->tx[i]));
                eh += static_cast<double>(std::norm(h[i]));
            }
            CBuffer y = matched_filter(res->rx, h, coherent);   // unity gain on the signal
            res->mf_power = magnitude_squared(y);
            res->mf_noise_floor_w = res->report.pn_w * eh / (coherent * coherent);
            res->report.mf_gain_db = lin_to_db(coherent * coherent / eh);
            if (p.run_detection && p.cfar_correlation_aware) {
                const auto c = noise_correlation_power(h, res->n_pulse,
                                                       static_cast<std::size_t>(2 * (p.cfar.guard_cells + p.cfar.train_cells)));
                neff = cfar_effective_cells(c, p.cfar.guard_cells, p.cfar.train_cells);
            }
        }
        res->has_range_profile = true;
        check();
        if (p.run_detection) {
            CfarParams cp = p.cfar;
            cp.effective_cells = neff;
            res->cfar_effective_cells = neff;
            res->cfar = ca_cfar(res->mf_power, cp);
            const double nbins = static_cast<double>(res->mf_power.size());
            for (auto d : res->cfar.detections) {
                const double b = std::fmod(static_cast<double>(d) + peak_offset(res->mf_power, d) + nbins, nbins);
                res->detection_bin.push_back(b);
                res->detection_range_m.push_back(range_from_bin(b, p.wf.fs));
            }
        }
    }
    return res;
}

}  // namespace radar
