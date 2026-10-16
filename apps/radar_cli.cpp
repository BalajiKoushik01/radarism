#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "radar/export_hdf5.hpp"
#include "radar/exporters.hpp"
#include "radar/simulation.hpp"

using namespace radar;

struct Args {
    std::map<std::string, std::string> kv;
    std::set<std::string> used;
    Args(int argc, char** argv) {
        for (int i = 1; i < argc; ++i) {
            const std::string s = argv[i];
            const auto pos = s.find('=');
            if (pos == std::string::npos || pos == 0) throw std::invalid_argument("expected key=value, got: " + s);
            kv[s.substr(0, pos)] = s.substr(pos + 1);
        }
    }
    std::string str(const std::string& k, const std::string& def) {
        used.insert(k);
        auto it = kv.find(k);
        return it == kv.end() ? def : it->second;
    }
    double num(const std::string& k, double def) {
        used.insert(k);
        auto it = kv.find(k);
        if (it == kv.end()) return def;
        char* e = nullptr;
        const double v = std::strtod(it->second.c_str(), &e);
        if (it->second.empty() || *e != '\0') throw std::invalid_argument("bad number for " + k + ": " + it->second);
        return v;
    }
    bool flag(const std::string& k, bool def) {
        const std::string v = str(k, def ? "on" : "off");
        if (v == "on" || v == "1" || v == "true") return true;
        if (v == "off" || v == "0" || v == "false") return false;
        throw std::invalid_argument(k + " must be on or off");
    }
    void check_unused() const {
        for (const auto& [k, v] : kv)
            if (!used.count(k)) throw std::invalid_argument("unknown option: " + k + " (run with 'help')");
    }
};

static void usage() {
    std::puts(
        "radar_cli [key=value ...]   (no arguments = default LFM demo)\n"
        "\n"
        "Waveform : type=cw|pulse|lfm|barker|frank  fs=1e6  pri=10e-3  tau=2e-3  bw=2e5  code=13\n"
        "Channel  : fc=10e9  range=150000  vel=0  nf=3  seed=12345\n"
        "SNR      : mode=manual|physical\n"
        "           manual   : snr=-15        (per-sample input SNR, dB)\n"
        "           physical : pt=1000 gt=30 gr=30 rcs=1 loss=0\n"
        "Phase nse: pn=on|off  pnfloor=-140  pncorner=1000  pnflicker=100  pnseed=1\n"
        "CFAR     : detect=on|off  pfa=1e-6  guard=20  train=64\n"
        "MF weight: mfw=none|hamming|hann\n"
        "Export   : out=PREFIX  fmt=float|int16  csv=on|off\n"
        "           writes PREFIX.h5, PREFIX_tx.bin(+.json), PREFIX_rx.bin(+.json), PREFIX_rx.csv (if N<=100000)\n"
        "\n"
        "Example: ./build/radar_cli type=lfm range=300000 snr=-18 pn=on out=/tmp/run1");
}

int main(int argc, char** argv) {
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string s = argv[i];
            if (s == "help" || s == "-h" || s == "--help") { usage(); return 0; }
        }
        Args a(argc, argv);
        SimulationParams p;

        const std::string type = a.str("type", "lfm");
        if (type == "cw") p.wf.type = WaveformType::CW;
        else if (type == "pulse") p.wf.type = WaveformType::Pulse;
        else if (type == "lfm") p.wf.type = WaveformType::LFM;
        else if (type == "barker") p.wf.type = WaveformType::Barker;
        else if (type == "frank") p.wf.type = WaveformType::Frank;
        else throw std::invalid_argument("type must be cw, pulse, lfm, barker or frank");
        p.wf.fs = a.num("fs", 1e6);
        p.wf.pri = a.num("pri", 10e-3);
        p.wf.tau = a.num("tau", 2e-3);
        p.wf.bandwidth = a.num("bw", 2e5);
        p.wf.code_length = static_cast<int>(a.num("code", 13));

        p.ch.fc = a.num("fc", 10e9);
        p.ch.range_m = a.num("range", 150e3);
        p.ch.velocity_mps = a.num("vel", 0.0);
        p.ch.noise_figure_db = a.num("nf", 3.0);
        p.ch.seed = static_cast<std::uint64_t>(a.num("seed", 12345));
        const std::string mode = a.str("mode", "manual");
        if (mode == "manual") p.ch.snr_mode = SnrMode::Manual;
        else if (mode == "physical") p.ch.snr_mode = SnrMode::Physical;
        else throw std::invalid_argument("mode must be manual or physical");
        p.ch.snr_db = a.num("snr", -15.0);
        p.ch.pt_w = a.num("pt", 1e3);
        p.ch.gt_db = a.num("gt", 30.0);
        p.ch.gr_db = a.num("gr", 30.0);
        p.ch.rcs_m2 = a.num("rcs", 1.0);
        p.ch.loss_db = a.num("loss", 0.0);

        p.phase_noise_enabled = a.flag("pn", false);
        p.pn.floor_dbc_hz = a.num("pnfloor", -140.0);
        p.pn.leeson_corner_hz = a.num("pncorner", 1e3);
        p.pn.flicker_corner_hz = a.num("pnflicker", 100.0);
        p.pn.seed = static_cast<std::uint64_t>(a.num("pnseed", 1));

        p.run_detection = a.flag("detect", true);
        const std::string mfw = a.str("mfw", "none");
        if (mfw == "hamming") p.mf_weighting = MfWeighting::Hamming;
        else if (mfw == "hann") p.mf_weighting = MfWeighting::Hann;
        else if (mfw != "none") throw std::invalid_argument("mfw must be none, hamming or hann");
        p.cfar.pfa = a.num("pfa", 1e-6);
        p.cfar.guard_cells = static_cast<int>(a.num("guard", 20));
        p.cfar.train_cells = static_cast<int>(a.num("train", 64));

        const std::string out = a.str("out", "");
        const std::string fmt = a.str("fmt", "float");
        if (fmt != "float" && fmt != "int16") throw std::invalid_argument("fmt must be float or int16");
        const bool csv_on = a.flag("csv", true);
        a.check_unused();

        const auto t0 = std::chrono::steady_clock::now();
        ResultPtr res = simulate(p);
        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        const SimulationResult& r = *res;
        const auto& rep = r.report;
        const std::size_t N = r.rx.size();

        std::printf("== Simulation (%.1f ms) ==\n", ms);
        std::printf("Waveform      : %s, N = %zu samples, pulse = %zu samples\n", waveform_name(p.wf.type), N, r.n_pulse);
        std::printf("Delay         : %.3f us = %.3f samples, Doppler = %.1f Hz\n", rep.td_s * 1e6, rep.delay_samples, rep.doppler_hz);
        std::printf("R_unambiguous : %.1f km, R_blind: %.1f m\n", rep.r_unambiguous_m / 1e3, rep.r_blind_m);
        std::printf("Pr = %.3e W, Pn/sample = %.3e W\n", rep.pr_w, rep.pn_w);
        std::printf("Input SNR     : %.2f dB, MF gain: +%.2f dB, expected output SNR: %.2f dB\n",
                    rep.snr_in_db, rep.mf_gain_db, rep.snr_in_db + rep.mf_gain_db);
        if (rep.in_blind_zone) std::puts("FLAG: echo overlaps the transmit pulse (blind zone).");
        if (rep.straddles_pri) std::puts("FLAG: echo is eclipsed / truncated by the PRI boundary.");
        if (rep.second_time_around) std::puts("FLAG: second-time-around echo (range is ambiguous).");

        if (!r.has_range_profile) {
            std::puts("CW: no range profile (matched filter disabled).");
        } else {
            std::vector<float> tmp(r.mf_power.begin(), r.mf_power.end());
            std::nth_element(tmp.begin(), tmp.begin() + static_cast<std::ptrdiff_t>(N / 2), tmp.end());
            const double floor_mean = static_cast<double>(tmp[N / 2]) / std::log(2.0);   // exponential: median = mean*ln2
            const auto pk = static_cast<std::size_t>(std::max_element(r.mf_power.begin(), r.mf_power.end()) - r.mf_power.begin());
            const double exp_bin = std::fmod(rep.delay_samples, static_cast<double>(N));
            std::printf("\n== Range profile ==\n");
            std::printf("Expected peak bin: %.2f (apparent range %.3f km)\n", exp_bin, bin_to_range_m(static_cast<std::size_t>(std::llround(exp_bin)) % N, p.wf.fs) / 1e3);
            std::printf("Strongest bin    : %zu (%.3f km), %.2f dB above estimated floor\n", pk, bin_to_range_m(pk, p.wf.fs) / 1e3,
                        10.0 * std::log10(static_cast<double>(r.mf_power[pk]) / floor_mean));
            if (p.run_detection) {
                std::printf("CFAR (Pfa=%.0e, alpha=%.2f): %zu cells over threshold, %zu detection(s)\n", p.cfar.pfa, r.cfar.alpha,
                            r.cfar.cells_over_threshold, r.cfar.detections.size());
                for (std::size_t i = 0; i < r.cfar.detections.size(); ++i) {
                    const auto d = r.cfar.detections[i];
                    std::printf("  #%zu bin %zu  range %.3f km  %.2f dB over threshold\n", i + 1, d, r.detection_range_m[i] / 1e3,
                                10.0 * std::log10(static_cast<double>(r.mf_power[d]) / static_cast<double>(r.cfar.threshold[d])));
                }
            }
        }

        if (!out.empty()) {
            std::puts("\n== Export ==");
            write_hdf5(out + ".h5", r);
            std::printf("wrote %s.h5\n", out.c_str());
            const BinFormat bf = (fmt == "int16") ? BinFormat::Int16 : BinFormat::Float32;
            auto ti = write_bin(out + "_tx.bin", r.tx, bf, p.wf.fs, p.ch.fc, "tx_clean_for_awg");
            auto ri = write_bin(out + "_rx.bin", r.rx, bf, p.wf.fs, p.ch.fc, "rx_for_sdr_dsp");
            std::printf("wrote %s (+.json), scale %.6g\n", ti.bin_path.c_str(), ti.scale);
            std::printf("wrote %s (+.json), scale %.6g\n", ri.bin_path.c_str(), ri.scale);
            if (csv_on) {
                if (N <= kCsvMaxSamples) {
                    write_csv(out + "_rx.csv", r.rx, p.wf.fs);
                    std::printf("wrote %s_rx.csv\n", out.c_str());
                } else {
                    std::printf("CSV skipped: N = %zu exceeds the %zu-sample limit.\n", N, kCsvMaxSamples);
                }
            }
        }
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
