#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include "radar/cpi.hpp"

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
        "radar_cpi [key=value ...]   (no arguments = demo: LFM, target at 30 km closing at 3.7474 m/s = +250 Hz)\n"
        "\n"
        "Waveform : type=pulse|lfm|barker|frank  fs=1e6  pri=1e-3  tau=100e-6  bw=2e5  code=13\n"
        "Channel  : fc=10e9  range=30000  vel=3.7474  nf=3  seed=12345   (vel > 0 = closing)\n"
        "SNR      : mode=manual|physical   snr=-10   pt=1000 gt=30 gr=30 rcs=1 loss=0\n"
        "MF       : mfw=none|hamming|hann\n"
        "CPI      : pulses=64  gate=256  dfft=0 (0 = next power of two)  win=hann|hamming|rect\n"
        "CFAR 2D  : detect=on|off  pfa=1e-6  grange=10  gdop=2  trange=20  tdop=4");
}

int main(int argc, char** argv) {
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string s = argv[i];
            if (s == "help" || s == "-h" || s == "--help") { usage(); return 0; }
        }
        Args a(argc, argv);
        SimulationParams p;
        CpiParams c;
        const std::string type = a.str("type", "lfm");
        if (type == "pulse") p.wf.type = WaveformType::Pulse;
        else if (type == "lfm") p.wf.type = WaveformType::LFM;
        else if (type == "barker") p.wf.type = WaveformType::Barker;
        else if (type == "frank") p.wf.type = WaveformType::Frank;
        else throw std::invalid_argument("type must be pulse, lfm, barker or frank");
        p.wf.fs = a.num("fs", 1e6);
        p.wf.pri = a.num("pri", 1e-3);
        p.wf.tau = a.num("tau", 100e-6);
        p.wf.bandwidth = a.num("bw", 2e5);
        p.wf.code_length = static_cast<int>(a.num("code", 13));
        p.ch.fc = a.num("fc", 10e9);
        p.ch.range_m = a.num("range", 30e3);
        p.ch.velocity_mps = a.num("vel", 3.7474);
        p.ch.noise_figure_db = a.num("nf", 3.0);
        p.ch.seed = static_cast<std::uint64_t>(a.num("seed", 12345));
        const std::string mode = a.str("mode", "manual");
        if (mode == "manual") p.ch.snr_mode = SnrMode::Manual;
        else if (mode == "physical") p.ch.snr_mode = SnrMode::Physical;
        else throw std::invalid_argument("mode must be manual or physical");
        p.ch.snr_db = a.num("snr", -10.0);
        p.ch.pt_w = a.num("pt", 1e3);
        p.ch.gt_db = a.num("gt", 30.0);
        p.ch.gr_db = a.num("gr", 30.0);
        p.ch.rcs_m2 = a.num("rcs", 1.0);
        p.ch.loss_db = a.num("loss", 0.0);
        const std::string mfw = a.str("mfw", "none");
        if (mfw == "hamming") p.mf_weighting = MfWeighting::Hamming;
        else if (mfw == "hann") p.mf_weighting = MfWeighting::Hann;
        else if (mfw != "none") throw std::invalid_argument("mfw must be none, hamming or hann");
        c.pulses = static_cast<int>(a.num("pulses", 64));
        c.gate_len = static_cast<int>(a.num("gate", 256));
        c.doppler_fft = static_cast<int>(a.num("dfft", 0));
        const std::string win = a.str("win", "hann");
        if (win == "hann") c.window = SlowWindow::Hann;
        else if (win == "hamming") c.window = SlowWindow::Hamming;
        else if (win == "rect") c.window = SlowWindow::None;
        else throw std::invalid_argument("win must be hann, hamming or rect");
        c.run_detection = a.flag("detect", true);
        c.cfar.pfa = a.num("pfa", 1e-6);
        c.cfar.guard_range = static_cast<int>(a.num("grange", 10));
        c.cfar.guard_doppler = static_cast<int>(a.num("gdop", 2));
        c.cfar.train_range = static_cast<int>(a.num("trange", 20));
        c.cfar.train_doppler = static_cast<int>(a.num("tdop", 4));
        a.check_unused();

        const auto t0 = std::chrono::steady_clock::now();
        CpiResultPtr res = simulate_cpi(p, c);
        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        const CpiResult& r = *res;
        const auto& rep = r.report;
        const double lambda = kSpeedOfLight / p.ch.fc;
        const std::size_t G = r.gate_len, L = r.doppler_rows;

        std::printf("== CPI (%.1f ms) ==\n", ms);
        std::printf("PRF %.1f Hz | CPI %.2f ms (%d pulses) | N = %zu | gate bins %zu..%zu (len %zu) | Doppler rows %zu\n",
                    rep.prf_hz, rep.cpi_s * 1e3, c.pulses, r.n, r.gate_start, (r.gate_start + G - 1) % r.n, G, L);
        std::printf("Doppler true %.2f Hz, apparent %.2f Hz (%.3f m/s) | unambiguous velocity +/- %.3f m/s\n",
                    rep.doppler_true_hz, rep.doppler_apparent_hz, rep.velocity_apparent_mps, rep.v_unambiguous_mps);
        std::printf("Resolution: %.2f Hz = %.3f m/s | range migration %.3f bins\n",
                    rep.doppler_resolution_hz, rep.velocity_resolution_mps, rep.migration_bins);
        std::printf("Input SNR %.2f dB, MF gain +%.2f dB, slow-time gain +%.2f dB => expected peak %.2f dB above the RD floor\n",
                    r.first_pulse.snr_in_db, r.first_pulse.mf_gain_db, rep.integration_gain_db,
                    10.0 * std::log10(r.first_pulse.pr_w / r.rd_noise_floor_w));
        if (rep.doppler_ambiguous) std::puts("FLAG: Doppler ambiguous (target velocity exceeds +/- v_unambiguous, so it folds).");
        if (rep.target_leaves_gate) std::puts("FLAG: target migrates out of, or to the edge of, the range gate.");
        if (r.first_pulse.in_blind_zone) std::puts("FLAG: echo overlaps the transmit pulse (blind zone).");
        if (r.first_pulse.straddles_pri) std::puts("FLAG: echo is eclipsed by the PRI boundary.");
        if (r.first_pulse.second_time_around) std::puts("FLAG: second-time-around echo (range is ambiguous).");

        const auto pk = static_cast<std::size_t>(std::max_element(r.rd_power.begin(), r.rd_power.end()) - r.rd_power.begin());
        const std::size_t prow = pk / G, pcol = pk % G;
        const double fd = (static_cast<double>(prow) - static_cast<double>(L / 2)) * rep.prf_hz / static_cast<double>(L);
        const std::size_t abs_bin = (r.gate_start + pcol) % r.n;
        std::printf("\nStrongest cell: Doppler row %zu = %.2f Hz = %.3f m/s | range bin %zu = %.3f km | %.2f dB above floor\n",
                    prow, fd, fd * lambda / 2.0, abs_bin, range_from_bin(static_cast<double>(abs_bin), p.wf.fs) / 1e3,
                    10.0 * std::log10(static_cast<double>(r.rd_power[pk]) / r.rd_noise_floor_w));
        if (c.run_detection) {
            std::printf("2D CFAR (Pfa %.0e): %zu cells over threshold, %zu detection(s)\n", c.cfar.pfa,
                        r.cfar.cells_over_threshold, r.detections.size());
            for (std::size_t i = 0; i < r.detections.size() && i < 8; ++i) {
                const auto& d = r.detections[i];
                std::printf("  #%zu Doppler %.2f Hz  %.3f m/s | range %.3f km | %.2f dB above floor\n", i + 1, d.doppler_hz,
                            d.velocity_mps, d.range_m / 1e3, d.snr_db);
            }
        }
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
