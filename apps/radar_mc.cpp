#include <chrono>
#include <cmath>
#include <cstdio>
#include "args.hpp"
#include "radar/montecarlo.hpp"

using namespace radar;

static void usage() {
    std::puts(
        "radar_mc [key=value ...]   Monte Carlo Pd / Pfa validation against Marcum-Q and Albersheim\n"
        "\n"
        "Waveform : type=pulse|lfm|barker|frank  fs=1e6  pri=10e-3  tau=500e-6  bw=2e5  code=13  mfw=none|hamming|hann\n"
        "Detector : pfa=1e-3  guard=20  train=64  tol=2 (bins)  corr=off|on (correlation-aware CFAR threshold)\n"
        "Sweep    : snr_lo= snr_hi= snr_step=1   (matched-filter OUTPUT SNR in dB; default = Albersheim(0.9) -4 .. +4)\n"
        "Monte    : trials=400 (per point)  ntrials=200 (noise-only runs)  bin=0 (0 = 30% of PRI)  threads=0 (all)  seed=12345\n"
        "\n"
        "Example: ./build/radar_mc pfa=1e-4 trials=1000 mfw=hamming");
}

int main(int argc, char** argv) {
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string s = argv[i];
            if (s == "help" || s == "-h" || s == "--help") { usage(); return 0; }
        }
        Args a(argc, argv);
        McParams mp;
        SimulationParams& p = mp.sim;
        const std::string type = a.str("type", "lfm");
        if (type == "pulse") p.wf.type = WaveformType::Pulse;
        else if (type == "lfm") p.wf.type = WaveformType::LFM;
        else if (type == "barker") p.wf.type = WaveformType::Barker;
        else if (type == "frank") p.wf.type = WaveformType::Frank;
        else throw std::invalid_argument("type must be pulse, lfm, barker or frank");
        p.wf.fs = a.num("fs", 1e6);
        p.wf.pri = a.num("pri", 10e-3);
        p.wf.tau = a.num("tau", 500e-6);
        p.wf.bandwidth = a.num("bw", 2e5);
        p.wf.code_length = static_cast<int>(a.num("code", 13));
        p.ch.fc = a.num("fc", 10e9);
        p.ch.noise_figure_db = a.num("nf", 3.0);
        p.ch.seed = static_cast<std::uint64_t>(a.num("seed", 12345));
        const std::string mfw = a.str("mfw", "none");
        if (mfw == "hamming") p.mf_weighting = MfWeighting::Hamming;
        else if (mfw == "hann") p.mf_weighting = MfWeighting::Hann;
        else if (mfw != "none") throw std::invalid_argument("mfw must be none, hamming or hann");
        p.cfar.guard_cells = static_cast<int>(a.num("guard", 20));
        p.cfar_correlation_aware = (a.str("corr", "off") == "on");
        p.cfar.train_cells = static_cast<int>(a.num("train", 64));
        mp.pfa = a.num("pfa", 1e-3);
        mp.trials = static_cast<int>(a.num("trials", 400));
        mp.noise_trials = static_cast<int>(a.num("ntrials", 200));
        mp.target_bin = static_cast<std::size_t>(a.num("bin", 0));
        mp.threads = static_cast<unsigned>(a.num("threads", 0));
        mp.match_tol_bins = a.num("tol", 2.0);

        const double alb = albersheim_snr_db(0.9, mp.pfa);
        const double lo = a.has("snr_lo") ? a.num("snr_lo", 0.0) : std::floor(alb - 4.0);
        const double hi = a.has("snr_hi") ? a.num("snr_hi", 0.0) : std::ceil(alb + 4.0);
        const double step = a.num("snr_step", 1.0);
        a.check_unused();
        if (!(step > 0.0) || hi < lo) throw std::invalid_argument("bad SNR sweep");
        std::vector<double> snrs;
        for (double v = lo; v <= hi + 1e-9; v += step) snrs.push_back(v);

        const auto t0 = std::chrono::steady_clock::now();
        const McResult r = run_monte_carlo(mp, snrs);
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

        std::printf("== Monte Carlo: %s, N = %zu, pulse = %zu samples, MF gain +%.2f dB, target bin %zu, design Pfa %.1e ==\n",
                    waveform_name(p.wf.type), r.n, r.n_pulse, r.mf_gain_db, r.target_bin, mp.pfa);
        std::printf("%d trials per point x %zu points + %d noise-only runs on %u threads: %.1f s\n\n", mp.trials,
                    snrs.size(), mp.noise_trials, r.threads, sec);

        const McPfa& f = r.pfa;
        std::printf("False alarms (noise only, %zu cells):\n", f.cells);
        std::printf("  fixed threshold (T = %.2f x noise floor): measured %.3e  design %.3e  ratio %.3f   <- validates noise scaling\n",
                    -std::log(mp.pfa), f.pfa_fixed, f.pfa_design, f.pfa_fixed / f.pfa_design);
        std::printf("  CA-CFAR (alpha %.2f, %d training cells):  measured %.3e  ratio %.2f   | %.2f clusters per run | ideal CFAR loss %.2f dB\n\n",
                    f.alpha, 2 * p.cfar.train_cells, f.pfa_cfar, f.pfa_cfar / f.pfa_design, f.clusters_per_run,
                    f.expected_cfar_loss_db);

        std::printf("  CFAR training cells: effective %.1f of %d (%s)\n\n",
                    f.effective_cells > 0.0 ? f.effective_cells : 2.0 * p.cfar.train_cells, 2 * p.cfar.train_cells,
                    p.cfar_correlation_aware ? "correlation-aware threshold" : "independence assumed");
        std::printf("SNR_out(dB)  SNR_in(dB)  Pd theory   Pd fixed-threshold [95%% CI]    Pd CA-CFAR [95%% CI]\n");
        int inside = 0;
        for (const auto& pt : r.curve) {
            const bool ok = pt.pd_theory >= pt.lo_fixed && pt.pd_theory <= pt.hi_fixed;
            inside += ok ? 1 : 0;
            std::printf("  %6.2f     %7.2f     %.3f      %.3f [%.3f, %.3f] %s     %.3f [%.3f, %.3f]\n", pt.snr_out_db,
                        pt.snr_in_db, pt.pd_theory, pt.pd_fixed, pt.lo_fixed, pt.hi_fixed, ok ? " " : "*", pt.pd_cfar,
                        pt.lo_cfar, pt.hi_cfar);
        }
        std::printf("\nRequired output SNR for Pd = 0.9:  theory (Marcum) %.2f dB | Albersheim %.2f dB | measured fixed %.2f dB | measured CA-CFAR %.2f dB\n",
                    r.theory_snr90, r.albersheim_snr90, r.measured_snr90_fixed, r.measured_snr90_cfar);
        std::printf("Agreement: theory inside the 95%% CI at %d of %zu points (expect about 95%%; '*' marks a miss)\n", inside,
                    r.curve.size());
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
