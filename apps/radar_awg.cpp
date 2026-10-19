#include <cmath>
#include <cstdio>
#include <string>
#include "args.hpp"
#include "radar/awg_profile.hpp"
#include "radar/simulation.hpp"

using namespace radar;

static void usage() {
    std::puts(
        "radar_awg [key=value ...]   Build an AWG / SDR playback file from the simulated Tx (or Rx)\n"
        "\n"
        "Waveform : type=cw|pulse|lfm|barker|frank  fs=1e6 (baseband)  pri=1e-3  tau=100e-6  bw=2e5  code=13\n"
        "Source   : src=tx|rx   (rx also uses fc=10e9 range=30000 snr=-10 seed=12345)\n"
        "AWG      : awg_fs=1e8 (integer multiple of fs)  mode=iq|if  if=0 (Hz, mode=if)  snap=on|off\n"
        "           gran=1  minlen=0  bits=16  justify=right|left  head=0.9  fmt=int16|float32\n"
        "Output   : out=PREFIX  (writes PREFIX.bin and PREFIX.bin.json; omit for a dry run)\n"
        "\n"
        "Example: ./build/radar_awg awg_fs=1e8 mode=if if=10e6 gran=32 bits=14 justify=left out=/tmp/tx_awg");
}

int main(int argc, char** argv) {
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string s = argv[i];
            if (s == "help" || s == "-h" || s == "--help") { usage(); return 0; }
        }
        Args a(argc, argv);
        SimulationParams p;
        p.run_detection = false;
        const std::string type = a.str("type", "lfm");
        if (type == "cw") p.wf.type = WaveformType::CW;
        else if (type == "pulse") p.wf.type = WaveformType::Pulse;
        else if (type == "lfm") p.wf.type = WaveformType::LFM;
        else if (type == "barker") p.wf.type = WaveformType::Barker;
        else if (type == "frank") p.wf.type = WaveformType::Frank;
        else throw std::invalid_argument("type must be cw, pulse, lfm, barker or frank");
        p.wf.fs = a.num("fs", 1e6);
        p.wf.pri = a.num("pri", 1e-3);
        p.wf.tau = a.num("tau", 100e-6);
        p.wf.bandwidth = a.num("bw", 2e5);
        p.wf.code_length = static_cast<int>(a.num("code", 13));
        p.ch.fc = a.num("fc", 10e9);
        p.ch.range_m = a.num("range", 30e3);
        p.ch.snr_db = a.num("snr", -10.0);
        p.ch.seed = static_cast<std::uint64_t>(a.num("seed", 12345));
        const std::string src = a.str("src", "tx");
        if (src != "tx" && src != "rx") throw std::invalid_argument("src must be tx or rx");

        AwgParams ap;
        ap.fs_awg = a.num("awg_fs", 1e8);
        const std::string mode = a.str("mode", "iq");
        if (mode == "iq") ap.mode = AwgMode::IQ;
        else if (mode == "if") ap.mode = AwgMode::RealIF;
        else throw std::invalid_argument("mode must be iq or if");
        ap.if_hz = a.num("if", 0.0);
        const std::string snap = a.str("snap", "on");
        ap.snap_if = (snap != "off");
        ap.granularity = static_cast<std::size_t>(a.num("gran", 1));
        ap.min_samples = static_cast<std::size_t>(a.num("minlen", 0));
        ap.dac_bits = static_cast<int>(a.num("bits", 16));
        const std::string just = a.str("justify", "right");
        if (just != "right" && just != "left") throw std::invalid_argument("justify must be right or left");
        ap.left_justify = (just == "left");
        ap.headroom = a.num("head", 0.9);
        const std::string fmt = a.str("fmt", "int16");
        if (fmt != "int16" && fmt != "float32") throw std::invalid_argument("fmt must be int16 or float32");
        ap.int16_out = (fmt == "int16");
        const std::string out = a.str("out", "");
        a.check_unused();

        const ResultPtr res = simulate(p);
        const CBuffer& source = (src == "rx") ? res->rx : res->tx;
        const AwgWaveform w = make_awg(source, p.wf.fs, ap);

        std::printf("== AWG profile (%s buffer) ==\n", src == "rx" ? "Rx" : "clean Tx");
        std::printf("Baseband  : %zu samples at %.6g Hz  (PRI %.6g s)\n", w.source_samples, w.fs_baseband, w.source_period_s);
        std::printf("Output    : %zu samples/channel x %zu channel(s) at %.6g Hz  (x%zu interpolation, %s)\n", w.samples, w.channels,
                    w.fs_awg, w.upsample, w.mode == AwgMode::IQ ? "I/Q" : "real digital IF");
        if (w.mode == AwgMode::RealIF) std::printf("IF        : %.9g Hz (whole cycles per PRI: loop is phase-continuous)\n", w.if_used_hz);
        std::printf("Loop      : period %.9g s = PRF %.6f Hz", w.loop_period_s, 1.0 / w.loop_period_s);
        if (w.padded > 0) std::printf("   (zero-padded by %zu samples; source PRF %.6f Hz)", w.padded, 1.0 / w.source_period_s);
        std::printf("\n");
        std::printf("Amplitude : peak %.4g sqrt(W) -> %.0f%% of full scale; %s%d-bit%s, %.6g counts per sqrt(W)\n", w.peak,
                    100.0 * w.headroom, w.int16 ? "" : "float32 (", w.bits, w.int16 ? (w.left_justified ? " left-justified" : " right-justified") : ")",
                    w.counts_per_unit);
        if (w.int16) std::printf("Ideal full-scale-sine SQNR at this resolution: %.1f dB\n", 6.02 * w.bits + 1.76);
        if (!out.empty()) {
            write_awg(out + ".bin", w, src == "rx" ? "rx_for_injection" : "tx_clean_for_awg");
            std::printf("wrote %s.bin and %s.bin.json\n", out.c_str(), out.c_str());
        } else {
            std::puts("(dry run: no out= given, nothing written)");
        }
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
