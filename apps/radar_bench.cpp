#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <thread>
#include <vector>
#include "args.hpp"
#include "radar/fft.hpp"
#include "radar/simulation.hpp"

using namespace radar;
using Clock = std::chrono::steady_clock;

static double peak_rss_mb() {
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line))
        if (line.rfind("VmHWM:", 0) == 0) return std::atof(line.c_str() + 6) / 1024.0;
    return 0.0;
}

int main(int argc, char** argv) {
    try {
        Args a(argc, argv);
        const int threads = static_cast<int>(a.num("threads", 0));   // 0 = leave RADAR_FFT_THREADS alone
        const double maxn = a.num("maxn", 10e6);
        a.check_unused();
        if (threads > 0) setenv("RADAR_FFT_THREADS", std::to_string(threads).c_str(), 1);
        const char* env = std::getenv("RADAR_FFT_THREADS");
#ifdef RADAR_FFTW_THREADS
        const char* lib = "compiled in";
#else
        const char* lib = "not compiled in";
#endif
        std::printf("hardware threads: %u | FFTW threads library: %s | RADAR_FFT_THREADS=%s\n",
                    std::thread::hardware_concurrency(), lib, env ? env : "(unset = 1)");

        std::printf("\n%-12s %-12s %-10s\n", "FFT size", "time (ms)", "GFLOP/s");
        for (std::size_t n : {std::size_t(1) << 16, std::size_t(1) << 18, std::size_t(1) << 20, std::size_t(1) << 22,
                              std::size_t(5000000), std::size_t(10000000)}) {
            if (static_cast<double>(n) > maxn) continue;
            CBuffer x(n), X(n);
            for (std::size_t i = 0; i < n; ++i) x[i] = {static_cast<float>(std::sin(0.001 * static_cast<double>(i))), 0.0f};
            Fft f(n, true);
            f.execute(x.data(), X.data());   // warm-up
            const auto t0 = Clock::now();
            for (int r = 0; r < 5; ++r) f.execute(x.data(), X.data());
            const double ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count() / 5.0;
            std::printf("%-12zu %-12.2f %-10.2f\n", n, ms, 5.0 * static_cast<double>(n) * std::log2(static_cast<double>(n)) / (ms * 1e6));
        }

        std::printf("\n%-12s %-16s %-14s\n", "N samples", "simulate (ms)", "peak RSS (MB)");
        for (std::size_t N : {std::size_t(1000000), std::size_t(5000000), std::size_t(10000000)}) {
            if (static_cast<double>(N) > maxn) continue;
            SimulationParams p;
            p.wf.type = WaveformType::LFM; p.wf.fs = 50e6; p.wf.bandwidth = 10e6; p.wf.tau = 100e-6;
            p.wf.pri = static_cast<double>(N) / p.wf.fs;
            p.ch.range_m = kSpeedOfLight * 0.3 * p.wf.pri / 2.0;
            p.ch.snr_db = -15.0;
            p.cfar_correlation_aware = true;
            const auto t0 = Clock::now();
            const ResultPtr r = simulate(p);
            const double ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
            std::printf("%-12zu %-16.0f %-14.0f (%zu detection(s))\n", N, ms, peak_rss_mb(), r->cfar.detections.size());
        }
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
