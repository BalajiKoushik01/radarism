#pragma once
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>
#include "radar/types.hpp"

namespace radar {

enum class WaveformType { CW, Pulse, LFM, Barker, Frank };

struct WaveformParams {
    WaveformType type = WaveformType::Pulse;
    double fs = 10e6;        // baseband sampling rate [Hz]
    double pri = 1e-3;       // pulse repetition interval [s]
    double tau = 100e-6;     // pulse width [s] (ignored for CW)
    double bandwidth = 0.0;  // LFM sweep bandwidth [Hz]
    int code_length = 13;    // Barker length or Frank L = M*M
};

struct Waveform {
    CBuffer samples;                  // length N = round(pri*fs), zero-padded
    std::size_t n_pulse = 0;          // active samples (0 .. n_pulse-1)
    std::size_t samples_per_chip = 0; // coded waveforms only
};

inline constexpr std::size_t kMaxSamples = 100'000'000;

inline std::size_t total_samples(const WaveformParams& p) {
    return static_cast<std::size_t>(std::llround(p.pri * p.fs));
}

inline const std::vector<int>& barker_code(int len) {
    static const std::vector<int> b2{1, -1};
    static const std::vector<int> b3{1, 1, -1};
    static const std::vector<int> b4{1, 1, -1, 1};
    static const std::vector<int> b5{1, 1, 1, -1, 1};
    static const std::vector<int> b7{1, 1, 1, -1, -1, 1, -1};
    static const std::vector<int> b11{1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1};
    static const std::vector<int> b13{1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1};
    static const std::vector<int> none{};
    switch (len) {
        case 2: return b2;   case 3: return b3;   case 4: return b4;
        case 5: return b5;   case 7: return b7;   case 11: return b11;
        case 13: return b13; default: return none;
    }
}

// Returns "" if parameters are valid, otherwise a human-readable error (shown in the UI later).
inline std::string validate(const WaveformParams& p) {
    if (!(p.fs > 0.0)) return "Sampling rate must be > 0.";
    if (!(p.pri > 0.0)) return "PRI must be > 0.";
    const double n_total = std::round(p.pri * p.fs);
    if (n_total < 2.0) return "PRI * Fs gives fewer than 2 samples.";
    if (n_total > static_cast<double>(kMaxSamples)) return "PRI * Fs exceeds the 100M-sample limit.";
    if (p.type == WaveformType::CW) return "";

    if (!(p.tau > 0.0)) return "Pulse width must be > 0.";
    if (p.tau > p.pri) return "Pulse width cannot exceed PRI.";
    if (std::round(p.tau * p.fs) < 2.0) return "Pulse is shorter than 2 samples; raise Fs or tau.";

    switch (p.type) {
        case WaveformType::LFM:
            if (!(p.bandwidth > 0.0)) return "LFM needs sweep bandwidth > 0.";
            if (p.fs < 1.2 * p.bandwidth) return "Fs must be at least 1.2 x sweep bandwidth.";
            break;
        case WaveformType::Barker:
        case WaveformType::Frank: {
            const int L = p.code_length;
            if (p.type == WaveformType::Barker && barker_code(L).empty())
                return "Barker length must be 2, 3, 4, 5, 7, 11 or 13.";
            if (p.type == WaveformType::Frank) {
                const int M = static_cast<int>(std::lround(std::sqrt(static_cast<double>(L))));
                if (L < 4 || M * M != L) return "Frank code length must be a perfect square (4, 9, 16, ...).";
            }
            const double spc = p.fs * p.tau / L;
            if (spc < 1.0 || std::fabs(spc - std::round(spc)) > 1e-6)
                return "Chip width (tau/L) must be an integer number of samples at this Fs.";
            break;
        }
        default: break;
    }
    return "";
}

inline Waveform generate(const WaveformParams& p) {
    if (auto err = validate(p); !err.empty()) throw std::invalid_argument(err);

    const std::size_t N = total_samples(p);
    Waveform w;
    w.samples.assign(N, cf32{0.0f, 0.0f});  // zero-padded listening window by default
    const double two_pi = 2.0 * std::numbers::pi;

    switch (p.type) {
        case WaveformType::CW:
            w.n_pulse = N;
            for (auto& s : w.samples) s = {1.0f, 0.0f};
            break;
        case WaveformType::Pulse:
            w.n_pulse = static_cast<std::size_t>(std::llround(p.tau * p.fs));
            for (std::size_t n = 0; n < w.n_pulse; ++n) w.samples[n] = {1.0f, 0.0f};
            break;
        case WaveformType::LFM: {
            w.n_pulse = static_cast<std::size_t>(std::llround(p.tau * p.fs));
            const double k = p.bandwidth / p.tau;  // chirp rate [Hz/s]
            for (std::size_t n = 0; n < w.n_pulse; ++n) {
                const double t = static_cast<double>(n) / p.fs - 0.5 * p.tau;  // sweep -B/2 .. +B/2
                const double ph = std::numbers::pi * k * t * t;
                w.samples[n] = {static_cast<float>(std::cos(ph)), static_cast<float>(std::sin(ph))};
            }
            break;
        }
        case WaveformType::Barker: {
            const auto& code = barker_code(p.code_length);
            w.samples_per_chip = static_cast<std::size_t>(std::llround(p.fs * p.tau / p.code_length));
            w.n_pulse = w.samples_per_chip * code.size();
            for (std::size_t c = 0; c < code.size(); ++c)
                for (std::size_t s = 0; s < w.samples_per_chip; ++s)
                    w.samples[c * w.samples_per_chip + s] = {static_cast<float>(code[c]), 0.0f};
            break;
        }
        case WaveformType::Frank: {
            const int L = p.code_length;
            const int M = static_cast<int>(std::lround(std::sqrt(static_cast<double>(L))));
            w.samples_per_chip = static_cast<std::size_t>(std::llround(p.fs * p.tau / L));
            w.n_pulse = w.samples_per_chip * static_cast<std::size_t>(L);
            for (int idx = 0; idx < L; ++idx) {
                const int i = idx / M, j = idx % M;
                const double ph = two_pi * static_cast<double>(i * j) / M;
                const cf32 v{static_cast<float>(std::cos(ph)), static_cast<float>(std::sin(ph))};
                for (std::size_t s = 0; s < w.samples_per_chip; ++s)
                    w.samples[static_cast<std::size_t>(idx) * w.samples_per_chip + s] = v;
            }
            break;
        }
    }
    return w;
}

}  // namespace radar
