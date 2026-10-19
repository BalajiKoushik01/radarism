#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>
#include "radar/exporters.hpp"
#include "radar/fft.hpp"
#include "radar/types.hpp"

namespace radar {

// AWG / SDR playback profile for a complex-baseband buffer that is one PRI long (so it loops seamlessly):
//   1. band-limited (FFT) interpolation to the instrument sample rate (integer multiple of the baseband rate)
//   2. optional digital IF -> a single real channel  s[m] = I cos(2 pi f m / Fs) - Q sin(2 pi f m / Fs);
//      f is snapped to a whole number of cycles per PRI so the loop point is phase-continuous
//   3. zero padding to the instrument's minimum length and granularity (changes the loop period; reported)
//   4. scaling to a fraction of DAC full scale, quantised to N bits (optionally left-justified in int16)
enum class AwgMode { IQ, RealIF };

struct AwgParams {
    double fs_awg = 1e8;
    AwgMode mode = AwgMode::IQ;
    double if_hz = 0.0;            // RealIF only
    bool snap_if = true;
    std::size_t granularity = 1;   // output length is a multiple of this
    std::size_t min_samples = 0;
    int dac_bits = 16;             // 8..16
    bool left_justify = false;     // int16 only: put the DAC code in the top bits
    double headroom = 0.9;         // peak as a fraction of full scale
    bool int16_out = true;         // false: float32 normalised to +/- headroom
};

struct AwgWaveform {
    std::vector<std::int16_t> i16;
    std::vector<float> f32;
    std::size_t channels = 1, samples = 0, padded = 0, source_samples = 0, upsample = 1;
    double fs_awg = 0.0, fs_baseband = 0.0, if_used_hz = 0.0, loop_period_s = 0.0, source_period_s = 0.0;
    double scale = 1.0, peak = 0.0, counts_per_unit = 1.0, headroom = 0.9;
    bool int16 = true, left_justified = false;
    int bits = 16;
    AwgMode mode = AwgMode::IQ;
};

// Exact for band-limited periodic signals (zero-padding in frequency; the Nyquist bin is split for even N).
inline CBuffer resample_fft(const CBuffer& x, std::size_t up) {
    const std::size_t N = x.size();
    if (N < 2 || up < 1) throw std::invalid_argument("resample_fft: bad arguments");
    if (up == 1) return x;
    const std::size_t M = N * up;
    CBuffer X(N), Y(M, cf32{0.0f, 0.0f}), y(M);
    Fft fwd(N, true), inv(M, false);
    fwd.execute(x.data(), X.data());
    const bool even = (N % 2 == 0);
    for (std::size_t k = 0; k < N; ++k) {
        if (even && k == N / 2) {
            Y[k] += X[k] * 0.5f;
            Y[M - N / 2] += X[k] * 0.5f;
        } else if (k < (N + 1) / 2) {
            Y[k] = X[k];
        } else {
            Y[M - (N - k)] = X[k];
        }
    }
    inv.execute(Y.data(), y.data());
    const float s = 1.0f / static_cast<float>(N);   // forward length N, inverse length M: y[m*up] = x[m]
    for (auto& v : y) v *= s;
    return y;
}

inline AwgWaveform make_awg(const CBuffer& x, double fs_bb, const AwgParams& p) {
    if (x.size() < 2 || !(fs_bb > 0.0)) throw std::invalid_argument("make_awg: empty input or bad baseband rate");
    if (!(p.fs_awg >= fs_bb)) throw std::invalid_argument("AWG sample rate must be at least the baseband rate.");
    const double ratio = p.fs_awg / fs_bb;
    const std::size_t up = static_cast<std::size_t>(std::llround(ratio));
    if (up < 1 || std::fabs(ratio - static_cast<double>(up)) > 1e-6 * ratio)
        throw std::invalid_argument("AWG sample rate must be an integer multiple of the baseband rate (ratio " +
                                    std::to_string(ratio) + ").");
    if (p.dac_bits < 8 || p.dac_bits > 16) throw std::invalid_argument("DAC resolution must be 8..16 bits.");
    if (!(p.headroom > 0.0 && p.headroom <= 1.0)) throw std::invalid_argument("Headroom must be in (0, 1].");
    const std::size_t N = x.size();
    if (N * up > 200'000'000ULL)
        throw std::invalid_argument("Output would exceed 200M samples; lower the AWG rate or shorten the PRI.");
    const double period = static_cast<double>(N) / fs_bb;

    double fif = 0.0;
    if (p.mode == AwgMode::RealIF) {
        fif = p.if_hz;
        if (p.snap_if) {
            const double k = std::round(fif * period);
            if (k < 1.0) throw std::invalid_argument("IF is below one cycle per PRI.");
            fif = k / period;
        }
        if (fif < 0.5 * fs_bb)
            throw std::invalid_argument("IF must be at least half the baseband rate so the signal does not overlap its image.");
        if (fif + 0.5 * fs_bb > 0.45 * p.fs_awg)
            throw std::invalid_argument("IF plus half the baseband rate exceeds 0.45 x the AWG sample rate.");
    }

    CBuffer y = resample_fft(x, up);
    const std::size_t L0 = y.size();
    const std::size_t g = std::max<std::size_t>(1, p.granularity);
    std::size_t L = std::max(L0, p.min_samples);
    L = (L + g - 1) / g * g;
    const std::size_t ch = (p.mode == AwgMode::IQ) ? 2 : 1;
    std::vector<float> raw(ch * L, 0.0f);
    if (ch == 2) {
        for (std::size_t m = 0; m < L0; ++m) {
            raw[2 * m] = y[m].real();
            raw[2 * m + 1] = y[m].imag();
        }
    } else {
        const double two_pi = 2.0 * std::numbers::pi;
        for (std::size_t m = 0; m < L0; ++m) {
            const double cyc = fif * static_cast<double>(m) / p.fs_awg;
            const double ph = two_pi * (cyc - std::floor(cyc));
            raw[m] = static_cast<float>(static_cast<double>(y[m].real()) * std::cos(ph) -
                                        static_cast<double>(y[m].imag()) * std::sin(ph));
        }
    }
    y = CBuffer();
    float peak = 0.0f;
    for (float v : raw) peak = std::max(peak, std::fabs(v));

    AwgWaveform w;
    w.channels = ch;
    w.samples = L;
    w.padded = L - L0;
    w.source_samples = N;
    w.upsample = up;
    w.fs_awg = p.fs_awg;
    w.fs_baseband = fs_bb;
    w.if_used_hz = fif;
    w.loop_period_s = static_cast<double>(L) / p.fs_awg;
    w.source_period_s = period;
    w.peak = static_cast<double>(peak);
    w.scale = (peak > 0.0f) ? p.headroom / static_cast<double>(peak) : 1.0;
    w.headroom = p.headroom;
    w.int16 = p.int16_out;
    w.bits = p.dac_bits;
    w.left_justified = p.left_justify;
    w.mode = p.mode;
    if (p.int16_out) {
        const long full = (1L << (p.dac_bits - 1)) - 1;
        const long mult = p.left_justify ? (1L << (16 - p.dac_bits)) : 1L;
        w.counts_per_unit = w.scale * static_cast<double>(full) * static_cast<double>(mult);
        w.i16.resize(raw.size());
        for (std::size_t i = 0; i < raw.size(); ++i) {
            long q = std::lround(static_cast<double>(raw[i]) * w.scale * static_cast<double>(full));
            q = std::clamp(q, -full, full);
            w.i16[i] = static_cast<std::int16_t>(q * mult);
        }
    } else {
        w.counts_per_unit = w.scale;
        w.f32.resize(raw.size());
        for (std::size_t i = 0; i < raw.size(); ++i) w.f32[i] = static_cast<float>(static_cast<double>(raw[i]) * w.scale);
    }
    return w;
}

// <path> = raw little-endian samples (interleaved I,Q for IQ mode), <path>.json = everything needed to play it.
inline void write_awg(const std::string& path, const AwgWaveform& w, const std::string& role) {
    {
        auto f = open_file(path, "wb");
        const std::size_t n = w.int16 ? w.i16.size() : w.f32.size();
        const std::size_t esz = w.int16 ? sizeof(std::int16_t) : sizeof(float);
        const void* data = w.int16 ? static_cast<const void*>(w.i16.data()) : static_cast<const void*>(w.f32.data());
        if (std::fwrite(data, esz, n, f.get()) != n) throw std::runtime_error("write failed: " + path);
        if (std::fflush(f.get()) != 0) throw std::runtime_error("write failed: " + path);
    }
    auto j = open_file(path + ".json", "w");
    const std::string fname = std::filesystem::path(path).filename().string();
    std::fprintf(j.get(),
        "{\n"
        "  \"file\": \"%s\",\n"
        "  \"role\": \"%s\",\n"
        "  \"mode\": \"%s\",\n"
        "  \"layout\": \"%s\",\n"
        "  \"format\": \"%s\",\n"
        "  \"endianness\": \"little\",\n"
        "  \"sample_rate_hz\": %.17g,\n"
        "  \"baseband_rate_hz\": %.17g,\n"
        "  \"upsample_factor\": %zu,\n"
        "  \"intermediate_frequency_hz\": %.17g,\n"
        "  \"channels\": %zu,\n"
        "  \"num_samples_per_channel\": %zu,\n"
        "  \"source_samples_per_pri\": %zu,\n"
        "  \"padded_samples\": %zu,\n"
        "  \"loop_period_s\": %.17g,\n"
        "  \"source_period_s\": %.17g,\n"
        "  \"dac_bits\": %d,\n"
        "  \"justification\": \"%s\",\n"
        "  \"peak_fraction_of_full_scale\": %.17g,\n"
        "  \"scale_counts_per_unit\": %.17g,\n"
        "  \"peak_abs_original\": %.17g,\n"
        "  \"units\": \"original amplitude in sqrt(watts); original = stored / scale_counts_per_unit\"\n"
        "}\n",
        fname.c_str(), role.c_str(), w.mode == AwgMode::IQ ? "iq" : "real_if",
        w.mode == AwgMode::IQ ? "interleaved I,Q" : "single real channel", w.int16 ? "int16" : "float32",
        w.fs_awg, w.fs_baseband, w.upsample, w.if_used_hz, w.channels, w.samples, w.source_samples, w.padded,
        w.loop_period_s, w.source_period_s, w.bits, w.int16 ? (w.left_justified ? "left" : "right") : "n/a",
        w.headroom, w.counts_per_unit, w.peak);
    if (std::fflush(j.get()) != 0) throw std::runtime_error("write failed: " + path + ".json");
}

}  // namespace radar
