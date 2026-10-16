#pragma once
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include "radar/types.hpp"

namespace radar {

static_assert(std::endian::native == std::endian::little, "Exporters assume a little-endian host.");
static_assert(sizeof(cf32) == 8, "cf32 must be two packed floats.");

inline constexpr std::size_t kCsvMaxSamples = 100'000;

struct FileCloser {
    void operator()(std::FILE* f) const { if (f) std::fclose(f); }
};
using FilePtr = std::unique_ptr<std::FILE, FileCloser>;

inline FilePtr open_file(const std::string& path, const char* mode) {
    FilePtr f(std::fopen(path.c_str(), mode));
    if (!f) throw std::runtime_error("cannot open file: " + path);
    return f;
}

// Hard safety limit: refuses before creating any file.
inline void write_csv(const std::string& path, const CBuffer& x, double fs) {
    if (x.size() > kCsvMaxSamples)
        throw std::length_error("CSV export disabled: " + std::to_string(x.size()) +
                                " samples exceeds the " + std::to_string(kCsvMaxSamples) + "-sample limit.");
    auto f = open_file(path, "w");
    std::fprintf(f.get(), "index,time_s,i,q\n");
    for (std::size_t i = 0; i < x.size(); ++i)
        std::fprintf(f.get(), "%zu,%.9e,%.9e,%.9e\n", i, static_cast<double>(i) / fs,
                     static_cast<double>(x[i].real()), static_cast<double>(x[i].imag()));
    if (std::fflush(f.get()) != 0) throw std::runtime_error("write failed: " + path);
}

enum class BinFormat { Float32, Int16 };

struct BinExportInfo {
    std::string bin_path;
    std::string json_path;
    double scale = 1.0;          // stored = original * scale   (original = stored / scale)
    double peak_abs = 0.0;       // original peak of |I| and |Q|
    std::size_t clipped = 0;
};

// Interleaved I,Q little-endian + "<path>.json" sidecar.
// normalize=true (default): peak |I|,|Q| maps to full_scale_fraction (1.0 for float32, 32767 for int16).
// normalize=false: only valid for Float32; raw amplitudes (sqrt-watts) are written with scale = 1.
inline BinExportInfo write_bin(const std::string& path, const CBuffer& x, BinFormat fmt, double fs, double fc,
                               const std::string& role, bool normalize = true, double full_scale_fraction = 0.9) {
    if (x.empty()) throw std::invalid_argument("write_bin: empty buffer");
    if (!(full_scale_fraction > 0.0 && full_scale_fraction <= 1.0))
        throw std::invalid_argument("write_bin: full_scale_fraction must be in (0, 1]");

    double peak = 0.0;
    for (const auto& s : x)
        peak = std::max({peak, static_cast<double>(std::fabs(s.real())), static_cast<double>(std::fabs(s.imag()))});

    BinExportInfo info;
    info.bin_path = path;
    info.json_path = path + ".json";
    info.peak_abs = peak;

    const bool raw_float = (fmt == BinFormat::Float32 && !normalize);
    const double full = (fmt == BinFormat::Int16) ? 32767.0 : 1.0;
    info.scale = raw_float ? 1.0 : (peak > 0.0 ? full_scale_fraction * full / peak : 1.0);

    constexpr std::size_t kChunk = 1u << 16;
    const std::size_t n = x.size();
    {
        auto f = open_file(path, "wb");
        if (raw_float) {
            if (std::fwrite(x.data(), sizeof(cf32), n, f.get()) != n) throw std::runtime_error("write failed: " + path);
        } else if (fmt == BinFormat::Float32) {
            std::vector<float> buf(2 * kChunk);
            for (std::size_t off = 0; off < n; off += kChunk) {
                const std::size_t m = std::min(kChunk, n - off);
                for (std::size_t j = 0; j < m; ++j) {
                    buf[2 * j] = static_cast<float>(x[off + j].real() * info.scale);
                    buf[2 * j + 1] = static_cast<float>(x[off + j].imag() * info.scale);
                }
                if (std::fwrite(buf.data(), sizeof(float), 2 * m, f.get()) != 2 * m) throw std::runtime_error("write failed: " + path);
            }
        } else {
            std::vector<std::int16_t> buf(2 * kChunk);
            auto q = [&](double v) {
                long long r = std::llround(v * info.scale);
                if (r > 32767) { r = 32767; ++info.clipped; }
                if (r < -32767) { r = -32767; ++info.clipped; }
                return static_cast<std::int16_t>(r);
            };
            for (std::size_t off = 0; off < n; off += kChunk) {
                const std::size_t m = std::min(kChunk, n - off);
                for (std::size_t j = 0; j < m; ++j) {
                    buf[2 * j] = q(x[off + j].real());
                    buf[2 * j + 1] = q(x[off + j].imag());
                }
                if (std::fwrite(buf.data(), sizeof(std::int16_t), 2 * m, f.get()) != 2 * m) throw std::runtime_error("write failed: " + path);
            }
        }
        if (std::fflush(f.get()) != 0) throw std::runtime_error("write failed: " + path);
    }

    auto j = open_file(info.json_path, "w");
    const std::string fname = std::filesystem::path(path).filename().string();
    std::fprintf(j.get(),
        "{\n"
        "  \"file\": \"%s\",\n"
        "  \"role\": \"%s\",\n"
        "  \"format\": \"%s\",\n"
        "  \"layout\": \"interleaved I,Q\",\n"
        "  \"endianness\": \"little\",\n"
        "  \"sample_rate_hz\": %.17g,\n"
        "  \"carrier_hz\": %.17g,\n"
        "  \"num_samples\": %zu,\n"
        "  \"bit_depth\": %d,\n"
        "  \"scale_counts_per_unit\": %.17g,\n"
        "  \"peak_abs_original\": %.17g,\n"
        "  \"full_scale_fraction\": %.17g,\n"
        "  \"clipped_samples\": %zu,\n"
        "  \"units\": \"original amplitude in sqrt(watts) (|x|^2 = watts); original = stored / scale_counts_per_unit\"\n"
        "}\n",
        fname.c_str(), role.c_str(), fmt == BinFormat::Int16 ? "int16" : "float32",
        fs, fc, n, fmt == BinFormat::Int16 ? 16 : 32, info.scale, peak,
        raw_float ? 0.0 : full_scale_fraction, info.clipped);
    if (std::fflush(j.get()) != 0) throw std::runtime_error("write failed: " + info.json_path);
    return info;
}

}  // namespace radar
