#include <gtest/gtest.h>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <type_traits>
#include "radar/export_hdf5.hpp"
#include "radar/exporters.hpp"
#include "radar/simulation.hpp"

using namespace radar;

static SimulationParams base() {
    SimulationParams p;
    p.wf.type = WaveformType::LFM; p.wf.fs = 1e6; p.wf.pri = 10e-3; p.wf.tau = 2e-3; p.wf.bandwidth = 2e5;
    p.ch.snr_db = -17.0;
    p.ch.range_m = kSpeedOfLight * (3000.0 / p.wf.fs) / 2.0;   // delay = 3000 samples
    return p;
}

static std::filesystem::path tmp(const std::string& name) {
    auto d = std::filesystem::temp_directory_path() / "radarsim_tests";
    std::filesystem::create_directories(d);
    return d / name;
}
static std::string slurp(const std::filesystem::path& p) {
    std::ifstream f(p); std::stringstream s; s << f.rdbuf(); return s.str();
}

TEST(Pipeline, DetectsTargetAndReturnsImmutableSharedResult) {
    auto p = base();
    auto res = simulate(p);
    static_assert(std::is_same_v<decltype(res), std::shared_ptr<const SimulationResult>>);
    ASSERT_TRUE(res->has_range_profile);
    EXPECT_EQ(res->rx.size(), 10000u);
    bool found = false;
    for (auto d : res->cfar.detections) if (std::llabs(static_cast<long long>(d) - 3000) <= 2) found = true;
    EXPECT_TRUE(found);
    EXPECT_EQ(res->detection_range_m.size(), res->cfar.detections.size());
    EXPECT_EQ(res->tx, generate(p.wf).samples);   // exported Tx is the clean waveform
}

TEST(Pipeline, PhaseNoiseChangesRxButNotCleanTxAndStillDetects) {
    auto p = base();
    auto clean = simulate(p);
    p.phase_noise_enabled = true; p.pn.floor_dbc_hz = -90.0;
    auto noisy = simulate(p);
    EXPECT_EQ(noisy->tx, clean->tx);
    std::size_t diff = 0;
    for (std::size_t n = 3000; n < 5000; ++n) if (noisy->rx[n] != clean->rx[n]) ++diff;
    EXPECT_GT(diff, 1000u);
    bool found = false;
    for (auto d : noisy->cfar.detections) if (std::llabs(static_cast<long long>(d) - 3000) <= 2) found = true;
    EXPECT_TRUE(found);
}

TEST(Pipeline, CwHasNoRangeProfile) {
    auto p = base(); p.wf.type = WaveformType::CW;
    auto res = simulate(p);
    EXPECT_FALSE(res->has_range_profile);
    EXPECT_TRUE(res->mf_power.empty());
    EXPECT_TRUE(res->cfar.detections.empty());
}

TEST(Pipeline, IsDeterministicForFixedSeeds) {
    auto p = base(); p.phase_noise_enabled = true;
    EXPECT_EQ(simulate(p)->rx, simulate(p)->rx);
}

TEST(Pipeline, CancellationAndValidation) {
    auto p = base();
    std::atomic<bool> stop{true};
    EXPECT_THROW(simulate(p, &stop), Cancelled);
    p.ch.range_m = -5.0;
    EXPECT_THROW(simulate(p), std::invalid_argument);
}

TEST(Export, Float32NormalizedRoundTripAndSidecar) {
    CBuffer x(1000);
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = {static_cast<float>(1e-7 * std::sin(0.1 * i)), static_cast<float>(5e-8 * std::cos(0.07 * i))};
    const auto path = tmp("f32.bin").string();
    auto info = write_bin(path, x, BinFormat::Float32, 1e6, 10e9, "test");
    EXPECT_EQ(std::filesystem::file_size(path), 8u * x.size());
    std::vector<float> v(2 * x.size());
    std::ifstream(path, std::ios::binary).read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(float)));
    float vmax = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        vmax = std::max({vmax, std::fabs(v[2 * i]), std::fabs(v[2 * i + 1])});
        EXPECT_NEAR(v[2 * i] / info.scale, x[i].real(), 1e-6 * info.peak_abs);
        EXPECT_NEAR(v[2 * i + 1] / info.scale, x[i].imag(), 1e-6 * info.peak_abs);
    }
    EXPECT_NEAR(vmax, 0.9f, 1e-4f);
    const std::string js = slurp(info.json_path);
    EXPECT_NE(js.find("\"format\": \"float32\""), std::string::npos);
    EXPECT_NE(js.find("\"sample_rate_hz\""), std::string::npos);
    EXPECT_NE(js.find("\"num_samples\": 1000"), std::string::npos);
}

TEST(Export, Int16ScalingAndRoundTrip) {
    CBuffer x(1000);
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = {static_cast<float>(3e-7 * std::sin(0.1 * i)), static_cast<float>(2e-7 * std::cos(0.07 * i))};
    const auto path = tmp("i16.bin").string();
    auto info = write_bin(path, x, BinFormat::Int16, 1e6, 10e9, "test");
    EXPECT_EQ(std::filesystem::file_size(path), 4u * x.size());
    EXPECT_EQ(info.clipped, 0u);
    std::vector<std::int16_t> v(2 * x.size());
    std::ifstream(path, std::ios::binary).read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(std::int16_t)));
    int vmax = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        vmax = std::max({vmax, std::abs(static_cast<int>(v[2 * i])), std::abs(static_cast<int>(v[2 * i + 1]))});
        EXPECT_NEAR(v[2 * i] / info.scale, x[i].real(), 0.51 / info.scale);
    }
    EXPECT_NEAR(vmax, 0.9 * 32767.0, 2.0);
    EXPECT_NE(slurp(info.json_path).find("\"format\": \"int16\""), std::string::npos);
}

TEST(Export, CsvHardLimit) {
    CBuffer big(kCsvMaxSamples + 1), ok(kCsvMaxSamples);
    const auto pbig = tmp("big.csv"), pok = tmp("ok.csv");
    std::filesystem::remove(pbig);
    EXPECT_THROW(write_csv(pbig.string(), big, 1e6), std::length_error);
    EXPECT_FALSE(std::filesystem::exists(pbig));
    EXPECT_NO_THROW(write_csv(pok.string(), ok, 1e6));
    std::ifstream f(pok); std::string line; std::size_t lines = 0;
    while (std::getline(f, line)) ++lines;
    EXPECT_EQ(lines, kCsvMaxSamples + 1);   // header + rows
}

TEST(Export, Hdf5RoundTripWithMetadata) {
    auto p = base();
    auto res = simulate(p);
    const auto path = tmp("sim.h5").string();
    write_hdf5(path, *res);

    hid_t f = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(f, 0);
    double fs = 0, fc = 0;
    hid_t a = H5Aopen(f, "sample_rate_hz", H5P_DEFAULT); ASSERT_GE(a, 0);
    H5Aread(a, H5T_NATIVE_DOUBLE, &fs); H5Aclose(a);
    a = H5Aopen(f, "carrier_hz", H5P_DEFAULT); ASSERT_GE(a, 0);
    H5Aread(a, H5T_NATIVE_DOUBLE, &fc); H5Aclose(a);
    EXPECT_DOUBLE_EQ(fs, 1e6);
    EXPECT_DOUBLE_EQ(fc, 10e9);

    H5Handle ct = h5_complex_type();
    hid_t d = H5Dopen2(f, "rx", H5P_DEFAULT); ASSERT_GE(d, 0);
    CBuffer back(res->rx.size());
    ASSERT_GE(H5Dread(d, ct.get(), H5S_ALL, H5S_ALL, H5P_DEFAULT, back.data()), 0);
    H5Dclose(d);
    EXPECT_EQ(back, res->rx);
    hid_t m = H5Dopen2(f, "mf_power", H5P_DEFAULT);
    EXPECT_GE(m, 0);
    if (m >= 0) H5Dclose(m);
    H5Fclose(f);
}

TEST(Pipeline, KeepsImpairedTxOnlyWhenPhaseNoiseEnabled) {
    auto p = base();
    EXPECT_TRUE(simulate(p)->tx_impaired.empty());
    p.phase_noise_enabled = true; p.pn.floor_dbc_hz = -60.0;
    auto r = simulate(p);
    ASSERT_EQ(r->tx_impaired.size(), r->tx.size());
    EXPECT_NE(r->tx_impaired, r->tx);
    for (std::size_t n = 0; n < r->n_pulse; ++n) ASSERT_NEAR(std::abs(r->tx_impaired[n]), 1.0f, 1e-4f);
}
