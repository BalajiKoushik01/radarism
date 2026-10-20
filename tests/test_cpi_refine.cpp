#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <vector>
#include "radar/cpi.hpp"
#include "radar/export_hdf5.hpp"

using namespace radar;

static SimulationParams base() {
    SimulationParams p;
    p.wf.type = WaveformType::LFM; p.wf.fs = 1e6; p.wf.pri = 1e-3; p.wf.tau = 100e-6; p.wf.bandwidth = 2e5;
    p.ch.fc = 10e9; p.ch.range_m = 30e3; p.ch.snr_db = 30.0;
    return p;
}

TEST(CpiRefine, RangeRefinementBeatsTheBinGrid) {
    auto p = base();
    p.ch.range_m = kSpeedOfLight * (200.4 / p.wf.fs) / 2.0;
    CpiParams c;
    c.pulses = 32; c.gate_len = 128;
    const auto r = simulate_cpi(p, c);
    ASSERT_FALSE(r->detections.empty());
    const auto& d = r->detections.front();
    const double binm = kSpeedOfLight / (2.0 * p.wf.fs);
    const double coarse = std::fabs(d.range_m - p.ch.range_m), refined = std::fabs(d.range_refined_m - p.ch.range_m);
    EXPECT_LT(refined, 0.2 * binm);
    EXPECT_LT(refined, coarse);
}

TEST(CpiRefine, DopplerRefinementBeatsTheBinGrid) {
    auto p = base();
    p.ch.velocity_mps = 258.0 * kSpeedOfLight / (2.0 * p.ch.fc);   // 16.5 bins at PRF 1 kHz and 64 rows
    CpiParams c;
    c.pulses = 64; c.gate_len = 128;
    const auto r = simulate_cpi(p, c);
    ASSERT_FALSE(r->detections.empty());
    const auto& d = r->detections.front();
    const double lambda = kSpeedOfLight / p.ch.fc;
    EXPECT_LT(std::fabs(d.doppler_refined_hz - 258.0), 2.4);                          // 0.15 bin
    EXPECT_LT(std::fabs(d.doppler_refined_hz - 258.0), std::fabs(d.doppler_hz - 258.0));
    EXPECT_NEAR(d.velocity_refined_mps, d.doppler_refined_hz * lambda / 2.0, 1e-9);
}

TEST(CpiExport, Hdf5HoldsTheMapDetectionsAndMetadata) {
    auto p = base();
    CpiParams c;
    c.pulses = 16; c.gate_len = 64;
    const auto r = simulate_cpi(p, c);
    const auto path = (std::filesystem::temp_directory_path() / "radarsim_cpi_test.h5").string();
    write_hdf5_cpi(path, *r);
    const hid_t f = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(f, 0);
    const hid_t d = H5Dopen2(f, "rd_power", H5P_DEFAULT);
    ASSERT_GE(d, 0);
    const hid_t sp = H5Dget_space(d);
    hsize_t dims[2] = {0, 0};
    ASSERT_EQ(H5Sget_simple_extent_ndims(sp), 2);
    H5Sget_simple_extent_dims(sp, dims, nullptr);
    EXPECT_EQ(dims[0], r->doppler_rows);
    EXPECT_EQ(dims[1], r->gate_len);
    std::vector<float> back(dims[0] * dims[1]);
    ASSERT_GE(H5Dread(d, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, back.data()), 0);
    for (std::size_t i = 0; i < back.size(); ++i) ASSERT_EQ(back[i], r->rd_power[i]) << i;
    const hid_t a = H5Aopen(f, "doppler_rows", H5P_DEFAULT);
    ASSERT_GE(a, 0);
    long long rows = 0;
    H5Aread(a, H5T_NATIVE_LLONG, &rows);
    EXPECT_EQ(static_cast<std::size_t>(rows), r->doppler_rows);
    H5Aclose(a); H5Sclose(sp); H5Dclose(d); H5Fclose(f);
}
