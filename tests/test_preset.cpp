#include <gtest/gtest.h>
#include <filesystem>
#include "radar/preset.hpp"

using namespace radar;

static Preset custom() {
    Preset p;
    p.cpi_mode = true; p.psd_window = 1;
    auto& s = p.sim;
    s.wf.type = WaveformType::Barker; s.wf.fs = 2.5e6; s.wf.pri = 2e-3; s.wf.tau = 65e-6; s.wf.bandwidth = 3e5; s.wf.code_length = 13;
    s.ch.fc = 9.4e9; s.ch.range_m = 41234.5; s.ch.velocity_mps = -12.25; s.ch.snr_mode = SnrMode::Physical;
    s.ch.snr_db = 3.5; s.ch.pt_w = 750; s.ch.gt_db = 31; s.ch.gr_db = 29; s.ch.rcs_m2 = 0.35; s.ch.loss_db = 2.5;
    s.ch.noise_figure_db = 4.5; s.ch.seed = 987654321;
    s.phase_noise_enabled = true; s.pn.floor_dbc_hz = -95.5; s.pn.leeson_corner_hz = 2500; s.pn.flicker_corner_hz = 321; s.pn.seed = 77;
    s.cfar.pfa = 3.3e-7; s.cfar.guard_cells = 11; s.cfar.train_cells = 77;
    s.mf_weighting = MfWeighting::Hamming; s.cfar_correlation_aware = true;
    s.rx_imp.adc_bits = 12; s.rx_imp.adc_fs_sigma = 5.5; s.rx_imp.iq_gain_db = 0.75; s.rx_imp.iq_phase_deg = -2.5;
    s.rx_imp.dc_i_sigma = 0.25; s.rx_imp.dc_q_sigma = -0.5;
    s.extra_targets = {{50e3, 1.5, -6.0}, {61e3, -2.5, -9.5}};
    p.cpi.pulses = 48; p.cpi.gate_len = 300; p.cpi.doppler_fft = 128; p.cpi.window = MfWeighting::Hamming;
    p.cpi.cfar.pfa = 1e-5; p.cpi.cfar.guard_range = 7; p.cpi.cfar.guard_doppler = 3; p.cpi.cfar.train_range = 18; p.cpi.cfar.train_doppler = 5;
    return p;
}

TEST(Preset, RoundTripKeepsEveryField) {
    const Preset a = custom(), b = preset_from_json(preset_to_json(a));
    const auto &s = a.sim, &t = b.sim;
    EXPECT_EQ(b.cpi_mode, a.cpi_mode); EXPECT_EQ(b.psd_window, a.psd_window);
    EXPECT_EQ(t.wf.type, s.wf.type); EXPECT_EQ(t.wf.code_length, s.wf.code_length);
    EXPECT_DOUBLE_EQ(t.wf.fs, s.wf.fs); EXPECT_DOUBLE_EQ(t.wf.pri, s.wf.pri); EXPECT_DOUBLE_EQ(t.wf.tau, s.wf.tau);
    EXPECT_DOUBLE_EQ(t.wf.bandwidth, s.wf.bandwidth);
    EXPECT_DOUBLE_EQ(t.ch.fc, s.ch.fc); EXPECT_DOUBLE_EQ(t.ch.range_m, s.ch.range_m); EXPECT_DOUBLE_EQ(t.ch.velocity_mps, s.ch.velocity_mps);
    EXPECT_EQ(t.ch.snr_mode, s.ch.snr_mode); EXPECT_DOUBLE_EQ(t.ch.snr_db, s.ch.snr_db); EXPECT_DOUBLE_EQ(t.ch.pt_w, s.ch.pt_w);
    EXPECT_DOUBLE_EQ(t.ch.gt_db, s.ch.gt_db); EXPECT_DOUBLE_EQ(t.ch.gr_db, s.ch.gr_db); EXPECT_DOUBLE_EQ(t.ch.rcs_m2, s.ch.rcs_m2);
    EXPECT_DOUBLE_EQ(t.ch.loss_db, s.ch.loss_db); EXPECT_DOUBLE_EQ(t.ch.noise_figure_db, s.ch.noise_figure_db); EXPECT_EQ(t.ch.seed, s.ch.seed);
    EXPECT_EQ(t.phase_noise_enabled, s.phase_noise_enabled); EXPECT_DOUBLE_EQ(t.pn.floor_dbc_hz, s.pn.floor_dbc_hz);
    EXPECT_DOUBLE_EQ(t.pn.leeson_corner_hz, s.pn.leeson_corner_hz); EXPECT_DOUBLE_EQ(t.pn.flicker_corner_hz, s.pn.flicker_corner_hz);
    EXPECT_EQ(t.pn.seed, s.pn.seed);
    EXPECT_DOUBLE_EQ(t.cfar.pfa, s.cfar.pfa); EXPECT_EQ(t.cfar.guard_cells, s.cfar.guard_cells); EXPECT_EQ(t.cfar.train_cells, s.cfar.train_cells);
    EXPECT_EQ(t.mf_weighting, s.mf_weighting); EXPECT_EQ(t.cfar_correlation_aware, s.cfar_correlation_aware);
    EXPECT_EQ(t.rx_imp.adc_bits, s.rx_imp.adc_bits); EXPECT_DOUBLE_EQ(t.rx_imp.adc_fs_sigma, s.rx_imp.adc_fs_sigma);
    EXPECT_DOUBLE_EQ(t.rx_imp.iq_gain_db, s.rx_imp.iq_gain_db); EXPECT_DOUBLE_EQ(t.rx_imp.iq_phase_deg, s.rx_imp.iq_phase_deg);
    EXPECT_DOUBLE_EQ(t.rx_imp.dc_i_sigma, s.rx_imp.dc_i_sigma); EXPECT_DOUBLE_EQ(t.rx_imp.dc_q_sigma, s.rx_imp.dc_q_sigma);
    ASSERT_EQ(t.extra_targets.size(), 2u);
    for (std::size_t i = 0; i < 2; ++i) {
        EXPECT_DOUBLE_EQ(t.extra_targets[i].range_m, s.extra_targets[i].range_m);
        EXPECT_DOUBLE_EQ(t.extra_targets[i].velocity_mps, s.extra_targets[i].velocity_mps);
        EXPECT_DOUBLE_EQ(t.extra_targets[i].rel_db, s.extra_targets[i].rel_db);
    }
    EXPECT_EQ(b.cpi.pulses, a.cpi.pulses); EXPECT_EQ(b.cpi.gate_len, a.cpi.gate_len); EXPECT_EQ(b.cpi.doppler_fft, a.cpi.doppler_fft);
    EXPECT_EQ(b.cpi.window, a.cpi.window); EXPECT_DOUBLE_EQ(b.cpi.cfar.pfa, a.cpi.cfar.pfa);
    EXPECT_EQ(b.cpi.cfar.guard_range, a.cpi.cfar.guard_range); EXPECT_EQ(b.cpi.cfar.guard_doppler, a.cpi.cfar.guard_doppler);
    EXPECT_EQ(b.cpi.cfar.train_range, a.cpi.cfar.train_range); EXPECT_EQ(b.cpi.cfar.train_doppler, a.cpi.cfar.train_doppler);
}

TEST(Preset, MissingKeysKeepDefaultsAndUnknownKeysAreIgnored) {
    const Preset p = preset_from_json("{ \"fs_hz\": 2e6,\n \"future_key\": \"x\" , \"extra_flag\": true }");
    const SimulationParams d;
    EXPECT_DOUBLE_EQ(p.sim.wf.fs, 2e6);
    EXPECT_DOUBLE_EQ(p.sim.wf.pri, d.wf.pri);
    EXPECT_EQ(p.sim.wf.type, d.wf.type);
    EXPECT_FALSE(p.cpi_mode);
    EXPECT_TRUE(p.sim.extra_targets.empty());
}

TEST(Preset, MalformedInputThrows) {
    for (const char* bad : {"", "{", "{\"a\" 1}", "{\"a\": }", "{\"waveform\": \"nope\"}", "{\"fs_hz\": \"fast\"}", "[1,2]"})
        EXPECT_THROW(preset_from_json(bad), std::invalid_argument) << bad;
}

TEST(Preset, FileRoundTrip) {
    const auto path = (std::filesystem::temp_directory_path() / "radarsim_preset_test.json").string();
    save_preset(path, custom());
    const Preset p = load_preset(path);
    EXPECT_DOUBLE_EQ(p.sim.ch.range_m, 41234.5);
    EXPECT_EQ(p.sim.extra_targets.size(), 2u);
    EXPECT_THROW(load_preset(path + ".missing"), std::runtime_error);
}
