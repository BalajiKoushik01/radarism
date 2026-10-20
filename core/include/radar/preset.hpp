#pragma once
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "radar/cpi.hpp"
#include "radar/simulation.hpp"

namespace radar {

// Everything the setup page holds, in a flat, human-readable JSON file. Missing keys keep their defaults and
// unknown keys are ignored, so presets survive version changes.
struct Preset {
    SimulationParams sim;
    CpiParams cpi;
    bool cpi_mode = false;
    int psd_window = 0;   // GUI: 0 = Hann, 1 = rectangular
};

namespace detail {
inline std::string jnum(double v) { char b[48]; std::snprintf(b, sizeof b, "%.17g", v); return b; }
inline std::string jstr(const std::string& s) {
    std::string o = "\"";
    for (char c : s) { if (c == '"' || c == '\\') o += '\\'; o += c; }
    return o + "\"";
}
inline const char* wf_name(WaveformType t) {
    switch (t) {
        case WaveformType::CW: return "cw";
        case WaveformType::Pulse: return "pulse";
        case WaveformType::LFM: return "lfm";
        case WaveformType::Barker: return "barker";
        case WaveformType::Frank: return "frank";
    }
    return "pulse";
}
inline WaveformType wf_from(const std::string& s) {
    if (s == "cw") return WaveformType::CW;
    if (s == "pulse") return WaveformType::Pulse;
    if (s == "lfm") return WaveformType::LFM;
    if (s == "barker") return WaveformType::Barker;
    if (s == "frank") return WaveformType::Frank;
    throw std::invalid_argument("preset: unknown waveform '" + s + "'");
}
inline const char* wt_name(MfWeighting w) {
    return w == MfWeighting::Hamming ? "hamming" : w == MfWeighting::Hann ? "hann" : "none";
}
inline MfWeighting wt_from(const std::string& s) {
    if (s == "none") return MfWeighting::None;
    if (s == "hamming") return MfWeighting::Hamming;
    if (s == "hann") return MfWeighting::Hann;
    throw std::invalid_argument("preset: unknown weighting '" + s + "'");
}

// Flat object of "key": number | "string" | true | false. Values are stored as "s:<text>" or "v:<token>".
inline std::map<std::string, std::string> parse_flat_json(const std::string& t) {
    std::map<std::string, std::string> kv;
    std::size_t i = 0;
    const std::size_t n = t.size();
    auto ws = [&] { while (i < n && std::isspace(static_cast<unsigned char>(t[i]))) ++i; };
    auto str = [&]() {
        if (i >= n || t[i] != '"') throw std::invalid_argument("preset: expected a quoted string");
        ++i;
        std::string o;
        while (i < n && t[i] != '"') { if (t[i] == '\\' && i + 1 < n) ++i; o += t[i++]; }
        if (i >= n) throw std::invalid_argument("preset: unterminated string");
        ++i;
        return o;
    };
    ws();
    if (i >= n || t[i] != '{') throw std::invalid_argument("preset: expected '{'");
    ++i;
    ws();
    if (i < n && t[i] == '}') return kv;
    for (;;) {
        ws();
        const std::string key = str();
        ws();
        if (i >= n || t[i] != ':') throw std::invalid_argument("preset: expected ':' after \"" + key + "\"");
        ++i;
        ws();
        std::string val;
        if (i < n && t[i] == '"') {
            val = "s:" + str();
        } else {
            std::size_t j = i;
            while (j < n && t[j] != ',' && t[j] != '}' && !std::isspace(static_cast<unsigned char>(t[j]))) ++j;
            if (j == i) throw std::invalid_argument("preset: empty value for \"" + key + "\"");
            val = "v:" + t.substr(i, j - i);
            i = j;
        }
        kv[key] = val;
        ws();
        if (i < n && t[i] == ',') { ++i; continue; }
        if (i < n && t[i] == '}') break;
        throw std::invalid_argument("preset: expected ',' or '}'");
    }
    return kv;
}
}  // namespace detail

inline std::string preset_to_json(const Preset& pr) {
    using namespace detail;
    const SimulationParams& s = pr.sim;
    const CpiParams& c = pr.cpi;
    std::vector<std::pair<std::string, std::string>> kv;
    auto N = [&](const std::string& k, double v) { kv.emplace_back(k, jnum(v)); };
    auto T = [&](const std::string& k, const std::string& v) { kv.emplace_back(k, jstr(v)); };
    auto B = [&](const std::string& k, bool v) { kv.emplace_back(k, v ? "true" : "false"); };
    N("version", 1);
    T("mode", pr.cpi_mode ? "cpi" : "single");
    T("waveform", wf_name(s.wf.type));
    N("fs_hz", s.wf.fs); N("pri_s", s.wf.pri); N("tau_s", s.wf.tau); N("bandwidth_hz", s.wf.bandwidth);
    N("code_length", s.wf.code_length);
    N("carrier_hz", s.ch.fc); N("range_m", s.ch.range_m); N("velocity_mps", s.ch.velocity_mps);
    T("snr_mode", s.ch.snr_mode == SnrMode::Physical ? "physical" : "manual");
    N("snr_db", s.ch.snr_db); N("pt_w", s.ch.pt_w); N("gt_db", s.ch.gt_db); N("gr_db", s.ch.gr_db);
    N("rcs_m2", s.ch.rcs_m2); N("loss_db", s.ch.loss_db); N("noise_figure_db", s.ch.noise_figure_db);
    N("seed", static_cast<double>(s.ch.seed));
    B("phase_noise", s.phase_noise_enabled);
    N("pn_floor_dbc_hz", s.pn.floor_dbc_hz); N("pn_leeson_hz", s.pn.leeson_corner_hz);
    N("pn_flicker_hz", s.pn.flicker_corner_hz); N("pn_seed", static_cast<double>(s.pn.seed));
    T("mf_weighting", wt_name(s.mf_weighting));
    N("cfar_pfa", s.cfar.pfa); N("cfar_guard", s.cfar.guard_cells); N("cfar_train", s.cfar.train_cells);
    B("cfar_correlation_aware", s.cfar_correlation_aware);
    N("adc_bits", s.rx_imp.adc_bits); N("adc_fs_sigma", s.rx_imp.adc_fs_sigma);
    N("iq_gain_db", s.rx_imp.iq_gain_db); N("iq_phase_deg", s.rx_imp.iq_phase_deg);
    N("dc_i_sigma", s.rx_imp.dc_i_sigma); N("dc_q_sigma", s.rx_imp.dc_q_sigma);
    N("extra_count", static_cast<double>(s.extra_targets.size()));
    for (std::size_t i = 0; i < s.extra_targets.size(); ++i) {
        const std::string p = "extra" + std::to_string(i + 1) + "_";
        N(p + "range_m", s.extra_targets[i].range_m);
        N(p + "velocity_mps", s.extra_targets[i].velocity_mps);
        N(p + "rel_db", s.extra_targets[i].rel_db);
    }
    N("cpi_pulses", c.pulses); N("cpi_gate", c.gate_len); N("cpi_doppler_fft", c.doppler_fft);
    T("cpi_window", wt_name(c.window));
    N("cpi_pfa", c.cfar.pfa); N("cpi_guard_range", c.cfar.guard_range); N("cpi_guard_doppler", c.cfar.guard_doppler);
    N("cpi_train_range", c.cfar.train_range); N("cpi_train_doppler", c.cfar.train_doppler);
    N("psd_window", pr.psd_window);
    std::string out = "{\n";
    for (std::size_t i = 0; i < kv.size(); ++i)
        out += "  " + jstr(kv[i].first) + ": " + kv[i].second + (i + 1 < kv.size() ? ",\n" : "\n");
    return out + "}\n";
}

inline Preset preset_from_json(const std::string& text) {
    using namespace detail;
    const auto kv = parse_flat_json(text);
    auto find = [&](const std::string& k) -> const std::string* { auto it = kv.find(k); return it == kv.end() ? nullptr : &it->second; };
    auto num = [&](const std::string& k, double def) {
        const std::string* v = find(k);
        if (!v) return def;
        if ((*v)[0] != 'v') throw std::invalid_argument("preset: '" + k + "' must be a number");
        char* e = nullptr;
        const double x = std::strtod(v->c_str() + 2, &e);
        if (e == v->c_str() + 2 || *e != '\0') throw std::invalid_argument("preset: bad number for '" + k + "'");
        return x;
    };
    auto flag = [&](const std::string& k, bool def) {
        const std::string* v = find(k);
        if (!v) return def;
        if (*v == "v:true") return true;
        if (*v == "v:false") return false;
        throw std::invalid_argument("preset: '" + k + "' must be true or false");
    };
    auto str = [&](const std::string& k, const std::string& def) {
        const std::string* v = find(k);
        if (!v) return def;
        if ((*v)[0] != 's') throw std::invalid_argument("preset: '" + k + "' must be a string");
        return v->substr(2);
    };
    Preset pr;
    SimulationParams& s = pr.sim;
    CpiParams& c = pr.cpi;
    pr.cpi_mode = (str("mode", "single") == "cpi");
    s.wf.type = wf_from(str("waveform", wf_name(s.wf.type)));
    s.wf.fs = num("fs_hz", s.wf.fs); s.wf.pri = num("pri_s", s.wf.pri); s.wf.tau = num("tau_s", s.wf.tau);
    s.wf.bandwidth = num("bandwidth_hz", s.wf.bandwidth);
    s.wf.code_length = static_cast<int>(num("code_length", s.wf.code_length));
    s.ch.fc = num("carrier_hz", s.ch.fc); s.ch.range_m = num("range_m", s.ch.range_m);
    s.ch.velocity_mps = num("velocity_mps", s.ch.velocity_mps);
    const std::string mode = str("snr_mode", "manual");
    if (mode != "manual" && mode != "physical") throw std::invalid_argument("preset: unknown snr_mode '" + mode + "'");
    s.ch.snr_mode = (mode == "physical") ? SnrMode::Physical : SnrMode::Manual;
    s.ch.snr_db = num("snr_db", s.ch.snr_db); s.ch.pt_w = num("pt_w", s.ch.pt_w);
    s.ch.gt_db = num("gt_db", s.ch.gt_db); s.ch.gr_db = num("gr_db", s.ch.gr_db);
    s.ch.rcs_m2 = num("rcs_m2", s.ch.rcs_m2); s.ch.loss_db = num("loss_db", s.ch.loss_db);
    s.ch.noise_figure_db = num("noise_figure_db", s.ch.noise_figure_db);
    s.ch.seed = static_cast<std::uint64_t>(num("seed", static_cast<double>(s.ch.seed)));
    s.phase_noise_enabled = flag("phase_noise", s.phase_noise_enabled);
    s.pn.floor_dbc_hz = num("pn_floor_dbc_hz", s.pn.floor_dbc_hz);
    s.pn.leeson_corner_hz = num("pn_leeson_hz", s.pn.leeson_corner_hz);
    s.pn.flicker_corner_hz = num("pn_flicker_hz", s.pn.flicker_corner_hz);
    s.pn.seed = static_cast<std::uint64_t>(num("pn_seed", static_cast<double>(s.pn.seed)));
    s.mf_weighting = wt_from(str("mf_weighting", wt_name(s.mf_weighting)));
    s.cfar.pfa = num("cfar_pfa", s.cfar.pfa);
    s.cfar.guard_cells = static_cast<int>(num("cfar_guard", s.cfar.guard_cells));
    s.cfar.train_cells = static_cast<int>(num("cfar_train", s.cfar.train_cells));
    s.cfar_correlation_aware = flag("cfar_correlation_aware", s.cfar_correlation_aware);
    s.rx_imp.adc_bits = static_cast<int>(num("adc_bits", s.rx_imp.adc_bits));
    s.rx_imp.adc_fs_sigma = num("adc_fs_sigma", s.rx_imp.adc_fs_sigma);
    s.rx_imp.iq_gain_db = num("iq_gain_db", s.rx_imp.iq_gain_db);
    s.rx_imp.iq_phase_deg = num("iq_phase_deg", s.rx_imp.iq_phase_deg);
    s.rx_imp.dc_i_sigma = num("dc_i_sigma", s.rx_imp.dc_i_sigma);
    s.rx_imp.dc_q_sigma = num("dc_q_sigma", s.rx_imp.dc_q_sigma);
    const int ne = std::min(8, std::max(0, static_cast<int>(num("extra_count", 0))));
    for (int i = 1; i <= ne; ++i) {
        const std::string p = "extra" + std::to_string(i) + "_";
        TargetSpec t;
        t.range_m = num(p + "range_m", 0.0); t.velocity_mps = num(p + "velocity_mps", 0.0); t.rel_db = num(p + "rel_db", 0.0);
        s.extra_targets.push_back(t);
    }
    c.pulses = static_cast<int>(num("cpi_pulses", c.pulses)); c.gate_len = static_cast<int>(num("cpi_gate", c.gate_len));
    c.doppler_fft = static_cast<int>(num("cpi_doppler_fft", c.doppler_fft));
    c.window = wt_from(str("cpi_window", wt_name(c.window)));
    c.cfar.pfa = num("cpi_pfa", c.cfar.pfa);
    c.cfar.guard_range = static_cast<int>(num("cpi_guard_range", c.cfar.guard_range));
    c.cfar.guard_doppler = static_cast<int>(num("cpi_guard_doppler", c.cfar.guard_doppler));
    c.cfar.train_range = static_cast<int>(num("cpi_train_range", c.cfar.train_range));
    c.cfar.train_doppler = static_cast<int>(num("cpi_train_doppler", c.cfar.train_doppler));
    pr.psd_window = static_cast<int>(num("psd_window", pr.psd_window));
    return pr;
}

inline void save_preset(const std::string& path, const Preset& pr) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << preset_to_json(pr);
    if (!f) throw std::runtime_error("write failed: " + path);
}

inline Preset load_preset(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot read " + path);
    std::stringstream ss;
    ss << f.rdbuf();
    return preset_from_json(ss.str());
}

}  // namespace radar
