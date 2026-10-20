#pragma once
#include <hdf5.h>
#include <stdexcept>
#include <string>
#include <vector>
#include "radar/cpi.hpp"
#include "radar/simulation.hpp"

namespace radar {

class H5Handle {
public:
    H5Handle(hid_t id, herr_t (*closer)(hid_t), const char* what) : id_(id), closer_(closer) {
        if (id_ < 0) throw std::runtime_error(std::string("HDF5 error: ") + what);
    }
    ~H5Handle() { if (id_ >= 0) closer_(id_); }
    H5Handle(const H5Handle&) = delete;
    H5Handle& operator=(const H5Handle&) = delete;
    H5Handle(H5Handle&& o) noexcept : id_(o.id_), closer_(o.closer_) { o.id_ = -1; }
    hid_t get() const { return id_; }
private:
    hid_t id_;
    herr_t (*closer_)(hid_t);
};

// Compound {r: float32, i: float32}: same memory layout as std::complex<float>
inline H5Handle h5_complex_type() {
    H5Handle t(H5Tcreate(H5T_COMPOUND, sizeof(cf32)), H5Tclose, "create complex type");
    if (H5Tinsert(t.get(), "r", 0, H5T_NATIVE_FLOAT) < 0 ||
        H5Tinsert(t.get(), "i", sizeof(float), H5T_NATIVE_FLOAT) < 0)
        throw std::runtime_error("HDF5 error: complex type fields");
    return t;
}

inline void h5_write_dataset(hid_t loc, const char* name, hid_t type, const void* data, std::size_t n) {
    hsize_t dims[1] = {static_cast<hsize_t>(n)};
    H5Handle space(H5Screate_simple(1, dims, nullptr), H5Sclose, "dataspace");
    H5Handle ds(H5Dcreate2(loc, name, type, space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose, name);
    if (H5Dwrite(ds.get(), type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0)
        throw std::runtime_error(std::string("HDF5 error: writing ") + name);
}

inline void h5_attr_double(hid_t loc, const char* name, double v) {
    H5Handle sp(H5Screate(H5S_SCALAR), H5Sclose, "scalar space");
    H5Handle a(H5Acreate2(loc, name, H5T_NATIVE_DOUBLE, sp.get(), H5P_DEFAULT, H5P_DEFAULT), H5Aclose, name);
    if (H5Awrite(a.get(), H5T_NATIVE_DOUBLE, &v) < 0) throw std::runtime_error(std::string("HDF5 error: attr ") + name);
}

inline void h5_attr_int(hid_t loc, const char* name, long long v) {
    H5Handle sp(H5Screate(H5S_SCALAR), H5Sclose, "scalar space");
    H5Handle a(H5Acreate2(loc, name, H5T_NATIVE_LLONG, sp.get(), H5P_DEFAULT, H5P_DEFAULT), H5Aclose, name);
    if (H5Awrite(a.get(), H5T_NATIVE_LLONG, &v) < 0) throw std::runtime_error(std::string("HDF5 error: attr ") + name);
}

inline void h5_attr_str(hid_t loc, const char* name, const std::string& s) {
    H5Handle t(H5Tcopy(H5T_C_S1), H5Tclose, "string type");
    if (H5Tset_size(t.get(), s.size() + 1) < 0) throw std::runtime_error("HDF5 error: string size");
    H5Handle sp(H5Screate(H5S_SCALAR), H5Sclose, "scalar space");
    H5Handle a(H5Acreate2(loc, name, t.get(), sp.get(), H5P_DEFAULT, H5P_DEFAULT), H5Aclose, name);
    if (H5Awrite(a.get(), t.get(), s.c_str()) < 0) throw std::runtime_error(std::string("HDF5 error: attr ") + name);
}

// Datasets: /tx_clean, /rx (complex64), /mf_power (float32, not for CW), /detections_bin (if any).
// All physical metadata are attributes on the root group.
inline void write_hdf5(const std::string& path, const SimulationResult& r) {
    H5Handle file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose, "create file");
    const hid_t f = file.get();
    H5Handle ct = h5_complex_type();
    h5_write_dataset(f, "tx_clean", ct.get(), r.tx.data(), r.tx.size());
    h5_write_dataset(f, "rx", ct.get(), r.rx.data(), r.rx.size());
    if (r.has_range_profile) h5_write_dataset(f, "mf_power", H5T_NATIVE_FLOAT, r.mf_power.data(), r.mf_power.size());
    if (!r.cfar.detections.empty()) {
        std::vector<unsigned long long> d(r.cfar.detections.begin(), r.cfar.detections.end());
        h5_write_dataset(f, "detections_bin", H5T_NATIVE_ULLONG, d.data(), d.size());
    }

    const auto& p = r.params;
    const auto& rep = r.report;
    h5_attr_str(f, "waveform", waveform_name(p.wf.type));
    h5_attr_str(f, "snr_mode", p.ch.snr_mode == SnrMode::Physical ? "physical" : "manual");
    h5_attr_str(f, "units", "amplitude in sqrt(watts); |x|^2 = watts");
    h5_attr_double(f, "sample_rate_hz", p.wf.fs);
    h5_attr_double(f, "carrier_hz", p.ch.fc);
    h5_attr_double(f, "pri_s", p.wf.pri);
    h5_attr_double(f, "prf_hz", 1.0 / p.wf.pri);
    h5_attr_double(f, "pulse_width_s", p.wf.tau);
    h5_attr_double(f, "sweep_bandwidth_hz", p.wf.bandwidth);
    h5_attr_int(f, "code_length", p.wf.code_length);
    h5_attr_double(f, "target_range_m", p.ch.range_m);
    h5_attr_double(f, "target_velocity_mps", p.ch.velocity_mps);
    h5_attr_double(f, "doppler_hz", rep.doppler_hz);
    h5_attr_double(f, "delay_s", rep.td_s);
    h5_attr_double(f, "delay_samples", rep.delay_samples);
    h5_attr_double(f, "r_unambiguous_m", rep.r_unambiguous_m);
    h5_attr_double(f, "r_blind_m", rep.r_blind_m);
    h5_attr_double(f, "rx_pulse_power_w", rep.pr_w);
    h5_attr_double(f, "noise_power_per_sample_w", rep.pn_w);
    h5_attr_double(f, "snr_input_db", rep.snr_in_db);
    h5_attr_double(f, "mf_gain_db", rep.mf_gain_db);
    h5_attr_double(f, "noise_figure_db", p.ch.noise_figure_db);
    h5_attr_int(f, "phase_noise_enabled", p.phase_noise_enabled ? 1 : 0);
    h5_attr_double(f, "phase_noise_floor_dbc_hz", p.pn.floor_dbc_hz);
    h5_attr_int(f, "noise_seed", static_cast<long long>(p.ch.seed));
    h5_attr_int(f, "num_samples", static_cast<long long>(r.rx.size()));
    h5_attr_int(f, "num_pulse_samples", static_cast<long long>(r.n_pulse));
}

inline void h5_write_dataset2d(hid_t loc, const char* name, hid_t type, const void* data, std::size_t rows, std::size_t cols) {
    hsize_t dims[2] = {static_cast<hsize_t>(rows), static_cast<hsize_t>(cols)};
    H5Handle space(H5Screate_simple(2, dims, nullptr), H5Sclose, "dataspace");
    H5Handle ds(H5Dcreate2(loc, name, type, space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose, name);
    if (H5Dwrite(ds.get(), type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0)
        throw std::runtime_error(std::string("HDF5 error: writing ") + name);
}

// Range-Doppler map (/rd_power, rows = Doppler fftshifted with row doppler_rows/2 = 0 Hz, columns = range gate bins),
// CFAR detections (/det_*), and all CPI metadata as attributes.
inline void write_hdf5_cpi(const std::string& path, const CpiResult& r) {
    H5Handle file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose, "create file");
    const hid_t f = file.get();
    h5_write_dataset2d(f, "rd_power", H5T_NATIVE_FLOAT, r.rd_power.data(), r.doppler_rows, r.gate_len);
    if (!r.detections.empty()) {
        std::vector<double> rng, rngr, dop, dopr, vel, velr, snr;
        for (const auto& d : r.detections) {
            rng.push_back(d.range_m); rngr.push_back(d.range_refined_m);
            dop.push_back(d.doppler_hz); dopr.push_back(d.doppler_refined_hz);
            vel.push_back(d.velocity_mps); velr.push_back(d.velocity_refined_mps);
            snr.push_back(d.snr_db);
        }
        h5_write_dataset(f, "det_range_m", H5T_NATIVE_DOUBLE, rng.data(), rng.size());
        h5_write_dataset(f, "det_range_refined_m", H5T_NATIVE_DOUBLE, rngr.data(), rngr.size());
        h5_write_dataset(f, "det_doppler_hz", H5T_NATIVE_DOUBLE, dop.data(), dop.size());
        h5_write_dataset(f, "det_doppler_refined_hz", H5T_NATIVE_DOUBLE, dopr.data(), dopr.size());
        h5_write_dataset(f, "det_velocity_mps", H5T_NATIVE_DOUBLE, vel.data(), vel.size());
        h5_write_dataset(f, "det_velocity_refined_mps", H5T_NATIVE_DOUBLE, velr.data(), velr.size());
        h5_write_dataset(f, "det_snr_db", H5T_NATIVE_DOUBLE, snr.data(), snr.size());
    }
    const auto& p = r.sim;
    const auto& rep = r.report;
    h5_attr_str(f, "waveform", waveform_name(p.wf.type));
    h5_attr_str(f, "units", "rd_power = |x|^2 in watts, normalized so an on-bin target peak equals its single-pulse compressed power");
    h5_attr_double(f, "sample_rate_hz", p.wf.fs);
    h5_attr_double(f, "carrier_hz", p.ch.fc);
    h5_attr_double(f, "pri_s", p.wf.pri);
    h5_attr_double(f, "prf_hz", rep.prf_hz);
    h5_attr_int(f, "pulses", r.cpi.pulses);
    h5_attr_int(f, "gate_start_bin", static_cast<long long>(r.gate_start));
    h5_attr_int(f, "gate_len", static_cast<long long>(r.gate_len));
    h5_attr_int(f, "doppler_rows", static_cast<long long>(r.doppler_rows));
    h5_attr_int(f, "zero_doppler_row", static_cast<long long>(r.doppler_rows / 2));
    h5_attr_double(f, "doppler_resolution_hz", rep.doppler_resolution_hz);
    h5_attr_double(f, "velocity_resolution_mps", rep.velocity_resolution_mps);
    h5_attr_double(f, "v_unambiguous_mps", rep.v_unambiguous_mps);
    h5_attr_double(f, "rd_noise_floor_w", r.rd_noise_floor_w);
    h5_attr_double(f, "input_snr_db", r.first_pulse.snr_in_db);
    h5_attr_double(f, "mf_gain_db", r.first_pulse.mf_gain_db);
    h5_attr_double(f, "slow_time_gain_db", rep.integration_gain_db);
    h5_attr_double(f, "target_range_m", p.ch.range_m);
    h5_attr_double(f, "target_velocity_mps", p.ch.velocity_mps);
    h5_attr_int(f, "extra_targets", static_cast<long long>(p.extra_targets.size()));
    h5_attr_int(f, "noise_seed", static_cast<long long>(p.ch.seed));
}

}  // namespace radar
