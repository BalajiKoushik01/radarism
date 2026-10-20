#include "SimWorker.h"
#include <algorithm>
#include <cmath>
#include "radar/export_hdf5.hpp"
#include "radar/exporters.hpp"

namespace {
radar::RBuffer to_db(const radar::RBuffer& v, double reference) {
    radar::RBuffer out(v.size());
    const double inv = 1.0 / reference;
    for (std::size_t i = 0; i < v.size(); ++i)
        out[i] = static_cast<float>(10.0 * std::log10(std::max(static_cast<double>(v[i]) * inv, 1e-14)));
    return out;
}
}  // namespace

void SimWorker::run(const RunRequest& rq) {
    cancel_.store(false);
    try {
        radar::ResultPtr res = radar::simulate(rq.sim, &cancel_);
        auto pd = std::make_shared<PlotData>();
        pd->result = res;
        const std::size_t n = res->rx.size();
        const double fs = rq.sim.wf.fs;

        // Matched-filter trace in dB relative to its noise floor (empty for CW)
        if (res->has_range_profile) {
            pd->ref = res->mf_noise_floor_w;
            pd->trace = radar::MinMaxPyramid(to_db(res->mf_power, pd->ref));
        }

        // Borrowed I/Q views (interleaved complex floats, stride 2)
        const float* rx = reinterpret_cast<const float*>(res->rx.data());
        const float* tx = reinterpret_cast<const float*>(res->tx.data());
        pd->rx_i = radar::MinMaxPyramid(rx, n, 2);
        pd->rx_q = radar::MinMaxPyramid(rx + 1, n, 2);
        pd->tx_i = radar::MinMaxPyramid(tx, n, 2);
        pd->tx_q = radar::MinMaxPyramid(tx + 1, n, 2);
        pd->sigma = std::sqrt(res->report.pn_w / 2.0);

        // Transmit spectra, both normalized to the ideal spectrum's peak
        radar::Psd ideal = radar::compute_pulse_psd(res->tx, res->n_pulse, fs, rq.psd_window);
        double peak = static_cast<double>(*std::max_element(ideal.density.begin(), ideal.density.end()));
        if (!(peak > 0.0)) peak = 1.0;
        pd->psd_df = ideal.df;
        pd->psd_f_start = ideal.f_start;
        pd->psd_ideal = radar::MinMaxPyramid(to_db(ideal.density, peak));
        if (!res->tx_impaired.empty()) {
            radar::Psd imp = radar::compute_pulse_psd(res->tx_impaired, res->n_pulse, fs, rq.psd_window);
            pd->psd_impaired = radar::MinMaxPyramid(to_db(imp.density, peak));
        }
        emit finished(pd);
    } catch (const radar::Cancelled&) {
        emit failed("Cancelled.");
    } catch (const std::exception& e) {
        emit failed(QString::fromUtf8(e.what()));
    }
}

void SimWorker::exportData(const ExportRequest& rq) {
    try {
        const radar::SimulationResult& r = *rq.result;
        const std::string path = rq.path.toStdString();
        const double fs = r.params.wf.fs, fc = r.params.ch.fc;
        QString msg = rq.path;
        switch (rq.format) {
            case ExportFormat::Hdf5:
                radar::write_hdf5(path, r);
                msg += " (tx_clean, rx, mf_power + metadata attributes)";
                break;
            case ExportFormat::BinFloat32:
            case ExportFormat::BinInt16: {
                const auto fmt = (rq.format == ExportFormat::BinInt16) ? radar::BinFormat::Int16 : radar::BinFormat::Float32;
                const auto info = radar::write_bin(path, rq.tx ? r.tx : r.rx, fmt, fs, fc,
                                                   rq.tx ? "tx_clean_for_awg" : "rx_for_sdr_dsp");
                msg += QString(" + .json sidecar (scale %1 stored units per sqrt-W)").arg(info.scale, 0, 'g', 6);
                break;
            }
            case ExportFormat::Csv:
                radar::write_csv(path, rq.tx ? r.tx : r.rx, fs);
                break;
        }
        emit exportDone(msg);
    } catch (const std::exception& e) {
        emit exportFailed(QString::fromUtf8(e.what()));
    }
}

void SimWorker::runCpi(const CpiRequest& rq) {
    cancel_.store(false);
    try {
        auto res = radar::simulate_cpi(rq.sim, rq.cpi, &cancel_, [this](int k, int m) {
            if (k == m || k % 8 == 0) emit cpiProgress(k, m);
        });
        emit cpiFinished(res);
    } catch (const radar::Cancelled&) {
        emit failed("Cancelled.");
    } catch (const std::exception& e) {
        emit failed(QString::fromUtf8(e.what()));
    }
}

void SimWorker::runMc(const McRequest& rq) {
    cancel_.store(false);
    try {
        auto res = std::make_shared<const radar::McResult>(radar::run_monte_carlo(rq.mc, rq.snr_out_db, &cancel_));
        emit mcFinished(res);
    } catch (const radar::Cancelled&) {
        emit failed("Cancelled.");
    } catch (const std::exception& e) {
        emit failed(QString::fromUtf8(e.what()));
    }
}

void SimWorker::exportCpi(const CpiExportRequest& rq) {
    try {
        radar::write_hdf5_cpi(rq.path.toStdString(), *rq.result);
        emit exportDone(rq.path + " (rd_power, detections and metadata)");
    } catch (const std::exception& e) {
        emit exportFailed(QString::fromUtf8(e.what()));
    }
}
