#pragma once
#include <QMetaType>
#include <QString>
#include <memory>
#include "radar/minmax_pyramid.hpp"
#include "radar/peak_detect.hpp"
#include "radar/psd.hpp"
#include "radar/simulation.hpp"

struct RunRequest {
    radar::SimulationParams sim;
    radar::Window psd_window = radar::Window::Hann;
};

enum class ExportFormat { Hdf5 = 0, BinFloat32 = 1, BinInt16 = 2, Csv = 3 };   // same order as the dropdown

struct ExportRequest {
    radar::ResultPtr result;
    QString path;
    ExportFormat format = ExportFormat::Hdf5;
    bool tx = false;   // binary / CSV: export S_tx instead of S_rx
};

// Everything the UI needs for one run; immutable once published.
struct PlotData {
    radar::ResultPtr result;            // keeps the borrowed I/Q views below alive
    radar::MinMaxPyramid trace;         // matched-filter output, dB re MF noise floor
    double ref = 1.0;                   // linear power of the MF noise floor (Pn / E)
    radar::MinMaxPyramid rx_i, rx_q;    // borrowed views into result->rx
    radar::MinMaxPyramid tx_i, tx_q;    // borrowed views into result->tx
    std::vector<radar::Peak> peaks;     // peak detection on the matched-filter profile, strongest first
    double peak_threshold_db = 0.0;     // detection threshold used (dB above the noise floor)
    double sigma = 1.0;                 // noise std per I/Q component
    radar::MinMaxPyramid psd_ideal;     // transmit PSD, dB re ideal peak
    radar::MinMaxPyramid psd_impaired;  // same with LO phase noise (empty when disabled)
    double psd_df = 0.0, psd_f_start = 0.0;
};
using PlotDataPtr = std::shared_ptr<const PlotData>;

Q_DECLARE_METATYPE(radar::SimulationParams)
Q_DECLARE_METATYPE(RunRequest)
Q_DECLARE_METATYPE(ExportRequest)
Q_DECLARE_METATYPE(PlotDataPtr)
