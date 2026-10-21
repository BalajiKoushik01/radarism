#include "MainWindow.h"
#include <QCheckBox>
#include <QComboBox>
#include <QFontDatabase>
#include <QPlainTextEdit>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QStackedWidget>
#include <QStandardItemModel>
#include <QThread>
#include <QVBoxLayout>
#include <QWidget>
#include <algorithm>
#include <cmath>
#include "PlotWidget.h"
#include "RdPlot.h"
#include "SeriesPlot.h"
#include "SimWorker.h"
#include "radar/exporters.hpp"
#include "radar/passband.hpp"

namespace {
QDoubleSpinBox* spin(double lo, double hi, double val, int dec, const QString& suffix) {
    auto* s = new QDoubleSpinBox;
    s->setRange(lo, hi);
    s->setDecimals(dec);
    s->setValue(val);
    s->setSuffix(suffix);
    s->setKeyboardTracking(false);
    return s;
}
QSpinBox* ispin(int lo, int hi, int val) {
    auto* s = new QSpinBox;
    s->setRange(lo, hi);
    s->setValue(val);
    s->setKeyboardTracking(false);
    return s;
}
} // namespace

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    setWindowTitle("RadarSim");
    stack_ = new QStackedWidget;
    stack_->addWidget(buildSetupPage());
    stack_->addWidget(buildDashboardPage());
    stack_->addWidget(buildCpiPage());
    stack_->addWidget(buildMcPage());
    setCentralWidget(stack_);

    thread_ = new QThread(this);
    worker_ = new SimWorker;
    worker_->moveToThread(thread_);
    connect(thread_, &QThread::finished, worker_, &QObject::deleteLater);
    connect(this, &MainWindow::requestRun, worker_, &SimWorker::run);
    connect(worker_, &SimWorker::finished, this, &MainWindow::onFinished);
    connect(worker_, &SimWorker::failed, this, &MainWindow::onFailed);
    connect(this, &MainWindow::requestExport, worker_, &SimWorker::exportData);
    connect(worker_, &SimWorker::exportDone, this, &MainWindow::onExportDone);
    connect(worker_, &SimWorker::exportFailed, this, &MainWindow::onExportFailed);
    connect(this, &MainWindow::requestCpi, worker_, &SimWorker::runCpi);
    connect(worker_, &SimWorker::cpiFinished, this, &MainWindow::onCpiFinished);
    connect(worker_, &SimWorker::cpiProgress, this, &MainWindow::onCpiProgress);
    connect(this, &MainWindow::requestMc, worker_, &SimWorker::runMc);
    connect(worker_, &SimWorker::mcFinished, this, &MainWindow::onMcFinished);
    connect(this, &MainWindow::requestCpiExport, worker_, &SimWorker::exportCpi);
    thread_->start();

    updateForm();
}

MainWindow::~MainWindow() {
    worker_->cancel();
    thread_->quit();
    thread_->wait();
}

QWidget* MainWindow::buildSetupPage() {
    type_ = new QComboBox;
    type_->addItems({"CW (Doppler / calibration)", "Unmodulated pulse", "LFM (chirp)", "Phase code: Barker",
                     "Phase code: Frank"});
    type_->setCurrentIndex(2);

    fc_ = spin(0.001, 300.0, 10.0, 3, " GHz");
    fbb_ = spin(-500.0, 500.0, 0.0, 3, " MHz");
    bwBb_ = spin(0.001, 1000.0, 1.0, 3, " MHz");
    fs_ = spin(0.001, 2000.0, 2.0, 3, " MHz");

    pri_ = spin(0.001, 100000.0, 10.0, 3, " ms");
    tau_ = spin(0.01, 1e7, 500.0, 2, " µs");
    bw_ = spin(0.1, 1e6, 200.0, 1, " kHz");
    code_ = new QComboBox;

    range_ = spin(0.1, 5000.0, 150.0, 2, " km");
    dynamics_ = new QComboBox;
    dynamics_->addItems({"Stationary", "Moving"});
    vel_ = spin(-5000.0, 5000.0, 100.0, 1, " m/s");
    snrMode_ = new QComboBox;
    snrMode_->addItems({"Manual SNR", "Physical (radar equation)"});
    snr_ = spin(-80.0, 80.0, -10.0, 1, " dB");
    pt_ = spin(0.001, 1e9, 1000.0, 3, " W");
    gt_ = spin(0.0, 80.0, 30.0, 1, " dBi");
    gr_ = spin(0.0, 80.0, 30.0, 1, " dBi");
    rcs_ = spin(0.0001, 1e6, 1.0, 4, " m²");
    loss_ = spin(0.0, 60.0, 0.0, 1, " dB");
    nf_ = spin(0.0, 30.0, 3.0, 1, " dB");

    noiseProfile_ = new QComboBox;
    noiseProfile_->addItems({"AWGN only", "AWGN + Tx phase noise (Leeson)"});
    pnFloor_ = spin(-160.0, -20.0, -70.0, 1, " dBc/Hz");
    pnLeeson_ = spin(0.0, 10000.0, 1.0, 3, " kHz");
    pnFlicker_ = spin(0.0, 1e6, 100.0, 1, " Hz");
    psdWindow_ = new QComboBox;
    psdWindow_->addItems({"Hann", "Rectangular"});
    mfWin_ = new QComboBox;
    mfWin_->addItems({"None (matched filter)", "Hamming (about -1.3 dB SNR)", "Hann (about -1.8 dB SNR)"});
    pfaExp_ = spin(1.0, 15.0, 6.0, 1, "");
    pfaExp_->setPrefix("1e-");
    guard_ = ispin(0, 2000, 20);
    train_ = ispin(4, 5000, 64);
    buildCpiControls();
    buildExtraControls();

    auto* signalBox = new QGroupBox("Signal and Baseband Architecture");
    signalForm_ = new QFormLayout(signalBox);
    signalForm_->addRow("Processing", mode_);
    signalForm_->addRow("Waveform", type_);
    signalForm_->addRow("Carrier frequency (RF)", fc_);
    signalForm_->addRow("Baseband IF offset (f_bb)", fbb_);
    signalForm_->addRow("Baseband filter BW (B_bb)", bwBb_);
    signalForm_->addRow("ADC sampling rate (f_s)", fs_);
    signalForm_->addRow("PRI", pri_);
    signalForm_->addRow("Pulse width τ", tau_);
    signalForm_->addRow("Sweep bandwidth", bw_);
    signalForm_->addRow("Code length", code_);

    auto* targetBox = new QGroupBox("Target and receiver");
    targetForm_ = new QFormLayout(targetBox);
    targetForm_->addRow("Target range", range_);
    targetForm_->addRow("Target dynamics", dynamics_);
    targetForm_->addRow("Radial velocity", vel_);
    targetForm_->addRow("SNR control", snrMode_);
    targetForm_->addRow("Input SNR", snr_);
    targetForm_->addRow("Tx power", pt_);
    targetForm_->addRow("Tx antenna gain", gt_);
    targetForm_->addRow("Rx antenna gain", gr_);
    targetForm_->addRow("Target RCS", rcs_);
    targetForm_->addRow("System loss", loss_);
    targetForm_->addRow("Receiver noise figure", nf_);

    auto* noiseBox = new QGroupBox("Noise, impairments and detection");
    noiseForm_ = new QFormLayout(noiseBox);
    noiseForm_->addRow("Noise profile", noiseProfile_);
    noiseForm_->addRow("Phase-noise floor", pnFloor_);
    noiseForm_->addRow("Leeson corner", pnLeeson_);
    noiseForm_->addRow("Flicker corner", pnFlicker_);
    noiseForm_->addRow("PSD window", psdWindow_);
    noiseForm_->addRow("Matched-filter weighting", mfWin_);
    cfarBox_ = new QGroupBox("CA-CFAR detector");
    auto* cf = new QFormLayout(cfarBox_);
    cf->addRow("False-alarm probability", pfaExp_);
    cf->addRow("Guard cells (each side)", guard_);
    cf->addRow("Training cells (each side)", train_);
    noiseForm_->addRow(cfarBox_);

    auto* cols = new QHBoxLayout;
    cols->addWidget(signalBox, 1, Qt::AlignTop);
    cols->addWidget(targetBox, 1, Qt::AlignTop);
    cols->addWidget(noiseBox, 1, Qt::AlignTop);

    validation_ = new QLabel;
    validation_->setWordWrap(true);
    info_ = new QLabel;
    info_->setWordWrap(true);
    generate_ = new QPushButton("GENERATE");
    generate_->setMinimumHeight(40);
    QFont gf = generate_->font();
    gf.setBold(true);
    generate_->setFont(gf);
    cancel_ = new QPushButton("Cancel");
    cancel_->setMinimumHeight(40);
    cancel_->setEnabled(false);
    auto* buttons = new QHBoxLayout;
    buttons->addWidget(generate_, 3);
    buttons->addWidget(cancel_, 1);
    buttons->addWidget(savePreset_);
    buttons->addWidget(loadPreset_);
    buttons->addWidget(mcBtn_);

    auto* page = new QWidget;
    auto* root = new QVBoxLayout(page);
    root->addLayout(cols);
    root->addWidget(cpiBox_);
    root->addWidget(extrasRow_);
    root->addWidget(validation_);
    root->addWidget(info_);
    root->addLayout(buttons);
    root->addStretch(1);

    connect(type_, &QComboBox::currentIndexChanged, this, [this](int) { populateCodes(); });
    populateCodes();
    for (auto* w : page->findChildren<QDoubleSpinBox*>())
        connect(w, &QDoubleSpinBox::valueChanged, this, [this](double) { updateForm(); });
    for (auto* w : page->findChildren<QSpinBox*>())
        connect(w, &QSpinBox::valueChanged, this, [this](int) { updateForm(); });
    for (auto* w : page->findChildren<QComboBox*>())
        connect(w, &QComboBox::currentIndexChanged, this, [this](int) { updateForm(); });
    for (auto* w : page->findChildren<QCheckBox*>())
        connect(w, &QCheckBox::toggled, this, [this](bool) { updateForm(); });
    connect(generate_, &QPushButton::clicked, this, &MainWindow::onGenerate);
    connect(cancel_, &QPushButton::clicked, this, [this] { worker_->cancel(); });
    return page;
}

QWidget* MainWindow::buildDashboardPage() {
    back_ = new QPushButton("← Setup");
    status_ = new QLabel("No data yet.");
    peakInfo_ = new QLabel;
    peakInfo_->setTextFormat(Qt::RichText);
    peakInfo_->setWordWrap(true);
    {
        QFont pf = peakInfo_->font();
        pf.setPointSizeF(pf.pointSizeF() + 3.0);
        peakInfo_->setFont(pf);
    }
    status_->setWordWrap(true);
    iqSource_ = new QComboBox;
    iqSource_->addItems({"Rx (noisy)", "Tx (ideal)"});
    plot_ = new PlotWidget;
    iq_ = new SeriesPlot;
    psd_ = new SeriesPlot;
    pb_ = new SeriesPlot;
    pbCenter_ = spin(0.0, 1e8, 0.0, 6, " µs");
    pbWidth_ = spin(0.001, 1e6, 2.0, 3, " ns");
    pbInfo_ = new QLabel;
    expBuffer_ = new QComboBox;
    expBuffer_->addItems({"Rx buffer (S_rx)", "Tx buffer (S_tx, for AWG)"});
    expBuffer_->setEnabled(false);
    expFormat_ = new QComboBox;
    expFormat_->addItems({"HDF5 (.h5): Tx + Rx + metadata", "Binary float32 (.bin + JSON)",
                          "Binary int16 (.bin + JSON)", "CSV (.csv)"});
    export_ = new QPushButton("Export…");
    export_->setEnabled(false);

    auto* top = new QHBoxLayout;
    top->addWidget(back_);
    top->addWidget(status_, 1);
    top->addWidget(new QLabel("I/Q source"));
    top->addWidget(iqSource_);
    top->addSpacing(16);
    top->addWidget(new QLabel("Export"));
    top->addWidget(expBuffer_);
    top->addWidget(expFormat_);
    top->addWidget(export_);

    auto* pbBar = new QHBoxLayout;
    pbBar->addWidget(new QLabel("Passband view (computed lazily for the visible window only): centre"));
    pbBar->addWidget(pbCenter_);
    pbBar->addWidget(new QLabel("window"));
    pbBar->addWidget(pbWidth_);
    pbBar->addWidget(pbInfo_, 1);

    auto* plots = new QGridLayout;
    plots->addWidget(plot_, 0, 0, 1, 2);
    plots->addWidget(iq_, 1, 0);
    plots->addWidget(psd_, 1, 1);
    plots->addLayout(pbBar, 2, 0, 1, 2);
    plots->addWidget(pb_, 3, 0, 1, 2);
    plots->setRowStretch(0, 3);
    plots->setRowStretch(1, 2);
    plots->setRowStretch(3, 2);

    auto* page = new QWidget;
    auto* root = new QVBoxLayout(page);
    root->addLayout(top);
    root->addWidget(peakInfo_);
    root->addLayout(plots, 1);

    connect(back_, &QPushButton::clicked, this, [this] { stack_->setCurrentIndex(0); });
    connect(iqSource_, &QComboBox::currentIndexChanged, this, [this](int) { refreshIq(); resetPassbandView(); });
    connect(pbCenter_, &QDoubleSpinBox::valueChanged, this, [this](double) { refreshPassband(); });
    connect(pbWidth_, &QDoubleSpinBox::valueChanged, this, [this](double) { refreshPassband(); });
    connect(expFormat_, &QComboBox::currentIndexChanged, this, [this](int i) { expBuffer_->setEnabled(i != 0); });
    connect(export_, &QPushButton::clicked, this, &MainWindow::onExport);
    return page;
}

void MainWindow::populateCodes() {
    const int t = type_->currentIndex();
    QSignalBlocker block(code_);
    code_->clear();
    if (t == 3) {
        for (int L : {2, 3, 4, 5, 7, 11, 13}) code_->addItem(QString::number(L), L);
        code_->setCurrentIndex(6);
    } else if (t == 4) {
        for (int L : {4, 9, 16, 25, 36, 49, 64}) code_->addItem(QString::number(L), L);
        code_->setCurrentIndex(2);
    }
}

radar::SimulationParams MainWindow::buildParams() const {
    using radar::WaveformType;
    radar::SimulationParams p;
    switch (type_->currentIndex()) {
        case 0: p.wf.type = WaveformType::CW; break;
        case 1: p.wf.type = WaveformType::Pulse; break;
        case 2: p.wf.type = WaveformType::LFM; break;
        case 3: p.wf.type = WaveformType::Barker; break;
        default: p.wf.type = WaveformType::Frank; break;
    }
    p.wf.fs = fs_->value() * 1e6;
    p.wf.pri = pri_->value() * 1e-3;
    p.wf.tau = tau_->value() * 1e-6;
    p.wf.bandwidth = bw_->value() * 1e3;
    p.wf.code_length = code_->currentData().isValid() ? code_->currentData().toInt() : 13;

    p.ch.fc = fc_->value() * 1e9;
    p.ch.range_m = range_->value() * 1e3;
    p.ch.velocity_mps = (dynamics_->currentIndex() == 1) ? vel_->value() : 0.0;
    p.ch.snr_mode = (snrMode_->currentIndex() == 1) ? radar::SnrMode::Physical : radar::SnrMode::Manual;
    p.ch.snr_db = snr_->value();
    p.ch.pt_w = pt_->value();
    p.ch.gt_db = gt_->value();
    p.ch.gr_db = gr_->value();
    p.ch.rcs_m2 = rcs_->value();
    p.ch.loss_db = loss_->value();
    p.ch.noise_figure_db = nf_->value();

    p.phase_noise_enabled = (noiseProfile_->currentIndex() == 1);
    p.pn.floor_dbc_hz = pnFloor_->value();
    p.pn.leeson_corner_hz = pnLeeson_->value() * 1e3;
    p.pn.flicker_corner_hz = pnFlicker_->value();

    p.cfar.pfa = std::pow(10.0, -pfaExp_->value());
    p.cfar.guard_cells = guard_->value();
    p.cfar.train_cells = train_->value();
    p.run_detection = (p.wf.type != WaveformType::CW);
    fillExtras(p);
    p.cfar_correlation_aware = true;   // threshold sized for the matched filter's correlated noise
    p.mf_weighting = (mfWin_->currentIndex() == 1)   ? radar::MfWeighting::Hamming
                     : (mfWin_->currentIndex() == 2) ? radar::MfWeighting::Hann
                                                     : radar::MfWeighting::None;
    return p;
}

void MainWindow::updateForm() {
    updateVisibility();
    updateValidation();
}

void MainWindow::updateVisibility() {
    const int t = type_->currentIndex();
    const bool cw = (t == 0), lfm = (t == 2), coded = (t >= 3);
    signalForm_->setRowVisible(tau_, !cw);
    signalForm_->setRowVisible(bw_, lfm);
    signalForm_->setRowVisible(code_, coded);
    targetForm_->setRowVisible(vel_, dynamics_->currentIndex() == 1);
    const bool phys = (snrMode_->currentIndex() == 1);
    targetForm_->setRowVisible(snr_, !phys);
    for (QWidget* w : {static_cast<QWidget*>(pt_), static_cast<QWidget*>(gt_), static_cast<QWidget*>(gr_),
                       static_cast<QWidget*>(rcs_), static_cast<QWidget*>(loss_)})
        targetForm_->setRowVisible(w, phys);
    const bool pn = (noiseProfile_->currentIndex() == 1);
    for (QWidget* w : {static_cast<QWidget*>(pnFloor_), static_cast<QWidget*>(pnLeeson_),
                       static_cast<QWidget*>(pnFlicker_)})
        noiseForm_->setRowVisible(w, pn);
    const bool cpiMode = (mode_->currentIndex() == 1);
    cfarBox_->setVisible(!cw && !cpiMode);
    noiseForm_->setRowVisible(mfWin_, !cw);
    cpiBox_->setVisible(cpiMode);
}

void MainWindow::updateValidation() {
    const radar::SimulationParams p = buildParams();
    const double b_bb_hz = bwBb_->value() * 1e6;
    const double fs_hz = fs_->value() * 1e6;
    const double f_if_hz = std::abs(fbb_->value()) * 1e6;
    (void)f_if_hz; // Marked for upcoming DSP heterodyne mixer integration

    if (fs_hz < b_bb_hz) {
        valid_ = false;
        validation_->setStyleSheet("color: #ff6b6b;");
        validation_->setText(QString("Nyquist Error: ADC sampling rate (%1 MHz) < Baseband Filter BW (%2 MHz). Aliasing will occur.")
                                 .arg(fs_->value(), 0, 'f', 2).arg(bwBb_->value(), 0, 'f', 2));
        info_->clear();
        generate_->setEnabled(false);
        return;
    }

    if (mode_->currentIndex() == 1) { updateCpiValidation(p); return; }
    const std::string err = radar::validate_all(p);
    if (!err.empty()) {
        valid_ = false;
        validation_->setStyleSheet("color: #ff6b6b;");
        validation_->setText(QString::fromStdString(err));
        info_->clear();
        generate_->setEnabled(false);
        return;
    }
    valid_ = true;
    validation_->setStyleSheet("color: #5fd38d;");
    validation_->setText("Parameters OK");

    const radar::ChannelReport rep = radar::analyze(p.wf, p.ch);
    const bool cw = (p.wf.type == radar::WaveformType::CW);
    const double wloss = cw ? 0.0 : radar::weighting_loss_db(p.mf_weighting, static_cast<std::size_t>(std::llround(p.wf.tau * p.wf.fs)));
    
    const double gateM = radar::kSpeedOfLight / (2.0 * fs_hz);
    const double rangeResM = radar::kSpeedOfLight / (2.0 * std::max(b_bb_hz, p.wf.bandwidth));

    QString s = QString("N = %1 ADC samples (f_s = %2 MHz) | Baseband BW B_bb = %3 MHz (f_IF = %4 MHz)")
                    .arg(radar::total_samples(p.wf))
                    .arg(fs_->value(), 0, 'f', 1)
                    .arg(bwBb_->value(), 0, 'f', 1)
                    .arg(fbb_->value(), 0, 'f', 1);
    s += QString("\nRange Gate Spacing (f_s) = %1 m | Physical Range Res (B_bb) = %2 m | Delay %3 µs")
             .arg(gateM, 0, 'f', 2)
             .arg(rangeResM, 0, 'f', 2)
             .arg(rep.td_s * 1e6, 0, 'f', 2);
    if (!cw) s += QString(" | R blind %1 km").arg(rep.r_blind_m / 1e3, 0, 'f', 2);
    if (rep.doppler_hz != 0.0) s += QString(" | Doppler %1 kHz").arg(rep.doppler_hz / 1e3, 0, 'f', 3);
    s += QString("\nInput SNR %1 dB").arg(rep.snr_in_db, 0, 'f', 1);
    if (!cw)
        s += QString(" | MF gain +%1 dB | expected output SNR %2 dB")
                 .arg(rep.mf_gain_db + wloss, 0, 'f', 1)
                 .arg(rep.snr_in_db + rep.mf_gain_db + wloss, 0, 'f', 1);
    info_->setText(s);
    generate_->setEnabled(!busy_);
}

void MainWindow::setBusy(bool busy) {
    busy_ = busy;
    generate_->setEnabled(!busy && valid_);
    cancel_->setEnabled(busy);
}

void MainWindow::onGenerate() {
    if (busy_ || !valid_) return;
    if (mode_->currentIndex() == 1) {
        CpiRequest cr;
        cr.sim = buildParams();
        cr.cpi = buildCpiParams();
        setBusy(true);
        validation_->setStyleSheet("");
        validation_->setText("Running CPI…");
        clock_.start();
        emit requestCpi(cr);
        return;
    }
    RunRequest rq;
    rq.sim = buildParams();
    rq.psd_window = (psdWindow_->currentIndex() == 0) ? radar::Window::Hann : radar::Window::Rectangular;
    setBusy(true);
    validation_->setStyleSheet("");
    validation_->setText("Running…");
    clock_.start();
    emit requestRun(rq);
}

void MainWindow::onFinished(PlotDataPtr data) {
    setBusy(false);
    const auto& r = *data->result;
    const auto& rep = r.report;
    QString s = QString("%1 ms | N = %2 | input SNR %3 dB")
                    .arg(clock_.elapsed())
                    .arg(r.rx.size())
                    .arg(rep.snr_in_db, 0, 'f', 1);
    if (!r.has_range_profile) {
        s += " | CW: no range profile (matched filter disabled)";
    } else {
        s += QString(" | MF gain +%1 dB, expected output SNR %2 dB | ")
                 .arg(rep.mf_gain_db, 0, 'f', 1)
                 .arg(rep.snr_in_db + rep.mf_gain_db, 0, 'f', 1);
        if (r.detection_range_m.empty()) {
            s += "no detection";
        } else {
            s += QString("%1 detection(s): ").arg(r.detection_range_m.size());
            for (std::size_t i = 0; i < r.detection_range_m.size() && i < 6; ++i)
                s += QString("%1 km  ").arg(r.detection_range_m[i] / 1e3, 0, 'f', 3);
        }
    }
    if (rep.in_blind_zone) s += " | FLAG: blind zone";
    if (rep.straddles_pri) s += " | FLAG: eclipsed";
    if (rep.second_time_around) s += " | FLAG: second-time-around";
    runSummary_ = s;
    status_->setText(s);

    data_ = std::move(data);
    plot_->setData(data_);
    updatePeakInfo();
    refreshIq();
    refreshPsd();
    resetPassbandView();
    updateExportOptions();
    updateValidation();
    stack_->setCurrentIndex(1);
}

void MainWindow::onFailed(const QString& message) {
    setBusy(false);
    validation_->setStyleSheet("color: #ff6b6b;");
    validation_->setText("Error: " + message);
}

void MainWindow::refreshIq() {
    if (!data_) return;
    const auto& wf = data_->result->params.wf;
    const bool rx = (iqSource_->currentIndex() == 0);
    const bool us = wf.pri < 2e-3;
    const double ys = rx ? 1.0 / data_->sigma : 1.0;
    std::vector<SeriesTrace> t;
    t.push_back({rx ? &data_->rx_i : &data_->tx_i, QColor(77, 179, 255, 150), ys, QStringLiteral("I")});
    t.push_back({rx ? &data_->rx_q : &data_->tx_q, QColor(255, 158, 38, 150), ys, QStringLiteral("Q")});
    iq_->setSeries(data_, std::move(t), 0.0, 1.0 / wf.fs, us ? 1e6 : 1e3,
                   us ? "Time (µs)" : "Time (ms)",
                   rx ? "amplitude (noise σ)" : "amplitude (ideal, unit)");
}

void MainWindow::refreshPsd() {
    if (!data_) return;
    const double fs = data_->result->params.wf.fs;
    const bool mhz = fs >= 2e5;
    std::vector<SeriesTrace> t;
    t.push_back({&data_->psd_ideal, QColor(77, 179, 255), 1.0, QStringLiteral("Tx ideal")});
    if (data_->psd_impaired.size() > 0)
        t.push_back({&data_->psd_impaired, QColor(255, 158, 38), 1.0, QStringLiteral("Tx with LO phase noise")});
    psd_->setSeries(data_, std::move(t), data_->psd_f_start, data_->psd_df, mhz ? 1e-6 : 1e-3,
                    mhz ? "Frequency offset (MHz)" : "Frequency offset (kHz)", "PSD (dB re ideal peak)");
}

void MainWindow::resetPassbandView() {
    if (!data_) return;
    const auto& res = *data_->result;
    const auto& wf = res.params.wf;
    const bool rx = (iqSource_->currentIndex() == 0);
    double centre_s;
    if (wf.type == radar::WaveformType::CW) centre_s = 0.5 * wf.pri;
    else if (rx) centre_s = std::fmod(res.report.td_s, wf.pri) + 0.5 * wf.tau;
    else centre_s = 0.5 * wf.tau;
    const double width_ns = std::clamp(20.0 / res.params.ch.fc * 1e9, 0.001, 1e6);
    {
        QSignalBlocker b1(pbCenter_), b2(pbWidth_);
        pbCenter_->setValue(centre_s * 1e6);
        pbWidth_->setValue(width_ns);
    }
    refreshPassband();
}

void MainWindow::refreshPassband() {
    if (!data_) return;
    const auto& res = *data_->result;
    const auto& wf = res.params.wf;
    const double fc = res.params.ch.fc;
    const bool rx = (iqSource_->currentIndex() == 0);
    const radar::CBuffer& x = rx ? res.rx : res.tx;
    constexpr std::size_t kPoints = 8192;
    constexpr double kMaxCycles = 1000.0;
    const double width_s = pbWidth_->value() * 1e-9;
    const double centre_s = pbCenter_->value() * 1e-6;
    pbCenter_->setSingleStep(std::max(pbWidth_->value() * 1e-3 * 0.5, 1e-6));
    const double cycles = width_s * fc;
    if (cycles > kMaxCycles) {
        pbInfo_->setStyleSheet("color: #ff6b6b;");
        pbInfo_->setText(QString("Window spans %1 carrier cycles (limit %2). Narrow the window.")
                             .arg(cycles, 0, 'f', 0).arg(kMaxCycles, 0, 'f', 0));
        pb_->clearSeries();
        return;
    }
    pbInfo_->setStyleSheet("");
    pbInfo_->setText(QString("%1 carrier cycles, %2 points per cycle")
                         .arg(cycles, 0, 'f', 1)
                         .arg(cycles > 0.0 ? static_cast<double>(kPoints) / cycles : 0.0, 0, 'f', 0));
    const double dt = width_s / static_cast<double>(kPoints);
    const double t0 = centre_s - 0.5 * width_s;
    radar::RBuffer rf;
    radar::reconstruct_passband(x, wf.fs, fc, t0, dt, kPoints, rf);

    struct Holder { radar::MinMaxPyramid pyr; };
    auto h = std::make_shared<Holder>();
    h->pyr = radar::MinMaxPyramid(std::move(rf));
    std::vector<SeriesTrace> t;
    t.push_back({&h->pyr, QColor(120, 220, 160), rx ? 1.0 / data_->sigma : 1.0, QStringLiteral("RF")});
    pb_->setSeries(h, std::move(t), -0.5 * width_s, dt, 1e9, "Time relative to centre (ns)",
                   rx ? "RF amplitude (noise σ)" : "RF amplitude (ideal)");
}

void MainWindow::updateExportOptions() {
    export_->setEnabled(data_ != nullptr && !busy_);
    if (!data_) return;
    const bool csvOk = data_->result->rx.size() <= radar::kCsvMaxSamples;
    if (auto* m = qobject_cast<QStandardItemModel*>(expFormat_->model())) {
        if (auto* it = m->item(3)) {
            it->setEnabled(csvOk);
            it->setText(csvOk ? "CSV (.csv)" : "CSV disabled (N > 100,000 samples)");
        }
    }
    if (!csvOk && expFormat_->currentIndex() == 3) expFormat_->setCurrentIndex(0);
}

void MainWindow::onExport() {
    if (!data_ || busy_) return;
    const int f = expFormat_->currentIndex();
    const bool tx = (f != 0) && expBuffer_->currentIndex() == 1;
    static const char* ext[] = {"h5", "bin", "bin", "csv"};
    static const char* filt[] = {"HDF5 (*.h5)", "Binary I/Q (*.bin)", "Binary I/Q (*.bin)", "CSV (*.csv)"};
    const QString base = (f == 0) ? "radar_sim" : (tx ? "tx_iq" : "rx_iq");
    QString path = QFileDialog::getSaveFileName(this, "Export", QDir::homePath() + "/" + base + "." + ext[f], filt[f]);
    if (path.isEmpty()) return;
    if (QFileInfo(path).suffix().isEmpty()) path += QString(".") + ext[f];
    ExportRequest rq;
    rq.result = data_->result;
    rq.path = path;
    rq.format = static_cast<ExportFormat>(f);
    rq.tx = tx;
    busy_ = true;
    export_->setEnabled(false);
    generate_->setEnabled(false);
    status_->setText("Exporting…");
    emit requestExport(rq);
}

void MainWindow::onExportDone(const QString& message) {
    busy_ = false;
    generate_->setEnabled(valid_);
    export_->setEnabled(data_ != nullptr);
    status_->setText(runSummary_ + "\nExported: " + message);
    cpiStatus_->setText("Exported: " + message);
    cpiExport_->setEnabled(cpi_ != nullptr);
}

void MainWindow::onExportFailed(const QString& message) {
    busy_ = false;
    generate_->setEnabled(valid_);
    export_->setEnabled(data_ != nullptr);
    status_->setText(runSummary_ + "\nExport failed: " + message);
    cpiStatus_->setText("Export failed: " + message);
    cpiExport_->setEnabled(cpi_ != nullptr);
}

void MainWindow::buildCpiControls() {
    mode_ = new QComboBox;
    mode_->addItems({"Single pulse (range profile)", "CPI (range-Doppler)"});
    cpiPulses_ = ispin(2, 4096, 64);
    cpiGate_ = ispin(8, 8192, 256);
    cpiDfft_ = ispin(0, 8192, 0);
    cpiDfft_->setSpecialValueText("auto");
    slowWin_ = new QComboBox;
    slowWin_->addItems({"Hann", "Hamming", "Rectangular"});
    cpiPfa_ = spin(1.0, 15.0, 6.0, 1, "");
    cpiPfa_->setPrefix("1e-");
    cGr_ = ispin(0, 500, 10);
    cGd_ = ispin(0, 100, 2);
    cTr_ = ispin(1, 1000, 20);
    cTd_ = ispin(0, 100, 4);
    auto pair = [](QWidget* a, QWidget* b) {
        auto* w = new QWidget;
        auto* l = new QHBoxLayout(w);
        l->setContentsMargins(0, 0, 0, 0);
        l->addWidget(a);
        l->addWidget(b);
        return w;
    };
    cpiBox_ = new QGroupBox("CPI and range-Doppler processing");
    auto* f = new QFormLayout(cpiBox_);
    f->addRow("Pulses per CPI", cpiPulses_);
    f->addRow("Range gate (bins kept per pulse)", cpiGate_);
    f->addRow("Slow-time window", slowWin_);
    f->addRow("Doppler FFT size", cpiDfft_);
    f->addRow("2D CFAR false-alarm probability", cpiPfa_);
    f->addRow("Guard cells: range / Doppler", pair(cGr_, cGd_));
    f->addRow("Training cells: range / Doppler", pair(cTr_, cTd_));
}

radar::CpiParams MainWindow::buildCpiParams() const {
    radar::CpiParams c;
    c.pulses = cpiPulses_->value();
    c.gate_len = cpiGate_->value();
    c.doppler_fft = cpiDfft_->value();
    c.window = (slowWin_->currentIndex() == 0)   ? radar::SlowWindow::Hann
               : (slowWin_->currentIndex() == 1) ? radar::SlowWindow::Hamming
                                                 : radar::SlowWindow::None;
    c.cfar.pfa = std::pow(10.0, -cpiPfa_->value());
    c.cfar.guard_range = cGr_->value();
    c.cfar.guard_doppler = cGd_->value();
    c.cfar.train_range = cTr_->value();
    c.cfar.train_doppler = cTd_->value();
    return c;
}

void MainWindow::updateCpiValidation(const radar::SimulationParams& p) {
    const radar::CpiParams c = buildCpiParams();
    const std::string err = radar::validate_cpi(p, c);
    if (!err.empty()) {
        valid_ = false;
        validation_->setStyleSheet("color: #ff6b6b;");
        validation_->setText(QString::fromStdString(err));
        info_->clear();
        generate_->setEnabled(false);
        return;
    }
    valid_ = true;
    validation_->setStyleSheet("color: #5fd38d;");
    validation_->setText("Parameters OK (CPI mode)");
    const radar::CpiReport rep = radar::analyze_cpi(p, c);
    const radar::ChannelReport ch = radar::analyze(p.wf, p.ch);
    const double wloss = radar::weighting_loss_db(p.mf_weighting, static_cast<std::size_t>(std::llround(p.wf.tau * p.wf.fs)));
    const std::size_t N = radar::total_samples(p.wf);
    const std::size_t G = std::min<std::size_t>(static_cast<std::size_t>(c.gate_len), N);
    QString s = QString("N = %1 samples/pulse | PRF %2 Hz | CPI %3 ms | map: %4 Doppler rows x %5 range bins")
                    .arg(N).arg(rep.prf_hz, 0, 'f', 1).arg(rep.cpi_s * 1e3, 0, 'f', 1)
                    .arg(radar::cpi_doppler_rows(c)).arg(G);
    s += QString("\nDoppler %1 Hz (apparent %2 Hz) | unambiguous velocity +/- %3 m/s | resolution %4 m/s (%5 Hz)")
             .arg(rep.doppler_true_hz, 0, 'f', 1).arg(rep.doppler_apparent_hz, 0, 'f', 1)
             .arg(rep.v_unambiguous_mps, 0, 'f', 2).arg(rep.velocity_resolution_mps, 0, 'f', 3)
             .arg(rep.doppler_resolution_hz, 0, 'f', 2);
    s += QString("\nInput SNR %1 dB | MF gain +%2 dB | slow-time gain +%3 dB | expected peak %4 dB above the RD floor")
             .arg(ch.snr_in_db, 0, 'f', 1).arg(ch.mf_gain_db + wloss, 0, 'f', 1).arg(rep.integration_gain_db, 0, 'f', 1)
             .arg(ch.snr_in_db + ch.mf_gain_db + wloss + rep.integration_gain_db, 0, 'f', 1);
    info_->setText(s);
    generate_->setEnabled(!busy_);
}

QWidget* MainWindow::buildCpiPage() {
    cpiBack_ = new QPushButton("← Setup");
    cpiStatus_ = new QLabel("No CPI data yet.");
    cpiStatus_->setWordWrap(true);
    rdRangeDb_ = spin(10.0, 120.0, 45.0, 0, " dB");
    cpiExport_ = new QPushButton("Export RD map (HDF5)…");
    rdInfo_ = new QLabel;
    rd_ = new RdPlot;
    rdRangeCut_ = new SeriesPlot;
    rdDopCut_ = new SeriesPlot;

    auto* top = new QHBoxLayout;
    top->addWidget(cpiBack_);
    top->addWidget(cpiStatus_, 1);
    top->addWidget(new QLabel("Colour range"));
    top->addWidget(rdRangeDb_);
    top->addWidget(cpiExport_);
    auto* cuts = new QHBoxLayout;
    cuts->addWidget(rdRangeCut_);
    cuts->addWidget(rdDopCut_);

    auto* page = new QWidget;
    auto* root = new QVBoxLayout(page);
    root->addLayout(top);
    root->addWidget(rd_, 3);
    root->addWidget(rdInfo_);
    root->addLayout(cuts, 2);

    connect(cpiBack_, &QPushButton::clicked, this, [this] { stack_->setCurrentIndex(0); });
    connect(cpiExport_, &QPushButton::clicked, this, &MainWindow::onCpiExport);
    connect(rdRangeDb_, &QDoubleSpinBox::valueChanged, this, [this](double v) { rd_->setColourRange(v); });
    connect(rd_, &RdPlot::cellSelected, this, &MainWindow::onCellSelected);
    connect(rd_, &RdPlot::cursorInfo, this, [this](const QString& s) { rdInfo_->setText(s); });
    return page;
}

void MainWindow::onCpiProgress(int done, int total) {
    validation_->setText(QString("Running CPI… pulse %1 of %2").arg(done).arg(total));
}

void MainWindow::onCpiFinished(radar::CpiResultPtr r) {
    setBusy(false);
    cpi_ = r;
    const auto& rep = r->report;
    QString s = QString("%1 ms | %2 pulses, %3 range bins, %4 Doppler rows | PRF %5 Hz | unambiguous +/- %6 m/s | resolution %7 m/s | ")
                    .arg(clock_.elapsed()).arg(r->cpi.pulses).arg(r->gate_len).arg(r->doppler_rows)
                    .arg(rep.prf_hz, 0, 'f', 1).arg(rep.v_unambiguous_mps, 0, 'f', 2).arg(rep.velocity_resolution_mps, 0, 'f', 3);
    if (r->detections.empty()) {
        s += "no detection";
    } else {
        s += QString("%1 detection(s): ").arg(r->detections.size());
        for (std::size_t i = 0; i < r->detections.size() && i < 3; ++i) {
            const auto& d = r->detections[i];
            s += QString("#%1 %2 m/s @ %3 km (%4 dB)  ").arg(i + 1).arg(d.velocity_refined_mps, 0, 'f', 3)
                     .arg(d.range_refined_m / 1e3, 0, 'f', 4).arg(d.snr_db, 0, 'f', 1);
        }
    }
    cpiStatus_->setText(s);
    rd_->setData(r, rdRangeDb_->value());

    int row = 0, col = 0;
    const std::size_t G = r->gate_len;
    if (!r->detections.empty()) {
        row = static_cast<int>(r->detections[0].doppler_bin);
        col = static_cast<int>(r->detections[0].gate_bin);
    } else {
        const std::size_t idx = static_cast<std::size_t>(std::max_element(r->rd_power.begin(), r->rd_power.end()) - r->rd_power.begin());
        row = static_cast<int>(idx / G);
        col = static_cast<int>(idx % G);
    }
    onCellSelected(row, col);
    updateValidation();
    stack_->setCurrentIndex(2);
}

void MainWindow::onCellSelected(int row, int col) {
    if (!cpi_) return;
    const auto& r = *cpi_;
    rd_->setSelected(row, col);
    updateCuts(row, col);
    const double lambda = radar::kSpeedOfLight / r.sim.ch.fc;
    const double hz = (row - static_cast<int>(r.doppler_rows / 2)) * r.report.prf_hz / static_cast<double>(r.doppler_rows);
    const std::size_t bin = (r.gate_start + static_cast<std::size_t>(col)) % r.n;
    const double level = 10.0 * std::log10(std::max(
        static_cast<double>(r.rd_power[static_cast<std::size_t>(row) * r.gate_len + static_cast<std::size_t>(col)]) / r.rd_noise_floor_w, 1e-14));
    rdInfo_->setText(QString("Selected: range %1 km | velocity %2 m/s (Doppler %3 Hz) | %4 dB re floor. Click any cell to move the cuts.")
                         .arg(radar::range_from_bin(static_cast<double>(bin), r.sim.wf.fs) / 1e3, 0, 'f', 3)
                         .arg(hz * lambda / 2.0, 0, 'f', 3).arg(hz, 0, 'f', 1).arg(level, 0, 'f', 1));
}

void MainWindow::updateCuts(int row, int col) {
    if (!cpi_) return;
    const auto& r = *cpi_;
    const std::size_t G = r.gate_len, L = r.doppler_rows;
    row = std::clamp(row, 0, static_cast<int>(L) - 1);
    col = std::clamp(col, 0, static_cast<int>(G) - 1);
    const double ref = r.rd_noise_floor_w;
    auto db = [&](float v) { return static_cast<float>(10.0 * std::log10(std::max(static_cast<double>(v) / ref, 1e-14))); };

    struct Holder { radar::MinMaxPyramid range, dop; };
    auto h = std::make_shared<Holder>();
    radar::RBuffer rc(G), dc(L);
    for (std::size_t g = 0; g < G; ++g) rc[g] = db(r.rd_power[static_cast<std::size_t>(row) * G + g]);
    for (std::size_t f = 0; f < L; ++f) dc[f] = db(r.rd_power[f * G + static_cast<std::size_t>(col)]);
    h->range = radar::MinMaxPyramid(std::move(rc));
    h->dop = radar::MinMaxPyramid(std::move(dc));

    const double binM = radar::kSpeedOfLight / (2.0 * r.sim.wf.fs);
    const double dv = r.report.prf_hz / static_cast<double>(L) * (radar::kSpeedOfLight / r.sim.ch.fc) / 2.0;
    std::vector<SeriesTrace> t1, t2;
    t1.push_back({&h->range, QColor(77, 179, 255), 1.0, QStringLiteral("range cut at the selected velocity")});
    t2.push_back({&h->dop, QColor(255, 158, 38), 1.0, QStringLiteral("velocity cut at the selected range")});
    rdRangeCut_->setSeries(h, std::move(t1), static_cast<double>(r.gate_start) * binM, binM, 1e-3,
                           "Range (km)", "dB re RD noise floor");
    rdDopCut_->setSeries(h, std::move(t2), -static_cast<double>(L / 2) * dv, dv, 1.0,
                         "Velocity (m/s)", "dB re RD noise floor");
}


// ---------------- Phase 2: impairments, extra targets, seed, presets, Monte Carlo ----------------

void MainWindow::buildExtraControls() {
    seed_ = ispin(0, 2147483647, 12345);
    adcBits_ = ispin(0, 16, 0);
    adcBits_->setSpecialValueText("off");
    adcFs_ = spin(0.0, 1000.0, 0.0, 2, " σ");
    adcFs_->setSpecialValueText("auto");
    iqGain_ = spin(-20.0, 20.0, 0.0, 2, " dB");
    iqPhase_ = spin(-45.0, 45.0, 0.0, 2, "°");
    dcI_ = spin(-100.0, 100.0, 0.0, 2, " σ");
    dcQ_ = spin(-100.0, 100.0, 0.0, 2, " σ");
    auto* impBox = new QGroupBox("Receiver impairments and noise seed");
    auto* f = new QFormLayout(impBox);
    f->addRow("Noise seed", seed_);
    f->addRow("ADC resolution", adcBits_);
    f->addRow("ADC full scale", adcFs_);
    f->addRow("I/Q gain imbalance", iqGain_);
    f->addRow("I/Q phase error", iqPhase_);
    f->addRow("DC offset I", dcI_);
    f->addRow("DC offset Q", dcQ_);

    auto* tgtBox = new QGroupBox("Additional targets: range / velocity / level re target 1");
    auto* g = new QGridLayout(tgtBox);
    for (int i = 0; i < 3; ++i) {
        tgtOn_[i] = new QCheckBox(QString("Target %1").arg(i + 2));
        tgtRange_[i] = spin(0.1, 5000.0, 160.0 + 10.0 * i, 3, " km");
        tgtVel_[i] = spin(-5000.0, 5000.0, 0.0, 2, " m/s");
        tgtRel_[i] = spin(-60.0, 40.0, -3.0, 1, " dB");
        g->addWidget(tgtOn_[i], i, 0);
        g->addWidget(tgtRange_[i], i, 1);
        g->addWidget(tgtVel_[i], i, 2);
        g->addWidget(tgtRel_[i], i, 3);
    }
    extrasRow_ = new QWidget;
    auto* h = new QHBoxLayout(extrasRow_);
    h->setContentsMargins(0, 0, 0, 0);
    h->addWidget(impBox, 1, Qt::AlignTop);
    h->addWidget(tgtBox, 1, Qt::AlignTop);

    savePreset_ = new QPushButton("Save preset…");
    loadPreset_ = new QPushButton("Load preset…");
    mcBtn_ = new QPushButton("Validate Pd (Monte Carlo)");
    for (auto* bt : {savePreset_, loadPreset_, mcBtn_}) bt->setMinimumHeight(40);
    connect(savePreset_, &QPushButton::clicked, this, &MainWindow::onSavePreset);
    connect(loadPreset_, &QPushButton::clicked, this, &MainWindow::onLoadPreset);
    connect(mcBtn_, &QPushButton::clicked, this, &MainWindow::onMonteCarlo);
}

void MainWindow::fillExtras(radar::SimulationParams& p) const {
    p.ch.seed = static_cast<std::uint64_t>(seed_->value());
    p.rx_imp.adc_bits = adcBits_->value();
    p.rx_imp.adc_fs_sigma = adcFs_->value();
    p.rx_imp.iq_gain_db = iqGain_->value();
    p.rx_imp.iq_phase_deg = iqPhase_->value();
    p.rx_imp.dc_i_sigma = dcI_->value();
    p.rx_imp.dc_q_sigma = dcQ_->value();
    p.extra_targets.clear();
    for (int i = 0; i < 3; ++i) {
        if (!tgtOn_[i]->isChecked()) continue;
        radar::TargetSpec t;
        t.range_m = tgtRange_[i]->value() * 1e3;
        t.velocity_mps = tgtVel_[i]->value();
        t.rel_db = tgtRel_[i]->value();
        p.extra_targets.push_back(t);
    }
}

radar::Preset MainWindow::currentPreset() const {
    radar::Preset pr;
    pr.sim = buildParams();
    pr.cpi = buildCpiParams();
    pr.cpi_mode = (mode_->currentIndex() == 1);
    pr.psd_window = psdWindow_->currentIndex();
    return pr;
}

void MainWindow::applyPreset(const radar::Preset& pr) {
    using radar::WaveformType;
    const radar::SimulationParams& s = pr.sim;
    int ti = 2;
    switch (s.wf.type) {
        case WaveformType::CW: ti = 0; break;
        case WaveformType::Pulse: ti = 1; break;
        case WaveformType::LFM: ti = 2; break;
        case WaveformType::Barker: ti = 3; break;
        case WaveformType::Frank: ti = 4; break;
    }
    type_->setCurrentIndex(ti);
    if (const int ci = code_->findData(s.wf.code_length); ci >= 0) code_->setCurrentIndex(ci);
    fc_->setValue(s.ch.fc / 1e9);
    fs_->setValue(s.wf.fs / 1e6);
    pri_->setValue(s.wf.pri * 1e3);
    tau_->setValue(s.wf.tau * 1e6);
    bw_->setValue(s.wf.bandwidth / 1e3);
    range_->setValue(s.ch.range_m / 1e3);
    dynamics_->setCurrentIndex(s.ch.velocity_mps != 0.0 ? 1 : 0);
    vel_->setValue(s.ch.velocity_mps);
    snrMode_->setCurrentIndex(s.ch.snr_mode == radar::SnrMode::Physical ? 1 : 0);
    snr_->setValue(s.ch.snr_db);
    pt_->setValue(s.ch.pt_w);
    gt_->setValue(s.ch.gt_db);
    gr_->setValue(s.ch.gr_db);
    rcs_->setValue(s.ch.rcs_m2);
    loss_->setValue(s.ch.loss_db);
    nf_->setValue(s.ch.noise_figure_db);
    noiseProfile_->setCurrentIndex(s.phase_noise_enabled ? 1 : 0);
    pnFloor_->setValue(s.pn.floor_dbc_hz);
    pnLeeson_->setValue(s.pn.leeson_corner_hz / 1e3);
    pnFlicker_->setValue(s.pn.flicker_corner_hz);
    mfWin_->setCurrentIndex(s.mf_weighting == radar::MfWeighting::Hamming ? 1 : s.mf_weighting == radar::MfWeighting::Hann ? 2 : 0);
    pfaExp_->setValue(-std::log10(s.cfar.pfa));
    guard_->setValue(s.cfar.guard_cells);
    train_->setValue(s.cfar.train_cells);
    psdWindow_->setCurrentIndex(pr.psd_window);
    seed_->setValue(static_cast<int>(s.ch.seed & 0x7fffffffULL));
    adcBits_->setValue(s.rx_imp.adc_bits);
    adcFs_->setValue(s.rx_imp.adc_fs_sigma);
    iqGain_->setValue(s.rx_imp.iq_gain_db);
    iqPhase_->setValue(s.rx_imp.iq_phase_deg);
    dcI_->setValue(s.rx_imp.dc_i_sigma);
    dcQ_->setValue(s.rx_imp.dc_q_sigma);
    for (int i = 0; i < 3; ++i) {
        const bool on = static_cast<std::size_t>(i) < s.extra_targets.size();
        tgtOn_[i]->setChecked(on);
        if (on) {
            tgtRange_[i]->setValue(s.extra_targets[static_cast<std::size_t>(i)].range_m / 1e3);
            tgtVel_[i]->setValue(s.extra_targets[static_cast<std::size_t>(i)].velocity_mps);
            tgtRel_[i]->setValue(s.extra_targets[static_cast<std::size_t>(i)].rel_db);
        }
    }
    cpiPulses_->setValue(pr.cpi.pulses);
    cpiGate_->setValue(pr.cpi.gate_len);
    cpiDfft_->setValue(pr.cpi.doppler_fft);
    slowWin_->setCurrentIndex(pr.cpi.window == radar::SlowWindow::Hann ? 0 : pr.cpi.window == radar::SlowWindow::Hamming ? 1 : 2);
    cpiPfa_->setValue(-std::log10(pr.cpi.cfar.pfa));
    cGr_->setValue(pr.cpi.cfar.guard_range);
    cGd_->setValue(pr.cpi.cfar.guard_doppler);
    cTr_->setValue(pr.cpi.cfar.train_range);
    cTd_->setValue(pr.cpi.cfar.train_doppler);
    mode_->setCurrentIndex(pr.cpi_mode ? 1 : 0);
    updateForm();
}

void MainWindow::onSavePreset() {
    QString path = QFileDialog::getSaveFileName(this, "Save preset", QDir::homePath() + "/radarsim_preset.json", "Preset (*.json)");
    if (path.isEmpty()) return;
    if (QFileInfo(path).suffix().isEmpty()) path += ".json";
    try {
        radar::save_preset(path.toStdString(), currentPreset());
        validation_->setStyleSheet("");
        validation_->setText("Preset saved: " + path);
    } catch (const std::exception& e) {
        validation_->setStyleSheet("color: #ff6b6b;");
        validation_->setText(QString("Could not save the preset: ") + e.what());
    }
}

void MainWindow::onLoadPreset() {
    const QString path = QFileDialog::getOpenFileName(this, "Load preset", QDir::homePath(), "Preset (*.json)");
    if (path.isEmpty()) return;
    try {
        applyPreset(radar::load_preset(path.toStdString()));
    } catch (const std::exception& e) {
        validation_->setStyleSheet("color: #ff6b6b;");
        validation_->setText(QString("Could not load the preset: ") + e.what());
    }
}

void MainWindow::onCpiExport() {
    if (!cpi_ || busy_) return;
    QString path = QFileDialog::getSaveFileName(this, "Export range-Doppler map", QDir::homePath() + "/radar_cpi.h5", "HDF5 (*.h5)");
    if (path.isEmpty()) return;
    if (QFileInfo(path).suffix().isEmpty()) path += ".h5";
    CpiExportRequest rq;
    rq.result = cpi_;
    rq.path = path;
    busy_ = true;
    cpiExport_->setEnabled(false);
    cpiStatus_->setText("Exporting…");
    emit requestCpiExport(rq);
}

QWidget* MainWindow::buildMcPage() {
    mcBack_ = new QPushButton("← Setup");
    mcStatus_ = new QLabel("Press 'Validate Pd (Monte Carlo)' on the setup page.");
    mcStatus_->setWordWrap(true);
    mcPlot_ = new SeriesPlot;
    mcText_ = new QPlainTextEdit;
    mcText_->setReadOnly(true);
    mcText_->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    mcText_->setMaximumHeight(280);
    auto* top = new QHBoxLayout;
    top->addWidget(mcBack_);
    top->addWidget(mcStatus_, 1);
    auto* page = new QWidget;
    auto* root = new QVBoxLayout(page);
    root->addLayout(top);
    root->addWidget(mcPlot_, 3);
    root->addWidget(mcText_, 2);
    connect(mcBack_, &QPushButton::clicked, this, [this] { stack_->setCurrentIndex(0); });
    return page;
}

void MainWindow::onMonteCarlo() {
    if (busy_) return;
    auto fail = [this](const QString& m) { validation_->setStyleSheet("color: #ff6b6b;"); validation_->setText(m); };
    if (mode_->currentIndex() == 1) { fail("Monte Carlo validation uses the single-pulse range profile: set Processing to 'Single pulse'."); return; }
    if (!valid_) return;   // the red validation message is already showing
    radar::SimulationParams p = buildParams();
    if (p.wf.type == radar::WaveformType::CW) { fail("Monte Carlo validation needs a pulsed waveform."); return; }
    const std::size_t N = radar::total_samples(p.wf);
    if (N > 2'000'000) { fail("The GUI limits Monte Carlo to 2M samples per PRI; use radar_mc for larger runs."); return; }
    McRequest rq;
    rq.mc.sim = p;
    rq.mc.sim.extra_targets.clear();          // validation uses the single-target model and an ideal receiver
    rq.mc.sim.rx_imp = radar::RxImpairments{};
    rq.mc.pfa = std::clamp(p.cfar.pfa, 1e-9, 0.1);
    rq.mc.trials = 300;
    rq.mc.noise_trials = static_cast<int>(std::clamp(std::ceil(2e7 / static_cast<double>(N)), 50.0, 2000.0));
    rq.mc.threads = 0;
    const double alb = radar::albersheim_snr_db(0.9, rq.mc.pfa);
    for (double v = std::floor(alb - 4.0); v <= std::ceil(alb + 4.0) + 1e-9; v += 1.0) rq.snr_out_db.push_back(v);
    setBusy(true);
    validation_->setStyleSheet("");
    validation_->setText("Running Monte Carlo (about 300 trials per SNR point)…");
    clock_.start();
    emit requestMc(rq);
}

void MainWindow::onMcFinished(McResultPtr r) {
    setBusy(false);
    const auto& c = r->curve;
    const auto& f = r->pfa;
    QString t;
    t += QString::asprintf("N = %zu samples (%zu active) | matched-filter gain +%.2f dB | design Pfa %.1e\n", r->n, r->n_pulse, r->mf_gain_db, f.pfa_design);
    t += QString::asprintf("Noise-only false-alarm rate over %zu cells:\n", f.cells);
    t += QString::asprintf("  fixed threshold : %.3e (%.2fx design)   <- validates the noise scaling\n", f.pfa_fixed, f.pfa_fixed / f.pfa_design);
    t += QString::asprintf("  CA-CFAR         : %.3e (%.2fx design), effective training cells %.1f\n\n", f.pfa_cfar, f.pfa_cfar / f.pfa_design, f.effective_cells);
    t += "  SNR_out(dB)  Pd theory  Pd fixed  Pd CA-CFAR\n";
    for (const auto& pt : c)
        t += QString::asprintf("  %9.1f    %.3f      %.3f     %.3f\n", pt.snr_out_db, pt.pd_theory, pt.pd_fixed, pt.pd_cfar);
    t += QString::asprintf("\nPd = 0.9 needs %.2f dB (Marcum-Q) | %.2f dB measured, fixed threshold | %.2f dB measured, CA-CFAR\n",
                           r->theory_snr90, r->measured_snr90_fixed, r->measured_snr90_cfar);
    mcText_->setPlainText(t);

    struct Holder { radar::MinMaxPyramid th, fx, cf; };
    auto h = std::make_shared<Holder>();
    radar::RBuffer a(c.size()), b(c.size()), d(c.size());
    for (std::size_t i = 0; i < c.size(); ++i) {
        a[i] = static_cast<float>(c[i].pd_theory);
        b[i] = static_cast<float>(c[i].pd_fixed);
        d[i] = static_cast<float>(c[i].pd_cfar);
    }
    h->th = radar::MinMaxPyramid(std::move(a));
    h->fx = radar::MinMaxPyramid(std::move(b));
    h->cf = radar::MinMaxPyramid(std::move(d));
    std::vector<SeriesTrace> tr;
    tr.push_back({&h->th, QColor(77, 179, 255), 1.0, QStringLiteral("theory (Marcum-Q)")});
    tr.push_back({&h->fx, QColor(255, 158, 38), 1.0, QStringLiteral("measured, fixed threshold")});
    tr.push_back({&h->cf, QColor(120, 220, 160), 1.0, QStringLiteral("measured, CA-CFAR")});
    const double step = c.size() > 1 ? c[1].snr_out_db - c[0].snr_out_db : 1.0;
    mcPlot_->setSeries(h, std::move(tr), c.front().snr_out_db, step, 1.0, "Matched-filter output SNR (dB)", "Probability of detection");
    mcStatus_->setText(QString::asprintf("Monte Carlo done in %lld ms | Pd = 0.9 at %.2f dB theory, %.2f dB measured (CA-CFAR)",
                                         static_cast<long long>(clock_.elapsed()), r->theory_snr90, r->measured_snr90_cfar));
    updateValidation();
    stack_->setCurrentIndex(3);
}


void MainWindow::updatePeakInfo() {
    if (!data_) { peakInfo_->clear(); return; }
    const auto& res = *data_->result;
    if (!res.has_range_profile) {
        peakInfo_->setText("CW carries no delay information: there is no range profile to detect a peak in.");
        return;
    }
    const auto& pk = data_->peaks;
    if (pk.empty()) { peakInfo_->setText("Peak detection: no peak found."); return; }
    const auto& t = pk.front();
    const double expect = res.report.snr_in_db + res.report.mf_gain_db;
    QString s = QString("<b>Target range (peak detection): %1 km</b> &nbsp;|&nbsp; bin %2 &nbsp;|&nbsp; %3 dB above the noise floor (expected %4 dB) &nbsp;|&nbsp; ")
                    .arg(t.range_m / 1e3, 0, 'f', 4).arg(t.bin_refined, 0, 'f', 2).arg(t.snr_db, 0, 'f', 1).arg(expect, 0, 'f', 1);
    s += t.above_threshold
             ? QString("<span style='color:#5fd38d'>above the %1 dB detection threshold</span>").arg(data_->peak_threshold_db, 0, 'f', 1)
             : QString("<span style='color:#ff9e4d'>below the %1 dB detection threshold: not a confident detection</span>").arg(data_->peak_threshold_db, 0, 'f', 1);
    const double truth = std::fmod(res.report.td_s, res.params.wf.pri) * radar::kSpeedOfLight / 2.0;
    if (res.params.extra_targets.empty())
        s += QString("<br>True range %1 km | error %2 m").arg(truth / 1e3, 0, 'f', 4).arg(t.range_m - truth, 0, 'f', 1);
    else
        s += QString("<br>Target 1 true range %1 km (this peak may belong to another target)").arg(truth / 1e3, 0, 'f', 4);
    if (pk.size() > 1) {
        s += "<br>Other peaks: ";
        for (std::size_t i = 1; i < pk.size() && i < 5; ++i)
            s += QString("%1 km (%2 dB)&nbsp;&nbsp;").arg(pk[i].range_m / 1e3, 0, 'f', 3).arg(pk[i].snr_db, 0, 'f', 1);
    }
    peakInfo_->setText(s);
}
