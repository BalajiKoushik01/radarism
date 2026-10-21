#pragma once
#include <QElapsedTimer>
#include <QMainWindow>
#include "CpiData.h"
#include "McData.h"
#include "radar/preset.hpp"
#include "PlotData.h"

class QCheckBox;
class QComboBox;
class QPlainTextEdit;
class QDoubleSpinBox;
class QFormLayout;
class QGroupBox;
class QLabel;
class QPushButton;
class QSpinBox;
class QStackedWidget;
class QThread;
class PlotWidget;
class RdPlot;
class SeriesPlot;
class SimWorker;

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget* parent = nullptr);
    ~MainWindow() override;

signals:
    void requestRun(const RunRequest& rq);
    void requestExport(const ExportRequest& rq);
    void requestCpi(const CpiRequest& rq);
    void requestMc(const McRequest& rq);
    void requestCpiExport(const CpiExportRequest& rq);

private slots:
    void onGenerate();
    void onFinished(PlotDataPtr data);
    void onFailed(const QString& message);
    void onExport();
    void onExportDone(const QString& message);
    void onExportFailed(const QString& message);
    void onCpiProgress(int done, int total);
    void onCpiFinished(radar::CpiResultPtr r);
    void onCellSelected(int row, int col);
    void onSavePreset();
    void onLoadPreset();
    void onCpiExport();
    void onMonteCarlo();
    void onMcFinished(McResultPtr r);

private:
    QWidget* buildSetupPage();
    QWidget* buildDashboardPage();
    QWidget* buildCpiPage();
    void buildCpiControls();
    radar::SimulationParams buildParams() const;
    radar::CpiParams buildCpiParams() const;
    void populateCodes();
    void updateForm();
    void updateVisibility();
    void updateValidation();
    void updateCpiValidation(const radar::SimulationParams& p);
    void setBusy(bool busy);
    void refreshIq();
    void refreshPsd();
    void refreshPassband();
    void resetPassbandView();
    void updateExportOptions();
    void updateCuts(int row, int col);
    void updatePeakInfo();
    void buildExtraControls();
    void fillExtras(radar::SimulationParams& p) const;
    radar::Preset currentPreset() const;
    void applyPreset(const radar::Preset& pr);
    QWidget* buildMcPage();

    QStackedWidget* stack_ = nullptr;
    QComboBox *type_ = nullptr, *code_ = nullptr, *dynamics_ = nullptr, *snrMode_ = nullptr;
    QComboBox *noiseProfile_ = nullptr, *psdWindow_ = nullptr, *iqSource_ = nullptr, *mfWin_ = nullptr;
    QDoubleSpinBox *fc_ = nullptr, *fbb_ = nullptr, *bwBb_ = nullptr, *fs_ = nullptr;
    QDoubleSpinBox *pri_ = nullptr, *tau_ = nullptr, *bw_ = nullptr;
    QDoubleSpinBox *range_ = nullptr, *vel_ = nullptr, *snr_ = nullptr, *pt_ = nullptr, *gt_ = nullptr;
    QDoubleSpinBox *gr_ = nullptr, *rcs_ = nullptr, *loss_ = nullptr, *nf_ = nullptr;
    QDoubleSpinBox *pnFloor_ = nullptr, *pnLeeson_ = nullptr, *pnFlicker_ = nullptr, *pfaExp_ = nullptr;
    QSpinBox *guard_ = nullptr, *train_ = nullptr;
    QFormLayout *signalForm_ = nullptr, *targetForm_ = nullptr, *noiseForm_ = nullptr;
    QGroupBox* cfarBox_ = nullptr;
    QLabel *validation_ = nullptr, *info_ = nullptr, *status_ = nullptr;
    QPushButton *generate_ = nullptr, *cancel_ = nullptr, *back_ = nullptr;
    PlotWidget* plot_ = nullptr;
    SeriesPlot *iq_ = nullptr, *psd_ = nullptr, *pb_ = nullptr;
    QDoubleSpinBox *pbCenter_ = nullptr, *pbWidth_ = nullptr;
    QLabel* pbInfo_ = nullptr;
    QComboBox *expBuffer_ = nullptr, *expFormat_ = nullptr;
    QPushButton* export_ = nullptr;

    QComboBox *mode_ = nullptr, *slowWin_ = nullptr;
    QSpinBox *cpiPulses_ = nullptr, *cpiGate_ = nullptr, *cpiDfft_ = nullptr;
    QSpinBox *cGr_ = nullptr, *cGd_ = nullptr, *cTr_ = nullptr, *cTd_ = nullptr;
    QDoubleSpinBox *cpiPfa_ = nullptr, *rdRangeDb_ = nullptr;
    QGroupBox* cpiBox_ = nullptr;
    QLabel *cpiStatus_ = nullptr, *rdInfo_ = nullptr;
    QLabel* peakInfo_ = nullptr;
    QPushButton* cpiBack_ = nullptr;
    RdPlot* rd_ = nullptr;
    SeriesPlot *rdRangeCut_ = nullptr, *rdDopCut_ = nullptr;
    radar::CpiResultPtr cpi_;
    QSpinBox *seed_ = nullptr, *adcBits_ = nullptr;
    QDoubleSpinBox *adcFs_ = nullptr, *iqGain_ = nullptr, *iqPhase_ = nullptr, *dcI_ = nullptr, *dcQ_ = nullptr;
    QCheckBox* tgtOn_[3] = {};
    QDoubleSpinBox *tgtRange_[3] = {}, *tgtVel_[3] = {}, *tgtRel_[3] = {};
    QWidget* extrasRow_ = nullptr;
    QPushButton *savePreset_ = nullptr, *loadPreset_ = nullptr, *mcBtn_ = nullptr, *mcBack_ = nullptr, *cpiExport_ = nullptr;
    QLabel* mcStatus_ = nullptr;
    SeriesPlot* mcPlot_ = nullptr;
    QPlainTextEdit* mcText_ = nullptr;

    QThread* thread_ = nullptr;
    SimWorker* worker_ = nullptr;
    QElapsedTimer clock_;
    PlotDataPtr data_;
    bool busy_ = false, valid_ = false;
    QString runSummary_;
};
