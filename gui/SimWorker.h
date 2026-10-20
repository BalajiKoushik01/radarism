#pragma once
#include <QObject>
#include <QString>
#include <atomic>
#include "CpiData.h"
#include "McData.h"
#include "PlotData.h"

class SimWorker : public QObject {
    Q_OBJECT
public:
    explicit SimWorker(QObject* parent = nullptr) : QObject(parent) {}
    void cancel() { cancel_.store(true); }   // thread-safe, callable from the UI thread

public slots:
    void run(const RunRequest& rq);
    void exportData(const ExportRequest& rq);
    void runCpi(const CpiRequest& rq);
    void runMc(const McRequest& rq);
    void exportCpi(const CpiExportRequest& rq);

signals:
    void finished(PlotDataPtr data);
    void failed(const QString& message);
    void exportDone(const QString& message);
    void exportFailed(const QString& message);
    void cpiFinished(radar::CpiResultPtr data);
    void cpiProgress(int done, int total);
    void mcFinished(McResultPtr r);

private:
    std::atomic<bool> cancel_{false};
};
