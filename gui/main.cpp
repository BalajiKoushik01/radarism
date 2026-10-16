#include <QApplication>
#include <QMetaType>
#include <QSurfaceFormat>
#include "MainWindow.h"
#include "CpiData.h"
#include "PlotData.h"

int main(int argc, char** argv) {
    QSurfaceFormat fmt;
    fmt.setVersion(3, 3);
    fmt.setProfile(QSurfaceFormat::CoreProfile);
    fmt.setSwapInterval(1);
    QSurfaceFormat::setDefaultFormat(fmt);

    QApplication app(argc, argv);
    qRegisterMetaType<radar::SimulationParams>("radar::SimulationParams");
    qRegisterMetaType<PlotDataPtr>("PlotDataPtr");
    qRegisterMetaType<RunRequest>("RunRequest");
    qRegisterMetaType<ExportRequest>("ExportRequest");
    qRegisterMetaType<CpiRequest>("CpiRequest");
    qRegisterMetaType<radar::CpiResultPtr>("radar::CpiResultPtr");

    MainWindow w;
    w.resize(1400, 1050);
    w.show();
    return app.exec();
}
