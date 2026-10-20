#pragma once
#include <QMetaType>
#include <QString>
#include <memory>
#include <vector>
#include "radar/cpi.hpp"
#include "radar/montecarlo.hpp"

struct McRequest {
    radar::McParams mc;
    std::vector<double> snr_out_db;
};
struct CpiExportRequest {
    radar::CpiResultPtr result;
    QString path;
};
using McResultPtr = std::shared_ptr<const radar::McResult>;

Q_DECLARE_METATYPE(McRequest)
Q_DECLARE_METATYPE(CpiExportRequest)
Q_DECLARE_METATYPE(McResultPtr)
