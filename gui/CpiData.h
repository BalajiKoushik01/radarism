#pragma once
#include <QMetaType>
#include "radar/cpi.hpp"

struct CpiRequest {
    radar::SimulationParams sim;
    radar::CpiParams cpi;
};

Q_DECLARE_METATYPE(CpiRequest)
Q_DECLARE_METATYPE(radar::CpiResultPtr)
