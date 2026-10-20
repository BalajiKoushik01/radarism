#pragma once
#include <cstdio>
#include <stdexcept>
#include <string>
#include "radar/simulation.hpp"

// Options shared by the command-line tools:
//   adc=BITS adcfs=SIGMA iqgain=DB iqphase=DEG dci=SIGMA dcq=SIGMA      receiver impairments
//   t2=RANGE_M,VEL_MPS,REL_DB   (also t3, t4)                           additional targets
template <class ArgsT>
inline void parse_extras(ArgsT& a, radar::SimulationParams& p) {
    p.rx_imp.adc_bits = static_cast<int>(a.num("adc", 0));
    p.rx_imp.adc_fs_sigma = a.num("adcfs", 0.0);
    p.rx_imp.iq_gain_db = a.num("iqgain", 0.0);
    p.rx_imp.iq_phase_deg = a.num("iqphase", 0.0);
    p.rx_imp.dc_i_sigma = a.num("dci", 0.0);
    p.rx_imp.dc_q_sigma = a.num("dcq", 0.0);
    for (const char* key : {"t2", "t3", "t4"}) {
        const std::string s = a.str(key, "");
        if (s.empty()) continue;
        double r = 0.0, v = 0.0, d = 0.0;
        if (std::sscanf(s.c_str(), "%lf,%lf,%lf", &r, &v, &d) != 3)
            throw std::invalid_argument(std::string(key) + " must be RANGE_M,VEL_MPS,REL_DB");
        p.extra_targets.push_back({r, v, d});
    }
}
