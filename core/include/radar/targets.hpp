#pragma once

namespace radar {

// Additional point target. rel_db shifts its received power relative to the primary target (manual SNR mode)
// or acts as an extra RCS offset (physical mode).
struct TargetSpec {
    double range_m = 0.0;
    double velocity_mps = 0.0;
    double rel_db = 0.0;
};

}  // namespace radar
