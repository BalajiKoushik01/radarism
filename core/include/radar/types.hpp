#pragma once
#include <complex>
#include <vector>
#include "radar/aligned_allocator.hpp"

namespace radar {

using cf32 = std::complex<float>;
using CBuffer = std::vector<cf32, AlignedAllocator<cf32, 64>>;  // I/Q samples
using RBuffer = std::vector<float, AlignedAllocator<float, 64>>; // real-valued data

}  // namespace radar
