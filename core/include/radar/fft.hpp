#pragma once
#include <fftw3.h>
#include <climits>
#include <cstddef>
#include <cstdlib>
#include <mutex>
#include <stdexcept>
#include <utility>
#include "radar/types.hpp"

namespace radar {

// FFTW planning/destruction is NOT thread-safe; execution is. One global lock guards the former.
inline std::mutex& fftw_planner_mutex() {
    static std::mutex m;
    return m;
}

// Opt-in multi-threaded FFTs for large transforms: set RADAR_FFT_THREADS=N (N >= 2). Needs the fftw3f_threads
// library (RADAR_FFTW_THREADS is defined by CMake when it is found). Default: single-threaded, as before.
inline int fft_threads_requested() {
    static const int n = [] {
        const char* e = std::getenv("RADAR_FFT_THREADS");
        const int v = e ? std::atoi(e) : 1;
        return (v >= 1 && v <= 256) ? v : 1;
    }();
    return n;
}

// Call with the planner mutex held. Returns the number of threads the next plan will use.
inline int fft_set_threads(std::size_t n) {
#ifdef RADAR_FFTW_THREADS
    if (fft_threads_requested() <= 1) return 1;
    static const bool ready = (fftwf_init_threads() != 0);
    if (!ready) return 1;
    const int t = (n >= (static_cast<std::size_t>(1) << 17)) ? fft_threads_requested() : 1;
    fftwf_plan_with_nthreads(t);
    return t;
#else
    (void)n;
    return 1;
#endif
}

// Out-of-place, unnormalized 1-D complex FFT (like FFTW: forward then backward gives N * x).
// Buffers passed to execute() must be CBuffer-aligned and in != out.
class Fft {
public:
    Fft(std::size_t n, bool forward) : n_(n) {
        if (n == 0 || n > static_cast<std::size_t>(INT_MAX)) throw std::invalid_argument("Fft: bad size");
        CBuffer a(n), b(n);  // scratch buffers only used so FFTW can plan with correct alignment
        std::lock_guard<std::mutex> lk(fftw_planner_mutex());
        fft_set_threads(n);
        plan_ = fftwf_plan_dft_1d(static_cast<int>(n),
                                  reinterpret_cast<fftwf_complex*>(a.data()),
                                  reinterpret_cast<fftwf_complex*>(b.data()),
                                  forward ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
        if (!plan_) throw std::runtime_error("Fft: plan creation failed");
    }
    ~Fft() {
        if (plan_) {
            std::lock_guard<std::mutex> lk(fftw_planner_mutex());
            fftwf_destroy_plan(plan_);
        }
    }
    Fft(const Fft&) = delete;
    Fft& operator=(const Fft&) = delete;
    Fft(Fft&& o) noexcept : n_(o.n_), plan_(std::exchange(o.plan_, nullptr)) {}

    std::size_t size() const { return n_; }

    void execute(const cf32* in, cf32* out) const {
        fftwf_execute_dft(plan_, const_cast<fftwf_complex*>(reinterpret_cast<const fftwf_complex*>(in)),
                          reinterpret_cast<fftwf_complex*>(out));
    }

private:
    std::size_t n_;
    fftwf_plan plan_ = nullptr;
};

}  // namespace radar
