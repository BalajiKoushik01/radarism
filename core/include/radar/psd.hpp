#pragma once
#include <cmath>
#include <stdexcept>
#include <volk/volk.h>
#include "radar/fft.hpp"
#include "radar/types.hpp"

namespace radar {

enum class Window { Rectangular, Hann };

inline RBuffer make_window(std::size_t n, Window w) {
    RBuffer win(n, 1.0f);
    if (w == Window::Hann)
        for (std::size_t i = 0; i < n; ++i)
            win[i] = static_cast<float>(0.5 - 0.5 * std::cos(2.0 * 3.14159265358979323846 * i / n));
    return win;
}

// Two-sided PSD, linear units (power/Hz), fftshifted: bin i has frequency f_start + i*df.
// Normalization: PSD = |FFT(x*w)|^2 / (Fs * sum(w^2)), so sum(PSD)*df = signal power.
struct Psd {
    RBuffer density;
    double df = 0.0;
    double f_start = 0.0;
};

namespace detail {
inline Psd shift_and_scale(const CBuffer& X, double fs, double sw2) {
    const std::size_t n = X.size();
    Psd out;
    out.density.resize(n);
    out.df = fs / static_cast<double>(n);
    out.f_start = -static_cast<double>(n / 2) * out.df;
    const double scale = 1.0 / (fs * sw2);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t k = (i + n - n / 2) % n;
        out.density[i] = static_cast<float>(static_cast<double>(std::norm(X[k])) * scale);
    }
    return out;
}
}  // namespace detail

inline Psd compute_psd(const CBuffer& x, double fs, Window w = Window::Hann) {
    const std::size_t n = x.size();
    if (n < 2 || !(fs > 0.0)) throw std::invalid_argument("compute_psd: bad arguments");
    const RBuffer win = make_window(n, w);
    CBuffer xw(n), X(n);
    volk_32fc_32f_multiply_32fc(xw.data(), x.data(), win.data(), static_cast<unsigned int>(n));
    double sw2 = 0.0;
    for (float v : win) sw2 += static_cast<double>(v) * v;
    Fft fft(n, true);
    fft.execute(xw.data(), X.data());
    return detail::shift_and_scale(X, fs, sw2);
}

// PSD of the ACTIVE part of a pulse-train buffer: the first n_active samples are windowed, the rest of
// the buffer is zero, and the FFT length stays x.size() (fine frequency grid). Use this for waveform
// spectra: a window across the whole PRI would multiply the pulse by values near zero.
inline Psd compute_pulse_psd(const CBuffer& x, std::size_t n_active, double fs, Window w = Window::Rectangular) {
    const std::size_t n = x.size();
    if (n < 2 || n_active < 2 || n_active > n || !(fs > 0.0))
        throw std::invalid_argument("compute_pulse_psd: bad arguments");
    const RBuffer win = make_window(n_active, w);
    CBuffer xw(n, cf32{0.0f, 0.0f}), X(n);
    volk_32fc_32f_multiply_32fc(xw.data(), x.data(), win.data(), static_cast<unsigned int>(n_active));
    double sw2 = 0.0;
    for (float v : win) sw2 += static_cast<double>(v) * v;
    Fft fft(n, true);
    fft.execute(xw.data(), X.data());
    return detail::shift_and_scale(X, fs, sw2);
}

}  // namespace radar
