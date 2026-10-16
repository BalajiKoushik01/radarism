#include <gtest/gtest.h>
#include <cmath>
#include <numbers>
#include "radar/waveforms.hpp"

using namespace radar;

static WaveformParams base() {
    WaveformParams p;
    p.fs = 10e6; p.pri = 1e-3; p.tau = 100e-6;
    return p;
}

TEST(Waveform, PulseLengthAndZeroPadding) {
    auto p = base(); p.type = WaveformType::Pulse;
    auto w = generate(p);
    EXPECT_EQ(w.samples.size(), 10000u);
    EXPECT_EQ(w.n_pulse, 1000u);
    EXPECT_EQ(w.samples[0], cf32(1.0f, 0.0f));
    EXPECT_EQ(w.samples[999], cf32(1.0f, 0.0f));
    for (std::size_t n = 1000; n < w.samples.size(); ++n) ASSERT_EQ(w.samples[n], cf32(0.0f, 0.0f));
}

TEST(Waveform, CwFillsWholePri) {
    auto p = base(); p.type = WaveformType::CW;
    auto w = generate(p);
    EXPECT_EQ(w.n_pulse, w.samples.size());
}

TEST(Waveform, LfmSweepsMinusBToPlusB) {
    auto p = base(); p.type = WaveformType::LFM; p.bandwidth = 2e6;
    auto w = generate(p);
    const double k = p.bandwidth / p.tau;
    // Frequency measured between samples n and n+1 belongs to the midpoint n+0.5
    auto inst_freq = [&](std::size_t n) {
        return std::arg(w.samples[n + 1] * std::conj(w.samples[n])) * p.fs / (2.0 * std::numbers::pi);
    };
    auto expected = [&](std::size_t n) { return k * ((n + 0.5) / p.fs - 0.5 * p.tau); };
    EXPECT_NEAR(inst_freq(0), expected(0), 1e3);
    EXPECT_NEAR(inst_freq(500), expected(500), 1e3);
    EXPECT_NEAR(inst_freq(w.n_pulse - 2), expected(w.n_pulse - 2), 1e3);
    // Endpoints land within one sample-step of the +/- B/2 edges
    const double step = k / p.fs;  // Hz per sample
    EXPECT_NEAR(inst_freq(0), -p.bandwidth / 2, 2 * step);
    EXPECT_NEAR(inst_freq(w.n_pulse - 2), p.bandwidth / 2, 2 * step);
    for (std::size_t n = 0; n < w.n_pulse; ++n) ASSERT_NEAR(std::abs(w.samples[n]), 1.0f, 1e-5f);
}

TEST(Waveform, Barker13AutocorrelationSidelobes) {
    auto p = base(); p.type = WaveformType::Barker; p.code_length = 13; p.tau = 13e-6;  // 10 samples/chip
    auto w = generate(p);
    ASSERT_EQ(w.samples_per_chip, 10u);
    ASSERT_EQ(w.n_pulse, 130u);
    const std::size_t spc = w.samples_per_chip;
    for (std::size_t lag = 0; lag < w.n_pulse; ++lag) {
        double r = 0;
        for (std::size_t i = 0; i + lag < w.n_pulse; ++i) r += w.samples[i].real() * w.samples[i + lag].real();
                if (lag == 0) {
            EXPECT_NEAR(r, 130.0, 1e-3);
        } else if (lag >= spc) {
            EXPECT_LE(std::fabs(r), 10.0 + 1e-3) << "lag " << lag;  // sidelobe <= 1 chip
        }
        // lags 1..spc-1 are the main-lobe skirt (linear ramp down from 130), so no sidelobe check there
    }
    // Main-lobe skirt is a clean linear ramp: r(1) = 130 - 13
    double r1 = 0;
    for (std::size_t i = 0; i + 1 < w.n_pulse; ++i) r1 += w.samples[i].real() * w.samples[i + 1].real();
    EXPECT_NEAR(r1, 117.0, 1e-3);
}

TEST(Waveform, FrankIsConstantModulusWithZeroFirstRow) {
    auto p = base(); p.type = WaveformType::Frank; p.code_length = 16; p.tau = 16e-6;
    auto w = generate(p);
    ASSERT_EQ(w.n_pulse, 160u);
    for (std::size_t n = 0; n < w.n_pulse; ++n) ASSERT_NEAR(std::abs(w.samples[n]), 1.0f, 1e-5f);
    for (std::size_t n = 0; n < 40; ++n) EXPECT_NEAR(std::arg(w.samples[n]), 0.0f, 1e-5f);  // first M chips: phase 0
}

TEST(Validation, RejectsBadInputs) {
    auto p = base();
    p.type = WaveformType::Barker; p.code_length = 6;
    EXPECT_FALSE(validate(p).empty());                       // no Barker-6
    p.code_length = 13; p.tau = 12.5e-6;
    EXPECT_FALSE(validate(p).empty());                       // 125 samples / 13 chips: not an integer
    p.type = WaveformType::Frank; p.code_length = 10; p.tau = 100e-6;
    EXPECT_FALSE(validate(p).empty());                       // 10 is not a perfect square
    p.type = WaveformType::LFM; p.bandwidth = 9e6;
    EXPECT_FALSE(validate(p).empty());                       // Fs < 1.2 B
    p.type = WaveformType::Pulse; p.tau = 2e-3;
    EXPECT_FALSE(validate(p).empty());                       // tau > PRI
    EXPECT_THROW(generate(p), std::invalid_argument);
}

TEST(Validation, AcceptsGoodInputs) {
    auto p = base();
    p.type = WaveformType::Barker; p.code_length = 13; p.tau = 13e-6;
    EXPECT_TRUE(validate(p).empty());
}
