# Numeric input contract

RadarSim is a numerical simulation, not a physical RF output interlock.

Before any sample allocation or Doppler/propagation calculation, waveform and channel parameters must be finite. Invalid `NaN` and infinity values are rejected with human-readable validation messages, rather than reaching rounding, integer casts, exponential power calculations or FFT allocation.

The simulator still has physical-model simplifications and does not prove performance or equipment safety. Core CTest validates selected signal cases; separate instrument-specific calibration and hardware safeguards are mandatory for any connection to real RF hardware.

Regression tests in `tests/test_waveforms.cpp` and `tests/test_channel.cpp` cover non-finite inputs. Run the headless build as documented in `README.md`.
