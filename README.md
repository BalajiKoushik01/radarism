# RadarSim

Baseband (I/Q) radar signal generation, propagation and DSP simulator with a Qt/OpenGL dashboard.
Built for Rocky Linux 9 (GCC 11.5, C++20, CMake >= 3.25, Qt 6.6).

## Features
- Waveforms: CW, unmodulated pulse, LFM, Barker, Frank (chip width validated against Fs)
- Tx phase noise (Leeson), fractional-delay target channel with carrier phase and Doppler
- SNR from a manual value or the radar equation; noise power per sample = k*T0*Fs*F
- FFT matched filter (optional Hamming/Hann taper), CA-CFAR (optionally sized for correlated noise), sub-bin range refinement
- CPI mode: per-pulse range compression, slow-time FFT, range-Doppler map, 2D CA-CFAR, ambiguity warnings
- Qt/OpenGL GUI: setup page with live validation; dashboard with matched filter, I/Q, PSD, lazy passband and range-Doppler heatmap
- Export: HDF5, binary float32/int16 with JSON sidecar, CSV (refused above 100,000 samples)
- Monte Carlo Pd/Pfa validation against Marcum-Q and Albersheim (radar_mc)

## Dependencies (Rocky 9)
    sudo dnf install -y epel-release
    sudo dnf config-manager --set-enabled crb
    sudo dnf install -y gcc-c++ pkgconf-pkg-config fftw-devel volk-devel gtest-devel \
                        hdf5-devel qt6-qtbase-devel mesa-libGL-devel

## Build and test
    cmake -S . -B build
    cmake --build build -j
    ctest --test-dir build --output-on-failure

## Run
    ./build/radarsim_gui           # GUI: Setup -> GENERATE -> dashboard
    ./build/radar_cli help         # single-pulse run, prints numbers, can export
    ./build/radar_cpi help         # CPI / range-Doppler run
    ./build/radar_mc help          # Monte Carlo Pd/Pfa validation
    ./build/gl_smoke 5             # OpenGL / GPU diagnostic

## Layout
    core/include/radar/   Qt-free DSP core (header-only)
    gui/                  Qt widgets, worker thread, OpenGL plots
    apps/                 command-line tools and the GL diagnostic
    tests/                GoogleTest suites (one per core module)
    packaging/            .desktop, icon, RPM spec

## Conventions
- Baseband samples are sqrt(watts): |x|^2 = watts. The Tx waveform has unit amplitude.
- Noise per sample = k*T0*Fs*F, split equally between I and Q.
- SNR in the UI is the per-sample INPUT SNR. Matched-filter gain is reported separately.
- Matched filter and CFAR are circular on purpose (periodic pulse train; second-time-around echoes wrap).
- Exports: S_tx is the clean waveform (what an AWG should play); S_rx is echo + AWGN.
  Binary files are normalized to full scale and carry a JSON sidecar with the scale factor.

## Validation (measured with radar_mc)
- Noise-only fixed-threshold false-alarm rate: 1.009e-3 against a design value of 1.000e-3.
- Pd against output SNR (Pfa 1e-3, LFM, 500-sample pulse): theory inside the 95% interval at 10 of 10 points.
  Pd = 0.9 needs 10.85 dB measured against 10.76 dB (Marcum-Q) and 10.72 dB (Albersheim).
- CA-CFAR on oversampled noise: with independent-cell sizing the false-alarm rate was 1.8x the design value
  at Pfa 1e-3 (64+64 training cells, B/Fs = 0.2) because neighbouring bins are correlated.
  The correlation-aware threshold (corr=on, default in the GUI and radar_cli) sizes alpha from the
  effective number of independent training cells and restores the design rate.

## Known simplifications
- Phase noise is applied to Tx only (pessimistic: a real monostatic radar partly cancels it with range correlation).
- Doppler is a constant frequency shift across one pulse plus a per-pulse phase in CPI mode (stop-and-hop).
- CPI mode does not model phase noise and uses an independence assumption in its 2D CFAR.
- Unweighted LFM has -13 dB range sidelobes, which can raise the CFAR threshold near strong targets.
- Range readout is quantized to c/(2Fs) before sub-bin refinement.

## Packaging
    git archive --format=tar.gz --prefix=radarsim-X.Y.Z/ -o ~/rpmbuild/SOURCES/radarsim-X.Y.Z.tar.gz vX.Y.Z
    rpmbuild -ba packaging/radarsim.spec

## Development
    git add -A && git commit -m "checkpoint before step N"   # before applying any patch
    git diff --stat                                          # what a patch changed
    git restore .                                            # undo a bad patch

## Phase 2 features
- Receiver impairments: ADC resolution/full scale/clipping, I/Q gain and phase imbalance, DC offset
  (GUI group "Receiver impairments and noise seed"; CLI: adc= adcfs= iqgain= iqphase= dci= dcq=)
- Additional targets in single-pulse and CPI modes (GUI "Additional targets"; CLI: t2=R_m,V_mps,dB, also t3, t4)
- CPI detections carry sub-cell range and Doppler estimates; the range-Doppler map exports to HDF5 from the CPI page
- Presets: every setup parameter saved/loaded as flat JSON; the noise seed is part of the preset
- "Validate Pd (Monte Carlo)": plots theory against measured Pd for the current waveform and detector
- Optional multi-threaded FFTW: set RADAR_FFT_THREADS=N (needs libfftw3f_threads); ./build/radar_bench measures it
