# RadarSim development notes

## Layout
    core/include/radar/   Qt-free DSP core (header-only); one test suite per module in tests/
    gui/                  Qt widgets, SimWorker (QThread), OpenGL plots (PlotWidget, SeriesPlot, RdPlot)
    apps/                 radar_cli, radar_cpi, radar_mc, radar_awg, radar_bench, gl_smoke
    packaging/            RPM spec, .desktop, icon

## Signal chain
Tx waveform -> Leeson phase noise (on a copy; the clean Tx is what an AWG plays) -> target channel
(fractional delay by windowed sinc, carrier phase, Doppler, radar equation or manual SNR) -> additional targets ->
AWGN (k*T0*Fs*F per sample) -> receiver impairments (I/Q imbalance, DC offset, ADC) -> FFT matched filter with the
clean (optionally tapered) replica -> peak detection and CA-CFAR. CPI mode: per-pulse compression, slow-time FFT,
range-Doppler map, 2D CA-CFAR.

## Conventions
- Samples are sqrt(watts): |x|^2 = watts. Tx has unit amplitude. SNR in the UI is the per-sample INPUT SNR.
- Matched filter and CFAR are circular on purpose (periodic pulse train, second-time-around echoes wrap).
- The CFAR threshold is sized for correlated oversampled noise (effective independent cells from the replica
  autocorrelation); the plain independent-cell formula gives 1.8x (Pfa 1e-3) to 8x (Pfa 1e-6) too many false alarms.
- Peak detection (peak_detect.hpp): local maxima, non-maximum suppression over one compressed-pulse width,
  log-parabola sub-bin refinement; the strongest peak is always reported, flagged when below the Pfa threshold.

## Validation (radar_mc)
Fixed-threshold Pfa ratio 1.008; Pd theory inside the 95% interval at 10 of 10 points; Pd = 0.9 needs 10.85 dB
(Marcum-Q 10.76 dB, Albersheim 10.72 dB); correlation-aware CA-CFAR Pfa ratio 0.97 at 1e-3 and 1e-6.

## Known limitations
- Phase noise is Tx-only (pessimistic) and not modelled in CPI mode.
- CPI 2D CFAR assumes independent cells (the range direction is correlated).
- Unweighted LFM has -13 dB range sidelobes; they can raise the CFAR threshold and appear as extra peaks.
- Doppler within one pulse is a constant shift; CPI mode is stop-and-hop.

## Workflow that worked
- git commit a checkpoint before every patch; patch scripts apply all edits or none (text anchors).
- Staged scripts with automatic rollback and a summary file for unattended runs (see ~/radarsim_stages3).
- Long pasted heredocs garble the terminal echo but still write correctly; a clean compile is the check.
