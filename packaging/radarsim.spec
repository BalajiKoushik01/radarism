%global debug_package %{nil}

Name:           radarsim
Version:        0.4.0
Release:        1%{?dist}
Summary:        Baseband radar signal generation, propagation and DSP simulator
License:        Proprietary
Source0:        %{name}-%{version}.tar.gz

BuildRequires:  gcc-c++
BuildRequires:  pkgconf-pkg-config
BuildRequires:  fftw-devel
BuildRequires:  volk-devel
BuildRequires:  hdf5-devel
BuildRequires:  gtest-devel
BuildRequires:  qt6-qtbase-devel
BuildRequires:  mesa-libGL-devel
BuildRequires:  desktop-file-utils

%description
Closed-loop-free baseband (I/Q) radar simulator: waveform generation (CW, pulse,
LFM, Barker, Frank), Leeson phase noise, fractional-delay target channel with
Doppler, radar-equation or manual SNR, FFT matched filter with CA-CFAR, an
OpenGL dashboard, and HDF5 / binary / CSV export for SDR and AWG hardware.

%prep
%autosetup

%build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix}
cmake --build build --parallel %{_smp_build_ncpus}

%check
ctest --test-dir build --output-on-failure

%install
DESTDIR=%{buildroot} cmake --install build
desktop-file-validate %{buildroot}%{_datadir}/applications/radarsim.desktop

%files
%{_bindir}/radar*
%{_datadir}/applications/radarsim.desktop
%{_datadir}/icons/hicolor/scalable/apps/radarsim.svg

%changelog
* Tue Oct 20 2026 Balaji Koushik - 0.4.0-1
- Receiver impairments, additional targets, CPI refinement and HDF5 export, presets
- Monte Carlo panel, optional multi-threaded FFTW, radar_bench

* Mon Oct 19 2026 Balaji Koushik - 0.3.0-1
- CPI / range-Doppler mode and heatmap, matched-filter weighting, sub-bin range refinement
- Instrument-style plots, Monte Carlo Pd/Pfa validation, correlation-aware CFAR, AWG/SDR export profile
- Installs the command-line tools (radar_cli, radar_cpi, radar_mc, radar_awg)

* Tue Oct 13 2026 Balaji Koushik - 0.1.0-1
- Initial package
