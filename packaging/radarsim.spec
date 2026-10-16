%global debug_package %{nil}

Name:           radarsim
Version:        0.1.0
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
%{_bindir}/radarsim_gui
%{_bindir}/radar_cli
%{_datadir}/applications/radarsim.desktop
%{_datadir}/icons/hicolor/scalable/apps/radarsim.svg

%changelog
* Tue Oct 13 2026 Balaji Koushik - 0.1.0-1
- Initial package
