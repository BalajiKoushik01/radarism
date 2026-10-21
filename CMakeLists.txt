cmake_minimum_required(VERSION 3.25)
project(RadarSim VERSION 0.4.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(RADARSIM_BUILD_TESTS "Build the unit tests" ON)
option(RADARSIM_BUILD_GUI "Build the Qt GUI" ON)

include(GNUInstallDirs)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3F REQUIRED IMPORTED_TARGET fftw3f)
find_package(Volk QUIET)
if(TARGET Volk::volk)
  set(RADAR_VOLK_TARGET Volk::volk)
else()
  pkg_check_modules(VOLK REQUIRED IMPORTED_TARGET volk)
  set(RADAR_VOLK_TARGET PkgConfig::VOLK)
endif()
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(Threads REQUIRED)

# Header-only, Qt-free core
add_library(radar_core INTERFACE)
target_include_directories(radar_core INTERFACE ${CMAKE_SOURCE_DIR}/core/include)
target_link_libraries(radar_core INTERFACE PkgConfig::FFTW3F ${RADAR_VOLK_TARGET} Threads::Threads)
find_library(FFTW3F_THREADS_LIB NAMES fftw3f_threads)
if(FFTW3F_THREADS_LIB)
  target_compile_definitions(radar_core INTERFACE RADAR_FFTW_THREADS=1)
  target_link_libraries(radar_core INTERFACE ${FFTW3F_THREADS_LIB})
  message(STATUS "FFTW threads: enabled (${FFTW3F_THREADS_LIB})")
else()
  message(STATUS "FFTW threads: library not found, single-threaded FFTs")
endif()
target_compile_options(radar_core INTERFACE -Wall -Wextra -Wpedantic)

# Core + exporters (HDF5)
add_library(radar_export INTERFACE)
target_link_libraries(radar_export INTERFACE radar_core HDF5::HDF5)

add_executable(radar_cli apps/radar_cli.cpp)
target_link_libraries(radar_cli PRIVATE radar_export)
install(TARGETS radar_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(radar_cpi apps/radar_cpi.cpp)   # CPI / range-Doppler command-line tool (not installed)
target_link_libraries(radar_cpi PRIVATE radar_core)

add_executable(radar_bench apps/radar_bench.cpp)
target_link_libraries(radar_bench PRIVATE radar_export)

add_executable(radar_awg apps/radar_awg.cpp)   # AWG / SDR playback-file tool
target_link_libraries(radar_awg PRIVATE radar_core)

add_executable(radar_mc apps/radar_mc.cpp)   # Monte Carlo Pd/Pfa validation tool (not installed)
target_link_libraries(radar_mc PRIVATE radar_core)

if(RADARSIM_BUILD_TESTS)
  find_package(GTest REQUIRED)
  enable_testing()
  foreach(t smoke waveforms spectral channel detection pipeline pyramid passband refine cpi montecarlo cfar_neff awg impairments targets cpi_refine preset fft_threads peak_detect)
    add_executable(test_${t} tests/test_${t}.cpp)
    target_link_libraries(test_${t} PRIVATE radar_export GTest::gtest_main)
    add_test(NAME ${t} COMMAND test_${t})
  endforeach()
endif()

if(RADARSIM_BUILD_GUI)
  find_package(Qt6 6.2 REQUIRED COMPONENTS Widgets OpenGLWidgets)
  set(CMAKE_AUTOMOC ON)
  add_executable(radarsim_gui
    gui/main.cpp
    gui/MainWindow.cpp gui/MainWindow.h
    gui/PlotWidget.cpp gui/PlotWidget.h
    gui/SeriesPlot.cpp gui/SeriesPlot.h
    gui/SimWorker.cpp gui/SimWorker.h
    gui/RdPlot.cpp gui/RdPlot.h
    gui/PlotData.h gui/CpiData.h gui/McData.h)
  target_include_directories(radarsim_gui PRIVATE gui)
  target_compile_definitions(radarsim_gui PRIVATE RADARSIM_VERSION="${PROJECT_VERSION}")
  target_link_libraries(radarsim_gui PRIVATE radar_export Qt6::Widgets Qt6::OpenGLWidgets)

  add_executable(gl_smoke apps/gl_smoke.cpp)   # GPU / OpenGL diagnostic, not installed
  target_link_libraries(gl_smoke PRIVATE Qt6::Widgets Qt6::OpenGLWidgets)

  install(TARGETS radarsim_gui RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES packaging/radarsim.desktop DESTINATION ${CMAKE_INSTALL_DATADIR}/applications)
  install(FILES packaging/radarsim.svg DESTINATION ${CMAKE_INSTALL_DATADIR}/icons/hicolor/scalable/apps)
endif()

foreach(tool radar_cpi radar_mc radar_awg)
  if(TARGET ${tool})
    install(TARGETS ${tool} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()
endforeach()
