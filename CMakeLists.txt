# SPDX-License-Identifier: GPL-3.0-or-later
cmake_minimum_required(VERSION 3.20)
project(transitscope VERSION 0.8.2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(TRANSIT_WITH_CFITSIO  "Build the CFITSIO image backend"     ON)
option(TRANSIT_WITH_XISFCONV "Build the libxisfconv image backend" OFF)
if(NOT TRANSIT_WITH_CFITSIO AND NOT TRANSIT_WITH_XISFCONV)
  message(FATAL_ERROR "Enable at least one image backend (TRANSIT_WITH_CFITSIO / TRANSIT_WITH_XISFCONV)")
endif()

find_package(PkgConfig REQUIRED)
pkg_check_modules(ERFA REQUIRED IMPORTED_TARGET erfa)

add_library(transit
  src/image.cpp
  src/stats.cpp
  src/calib.cpp
  src/photometry.cpp
  src/timing.cpp
  src/transit_model.cpp
  src/phot_table.cpp
  src/lsq.cpp
  src/ensemble.cpp
  src/mcmc.cpp
  src/fit.cpp
  src/limb_darkening.cpp
  src/reduce.cpp
  src/project.cpp
  src/io/reader_factory.cpp)
target_include_directories(transit PUBLIC include)
target_include_directories(transit SYSTEM PRIVATE third_party/tomlplusplus)
target_compile_definitions(transit PRIVATE TRANSIT_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
find_package(Threads REQUIRED)
target_link_libraries(transit PUBLIC PkgConfig::ERFA Threads::Threads)
target_compile_options(transit PRIVATE -Wall -Wextra -Wpedantic)

if(TRANSIT_WITH_CFITSIO)
  pkg_check_modules(CFITSIO REQUIRED IMPORTED_TARGET cfitsio)
  target_sources(transit PRIVATE src/io/cfitsio_reader.cpp)
  target_link_libraries(transit PUBLIC PkgConfig::CFITSIO)
  target_compile_definitions(transit PUBLIC TRANSIT_WITH_CFITSIO=1)
endif()

if(TRANSIT_WITH_XISFCONV)
  pkg_check_modules(XISFCONV REQUIRED IMPORTED_TARGET xisfconv)
  target_sources(transit PRIVATE src/io/xisfconv_reader.cpp)
  target_link_libraries(transit PUBLIC PkgConfig::XISFCONV)
  # A static libxisfconv (its default build; xisfconv.pc then sets -DXISFCONV_STATIC) needs the
  # Libs.private dependencies (zlib, zstd) on the link line as well.
  if("-DXISFCONV_STATIC" IN_LIST XISFCONV_CFLAGS_OTHER)
    message(STATUS "libxisfconv is static: linking ${XISFCONV_STATIC_LIBRARIES}")
    target_link_directories(transit PUBLIC ${XISFCONV_STATIC_LIBRARY_DIRS})
    target_link_libraries(transit PUBLIC ${XISFCONV_STATIC_LIBRARIES})
  endif()
  target_compile_definitions(transit PUBLIC TRANSIT_WITH_XISFCONV=1)
endif()

add_executable(transitscope tools/transitscope.cpp)
target_link_libraries(transitscope PRIVATE transit)

add_executable(fitscompare tools/fitscompare.cpp)
target_link_libraries(fitscompare PRIVATE transit)

add_executable(model_dump tests/unit/model_dump.cpp)
target_link_libraries(model_dump PRIVATE transit)

add_executable(bench_reduce tests/unit/bench_reduce.cpp)
target_link_libraries(bench_reduce PRIVATE transit)

add_executable(reduce_api_test tests/unit/reduce_api_test.cpp)
target_link_libraries(reduce_api_test PRIVATE transit)

add_executable(fit_api_test tests/unit/fit_api_test.cpp)
target_link_libraries(fit_api_test PRIVATE transit)

add_executable(project_test tests/unit/project_test.cpp)
target_link_libraries(project_test PRIVATE transit)

option(TRANSIT_WITH_GUI "Build the Qt 6 desktop application transitscope-gui" OFF)
if(TRANSIT_WITH_GUI)
  find_package(Qt6 6.2 REQUIRED COMPONENTS Widgets Concurrent)
  set(CMAKE_AUTOMOC ON)
  add_executable(transitscope-gui
    gui/main.cpp gui/MainWindow.cpp gui/SettingsEditor.cpp gui/ImageView.cpp gui/PlotWidget.cpp gui/ResultsView.cpp)
  target_include_directories(transitscope-gui SYSTEM PRIVATE third_party/tomlplusplus)
  target_compile_definitions(transitscope-gui PRIVATE TRANSITSCOPE_VERSION="${PROJECT_VERSION}")
  target_link_libraries(transitscope-gui PRIVATE transit Qt6::Widgets Qt6::Concurrent)
  target_compile_options(transitscope-gui PRIVATE -Wall -Wextra)
  install(TARGETS transitscope-gui RUNTIME DESTINATION bin)

  find_package(Qt6 COMPONENTS Test)
  if(Qt6Test_FOUND)
    add_executable(gui_test tests/gui/gui_test.cpp
      gui/MainWindow.cpp gui/SettingsEditor.cpp gui/ImageView.cpp gui/PlotWidget.cpp gui/ResultsView.cpp)
    target_include_directories(gui_test SYSTEM PRIVATE third_party/tomlplusplus)
    target_compile_definitions(gui_test PRIVATE TRANSITSCOPE_VERSION="${PROJECT_VERSION}")
    target_link_libraries(gui_test PRIVATE transit Qt6::Widgets Qt6::Concurrent Qt6::Test)
  endif()
endif()
