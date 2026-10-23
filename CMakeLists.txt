# SPDX-License-Identifier: GPL-3.0-or-later
cmake_minimum_required(VERSION 3.20)
project(transitscope VERSION 0.11.0 LANGUAGES CXX)

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

# ERFA: the system library if pkg-config finds it, otherwise built from the liberfa sources
# (fetched at configure time; needed on Windows/MSYS2, which has no ERFA package)
option(TRANSIT_ERFA_FROM_SOURCE "Always build ERFA from source instead of using the system library" OFF)
if(NOT TRANSIT_ERFA_FROM_SOURCE)
  pkg_check_modules(ERFA IMPORTED_TARGET erfa)
endif()
if(TARGET PkgConfig::ERFA)
  set(TRANSIT_ERFA_TARGET PkgConfig::ERFA)
else()
  message(STATUS "ERFA: building liberfa 2.0.1 from source")
  enable_language(C)
  include(FetchContent)
  FetchContent_Declare(liberfa
    GIT_REPOSITORY https://github.com/liberfa/erfa.git
    GIT_TAG 9915ba38c9365f8b0738269b8c2ac1fdd5f8dee3   # v2.0.1
    SOURCE_SUBDIR no-cmake)   # liberfa has no CMakeLists.txt: only download it, build it below
  FetchContent_MakeAvailable(liberfa)
  file(GLOB ERFA_SOURCES ${liberfa_SOURCE_DIR}/src/*.c)
  list(FILTER ERFA_SOURCES EXCLUDE REGEX "/t_erfa_c")
  add_library(erfa_bundled STATIC ${ERFA_SOURCES})
  target_include_directories(erfa_bundled PUBLIC ${liberfa_SOURCE_DIR}/src)
  target_compile_definitions(erfa_bundled PRIVATE PACKAGE_VERSION="2.0.1" PACKAGE_VERSION_MAJOR=2
                             PACKAGE_VERSION_MINOR=0 PACKAGE_VERSION_MICRO=1 SOFA_VERSION="20231011")
  set_target_properties(erfa_bundled PROPERTIES POSITION_INDEPENDENT_CODE ON)
  if(UNIX)
    target_link_libraries(erfa_bundled PUBLIC m)
  endif()
  set(TRANSIT_ERFA_TARGET erfa_bundled)
endif()

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
  src/export.cpp
  src/io/reader_factory.cpp)
target_include_directories(transit PUBLIC include)
target_include_directories(transit SYSTEM PRIVATE third_party/tomlplusplus)
target_compile_definitions(transit PRIVATE TRANSIT_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
find_package(Threads REQUIRED)
target_link_libraries(transit PUBLIC ${TRANSIT_ERFA_TARGET} Threads::Threads)
target_compile_options(transit PRIVATE -Wall -Wextra -Wpedantic)

if(TRANSIT_WITH_CFITSIO)
  # pkg-config first; its paths are checked because MSYS2's cfitsio.pc gives MSYS paths
  # (/ucrt64/include) that a native CMake cannot use. Otherwise look for the header and library.
  pkg_check_modules(CFITSIO IMPORTED_TARGET cfitsio)
  set(CFITSIO_PC_USABLE FALSE)
  if(TARGET PkgConfig::CFITSIO)
    set(CFITSIO_PC_USABLE TRUE)
    get_target_property(CFITSIO_PC_INCLUDES PkgConfig::CFITSIO INTERFACE_INCLUDE_DIRECTORIES)
    if(CFITSIO_PC_INCLUDES)
      foreach(dir IN LISTS CFITSIO_PC_INCLUDES)
        if(NOT EXISTS "${dir}")
          set(CFITSIO_PC_USABLE FALSE)
        endif()
      endforeach()
    endif()
  endif()
  if(CFITSIO_PC_USABLE)
    target_link_libraries(transit PUBLIC PkgConfig::CFITSIO)
  else()
    find_path(CFITSIO_INCLUDE_DIR fitsio.h PATH_SUFFIXES cfitsio)
    find_library(CFITSIO_LIBRARY NAMES cfitsio libcfitsio)
    if(NOT CFITSIO_INCLUDE_DIR OR NOT CFITSIO_LIBRARY)
      message(FATAL_ERROR "CFITSIO not found: install it (libcfitsio-dev, brew cfitsio, MSYS2 mingw-w64-ucrt-x86_64-cfitsio) "
                          "or set CFITSIO_INCLUDE_DIR and CFITSIO_LIBRARY")
    endif()
    message(STATUS "CFITSIO: ${CFITSIO_LIBRARY} (headers in ${CFITSIO_INCLUDE_DIR})")
    target_include_directories(transit PUBLIC ${CFITSIO_INCLUDE_DIR})
    target_link_libraries(transit PUBLIC ${CFITSIO_LIBRARY})
  endif()
  target_sources(transit PRIVATE src/io/cfitsio_reader.cpp)
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

include(GNUInstallDirs)
install(TARGETS transitscope RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
# limb-darkening tables and passbands, found by the programs at <prefix>/share/transitscope/data
install(DIRECTORY data/ DESTINATION ${CMAKE_INSTALL_DATADIR}/transitscope/data)
install(FILES README.md LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/transitscope)
install(DIRECTORY examples/ DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/transitscope/examples)

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
  set(GUI_SOURCES gui/main.cpp gui/MainWindow.cpp gui/SettingsEditor.cpp gui/ImageView.cpp gui/PlotWidget.cpp
                  gui/ResultsView.cpp)
  if(WIN32)
    list(APPEND GUI_SOURCES packaging/transitscope.rc)               # the .exe icon
  elseif(APPLE)
    list(APPEND GUI_SOURCES packaging/transitscope.icns)
    set_source_files_properties(packaging/transitscope.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
  endif()
  # a Windows GUI program (no console window) and a macOS application bundle
  add_executable(transitscope-gui WIN32 MACOSX_BUNDLE ${GUI_SOURCES})
  set_target_properties(transitscope-gui PROPERTIES
    MACOSX_BUNDLE_BUNDLE_NAME TransitScope
    MACOSX_BUNDLE_GUI_IDENTIFIER io.github.jkobierczynski.transitscope
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION}
    MACOSX_BUNDLE_ICON_FILE transitscope.icns)
  qt_add_resources(transitscope-gui "icon" PREFIX "/" BASE packaging FILES packaging/transitscope.png)
  target_include_directories(transitscope-gui SYSTEM PRIVATE third_party/tomlplusplus)
  target_compile_definitions(transitscope-gui PRIVATE TRANSITSCOPE_VERSION="${PROJECT_VERSION}")
  target_link_libraries(transitscope-gui PRIVATE transit Qt6::Widgets Qt6::Concurrent)
  target_compile_options(transitscope-gui PRIVATE -Wall -Wextra)
  install(TARGETS transitscope-gui RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} BUNDLE DESTINATION .)
  if(UNIX AND NOT APPLE)
    install(FILES packaging/transitscope.desktop DESTINATION ${CMAKE_INSTALL_DATADIR}/applications)
    install(FILES packaging/transitscope.png DESTINATION ${CMAKE_INSTALL_DATADIR}/icons/hicolor/256x256/apps)
  endif()

  find_package(Qt6 COMPONENTS Test)
  if(Qt6Test_FOUND)
    add_executable(gui_test tests/gui/gui_test.cpp
      gui/MainWindow.cpp gui/SettingsEditor.cpp gui/ImageView.cpp gui/PlotWidget.cpp gui/ResultsView.cpp)
    target_include_directories(gui_test SYSTEM PRIVATE third_party/tomlplusplus)
    target_compile_definitions(gui_test PRIVATE TRANSITSCOPE_VERSION="${PROJECT_VERSION}")
    target_link_libraries(gui_test PRIVATE transit Qt6::Widgets Qt6::Concurrent Qt6::Test)
  endif()
endif()

add_executable(comp_stability_test tests/unit/comp_stability_test.cpp)
target_link_libraries(comp_stability_test PRIVATE transit)
