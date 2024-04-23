cmake_minimum_required(VERSION 3.12...3.28)

project(cairowindow
  DESCRIPTION "Window to draw in with cairo"
  LANGUAGES CXX C       # Language "C" is required for find_package(Threads).
)

include(AICxxProject)

#==============================================================================
# DEPENDENCIES
#

find_package(PkgConfig REQUIRED)
pkg_check_modules(CAIRO REQUIRED cairo)
pkg_check_modules(X11 REQUIRED x11)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

pkg_check_modules(PANGOCAIRO pangocairo)

find_package_handle_standard_args(
  PANGOCAIRO
  DEFAULT_MSG
  PANGOCAIRO_INCLUDE_DIRS
  PANGOCAIRO_LIBRARIES
)

mark_as_advanced(
  PANGOCAIRO_INCLUDE_DIRS
  PANGOCAIRO_LIBRARIES
)

#==============================================================================
# OPTIONS

# Option 'EnableCairoWindowTests' is OFF by default.
cw_option(EnableCairoWindowTests
  "Compile the tests" OFF
  "" ON
)

# Option 'EnableDebugWindow' is OFF by default.
cw_option(EnableDebugWindow
  "Enable debug windows that show the content of Layer surfaces and the off-screen surface" OFF
  "CW_BUILD_TYPE_IS_NOT_RELEASE" OFF
)

if (CW_BUILD_TYPE_IS_NOT_RELEASE)
  message(DEBUG "OptionEnableDebugWindow is ${OptionEnableDebugWindow}")
endif ()
if (OptionEnableDebugWindow)
  set(CAIROWINDOW_DEBUGWINDOW 1)
endif ()

#==============================================================================
# PLATFORM SPECIFIC CHECKS
#

# Specify utils specific configure file.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/config.h
               @ONLY)

#==============================================================================
# BUILD OBJECT LIBRARY
#

# This project is an OBJECT-library, used by other git submodules and the main project.
add_library(cairowindow_ObjLib OBJECT)

# Tell cwds/sys.h that we have a config.h.
target_compile_definitions(cairowindow_ObjLib
  PUBLIC
    HAVE_CAIROWINDOW_CONFIG_H
)

target_compile_options(cairowindow_ObjLib
  INTERFACE
    ${CAIRO_CFLAGS} ${X11_CFLAGS}
)

file(GLOB HEADER_FILES "*.h" "draw/*.h")

# The list of source files.
target_sources(cairowindow_ObjLib
  PRIVATE
    "BezierFitter.cxx"
    "BezierCurve.cxx"
    "Color.cxx"
    "Direction.cxx"
    "EventLoop.cxx"
    "FontPangoContext.cxx"
    "Layer.cxx"
    "LayerRegion.cxx"
    "Line.cxx"
    "LinePiece.cxx"
    "Point.cxx"
    "Pixel.cxx"
    "Plot.cxx"
    "Printable.cxx"
    "QuadraticEnergy.cxx"
    "Window.cxx"
    "draw/ArrowHead.cxx"
    "draw/BezierFitter.cxx"
    "draw/Connector.cxx"
    "draw/Curve.cxx"
    "draw/Grid.cxx"
    "draw/Line.cxx"
    "draw/PlotArea.cxx"
    "draw/Shape.cxx"
    "draw/Slider.cxx"
    "draw/ChessDiagram.cxx"

    ${HEADER_FILES}
)

if (OptionEnableDebugWindow)
  target_sources(cairowindow_ObjLib
    PRIVATE
      "DebugWindow.cxx"
      "DebugWindow.h"
  )
endif ()

if (OptionEnableLibcwd)
  target_sources(cairowindow_ObjLib
    PRIVATE
      "debugcairo.cxx"
      "Message.cxx"
  )
endif ()

# Required include search-paths.
get_target_property(CWDS_INTERFACE_INCLUDE_DIRECTORIES AICxx::cwds INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(cairowindow_ObjLib
  PUBLIC
    "${CWDS_INTERFACE_INCLUDE_DIRECTORIES}" # For sys.h and debug.h.
    "${top_objdir}"                         # For cairowindow/config.h.
    "${PANGOCAIRO_INCLUDE_DIRS}"
)

# Require support for C++20.
target_compile_features(cairowindow_ObjLib
  PUBLIC cxx_std_20
)

# Set link dependencies.
target_link_libraries(cairowindow_ObjLib
  PUBLIC
    AICxx::symbolic
    AICxx::chess
    AICxx::utils
    Eigen3::Eigen
  INTERFACE
    ${PANGOCAIRO_LIBRARIES}
    ${X11_LDFLAGS}
)

# Create an ALIAS target.
add_library(AICxx::cairowindow ALIAS cairowindow_ObjLib)

add_subdirectory(symbolic)
add_subdirectory(chess)

if (OptionEnableCairoWindowTests)
add_subdirectory(tests)
endif ()
