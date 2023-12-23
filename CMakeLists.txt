cmake_minimum_required(VERSION 3.12...3.27)

project(cairowindow
  DESCRIPTION "Window to draw in with cairo"
  LANGUAGES CXX C       # Language "C" is required for find_package(Threads).
)

include(AICxxProject)

#==============================================================================
# DEPENDENCIES
#

#==============================================================================
# OPTIONS

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

# Set SOURCE_DIR to the current directory.
set_property(
  TARGET cairowindow_ObjLib
  PROPERTY SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
)

# The list of source files.
target_sources(cairowindow_ObjLib
  PRIVATE
    "Layer.cxx"
    "Window.cxx"
    "EventLoop.cxx"
    "LayerRegion.cxx"
    "Plot.cxx"
    "draw/PlotArea.cxx"
    "draw/Grid.cxx"

    "Color.h"
    "Rectangle.h"
    "Layer.h"
    "Window.h"
    "EventLoop.h"
    "LayerRegion.h"
)

if (OptionEnableDebugWindow)
  target_sources(cairowindow_ObjLib
    PRIVATE
      "DebugWindow.cxx"
      "DebugWindow.h"
  )
endif ()

# Required include search-paths.
get_target_property(CWDS_INTERFACE_INCLUDE_DIRECTORIES AICxx::cwds INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(cairowindow_ObjLib
  PUBLIC
    "${CWDS_INTERFACE_INCLUDE_DIRECTORIES}" # For sys.h and debug.h.
    "${top_objdir}"                         # For cairowindow/config.h.
)

# Require support for C++20.
target_compile_features(cairowindow_ObjLib
  PUBLIC cxx_std_20
)

# Set link dependencies.
target_link_libraries(cairowindow_ObjLib
  PUBLIC
    AICxx::utils
)

# Create an ALIAS target.
add_library(AICxx::cairowindow ALIAS cairowindow_ObjLib)

# Prepend this object library to the list.
set(AICXX_OBJECTS_LIST AICxx::cairowindow ${AICXX_OBJECTS_LIST} CACHE INTERNAL "List of OBJECT libaries that this project uses.")
