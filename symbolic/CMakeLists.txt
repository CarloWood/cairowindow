cmake_minimum_required(VERSION 3.12...3.28)

include(AICxxProject)

#==============================================================================
# OPTIONS

# Option 'EnableSymbolicPrinting' is ${OptionEnableLibcwd} by default.
cw_option(EnableSymbolicPrinting
  "Enable printing of symbolic expressions to an ostream" ${OptionEnableLibcwd}
  "" OFF
)

message(DEBUG "OptionEnableSymbolicPrinting is ${OptionEnableSymbolicPrinting}")
if (OptionEnableSymbolicPrinting)
  set(SYMBOLIC_PRINTING 1)
endif ()

# Specify symbolic specific configure file.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/config.h
               @ONLY)

#==============================================================================
# BUILD OBJECT LIBRARY
#

# This project is an OBJECT-library, used by other git submodules and the main project.
add_library(symbolic_ObjLib OBJECT)

# Tell sys.h that we have a config.h.
target_compile_definitions(symbolic_ObjLib
  PUBLIC
    HAVE_SYMBOLIC_CONFIG_H
)

file(GLOB HEADER_FILES "*.h")

# The list of source files.
target_sources(symbolic_ObjLib
  PRIVATE
    "Expression.cxx"
    "Constant.cxx"
    "Symbol.cxx"
    "Sum.cxx"
    "Product.cxx"
    "Power.cxx"
    "Exponential.cxx"
    "Sin.cxx"
    "Cos.cxx"
    "Atan.cxx"
    "Log.cxx"
    "iomanip_fulldef.cxx"

    ${HEADER_FILES}
)

# Required include search-paths.
get_target_property(CWDS_INTERFACE_INCLUDE_DIRECTORIES AICxx::cwds INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(symbolic_ObjLib
  PUBLIC
    "${CWDS_INTERFACE_INCLUDE_DIRECTORIES}" # For sys.h and debug.h.
)

# Require support for C++20.
target_compile_features(symbolic_ObjLib
  PUBLIC cxx_std_20
)

# Set link dependencies.
target_link_libraries(symbolic_ObjLib
  PRIVATE
    AICxx::utils
)

# Create an ALIAS target.
add_library(AICxx::symbolic ALIAS symbolic_ObjLib)

# Testsuite.
add_executable(symbolic_testsuite symbolic_testsuite.cxx)
target_link_libraries(symbolic_testsuite AICxx::symbolic ${AICXX_OBJECTS_LIST})
