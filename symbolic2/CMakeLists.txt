cmake_minimum_required(VERSION 3.12...3.28)

include(AICxxProject)

#==============================================================================
# OPTIONS

# Option 'EnableSymbolic2Printing' is ${OptionEnableLibcwd} by default.
cw_option(EnableSymbolic2Printing
  "Enable printing of symbolic expressions to an ostream" ${OptionEnableLibcwd}
  "" OFF
)

message(DEBUG "OptionEnableSymbolic2Printing is ${OptionEnableSymbolic2Printing}")
if (OptionEnableSymbolic2Printing)
  set(SYMBOLIC2_PRINTING 1)
endif ()

# Specify utils specific configure file.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/config.h
               @ONLY)

#==============================================================================
# BUILD OBJECT LIBRARY
#

# This project is an OBJECT-library, used by other git submodules and the main project.
add_library(symbolic2_ObjLib OBJECT)

# Tell sys.h that we have a config.h.
target_compile_definitions(symbolic2_ObjLib
  PUBLIC
    HAVE_SYMBOLIC2_CONFIG_H
)

file(GLOB HEADER_FILES "*.h")

# The list of source files.
target_sources(symbolic2_ObjLib
  PRIVATE
    "Expression.cxx"
    "Constant.cxx"
    "Symbol.cxx"

    ${HEADER_FILES}
)

# Required include search-paths.
get_target_property(CWDS_INTERFACE_INCLUDE_DIRECTORIES AICxx::cwds INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(symbolic2_ObjLib
  PUBLIC
    "${CWDS_INTERFACE_INCLUDE_DIRECTORIES}" # For sys.h and debug.h.
)

# Require support for C++20.
target_compile_features(symbolic2_ObjLib
  PUBLIC cxx_std_20
)

# Set link dependencies.
target_link_libraries(symbolic2_ObjLib
  PUBLIC
    AICxx::utils
)

# Create an ALIAS target.
add_library(AICxx::symbolic2 ALIAS symbolic2_ObjLib)

# Testsuite.
add_executable(symbolic2_testsuite symbolic2_testsuite.cxx)
target_link_libraries(symbolic2_testsuite AICxx::symbolic2 ${AICXX_OBJECTS_LIST})
