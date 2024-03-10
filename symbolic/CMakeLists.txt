cmake_minimum_required(VERSION 3.12...3.28)

include(AICxxProject)

#==============================================================================
# BUILD OBJECT LIBRARY
#

# This project is an OBJECT-library, used by other git submodules and the main project.
add_library(symbolic_ObjLib OBJECT)

# Set SOURCE_DIR to the current directory.
set_property(
  TARGET symbolic_ObjLib
  PROPERTY SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
)

file(GLOB HEADER_FILES "*.h")

# The list of source files.
target_sources(symbolic_ObjLib
  PRIVATE
    "Constant.cxx"

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
  PUBLIC
    AICxx::utils
)

# Create an ALIAS target.
add_library(AICxx::symbolic ALIAS symbolic_ObjLib)

# Testsuite.
add_executable(symbolic_testsuite symbolic_testsuite.cxx)
target_link_libraries(symbolic_testsuite AICxx::symbolic ${AICXX_OBJECTS_LIST})
