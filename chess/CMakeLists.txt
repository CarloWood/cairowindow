cmake_minimum_required(VERSION 3.12...3.28)

include(AICxxProject)

#==============================================================================
# BUILD OBJECT LIBRARY
#

# This project is an OBJECT-library, used by other git submodules and the main project.
add_library(chess_ObjLib OBJECT)

file(GLOB HEADER_FILES "*.h")

# The list of source files.
target_sources(chess_ObjLib
  PRIVATE
    "Diagram.cxx"

    ${HEADER_FILES}
)

# Required include search-paths.
get_target_property(CWDS_INTERFACE_INCLUDE_DIRECTORIES AICxx::cwds INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(chess_ObjLib
  PUBLIC
    "${CWDS_INTERFACE_INCLUDE_DIRECTORIES}" # For sys.h and debug.h.
)

# Require support for C++20.
target_compile_features(chess_ObjLib
  PUBLIC cxx_std_20
)

# Set link dependencies.
target_link_libraries(chess_ObjLib
  PUBLIC
    AICxx::math
  PRIVATE
    AICxx::utils
)

# Create an ALIAS target.
add_library(AICxx::chess ALIAS chess_ObjLib)

# Testsuite.
#add_executable(chess_testsuite chess_testsuite.cxx)
#target_link_libraries(chess_testsuite AICxx::chess ${AICXX_OBJECTS_LIST})
