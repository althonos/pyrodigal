cmake_minimum_required(VERSION 3.20)
project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES C)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set_property(GLOBAL PROPERTY PYTHON_EXTENSIONS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# --- Detect SIMD --------------------------------------------------------------

include("src/scripts/cmake/FindAVX2.cmake")
include("src/scripts/cmake/FindNEON.cmake")
include("src/scripts/cmake/FindSSE2.cmake")
include("src/scripts/cmake/CythonExtension.cmake")

# --- Compile code -------------------------------------------------------------

add_subdirectory("src")
