cmake_minimum_required(VERSION 3.20)

if(DEFINED SKBUILD_PROJECT_NAME)
    project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION})
else()
    project(pyncbitk)
endif()

enable_language(C)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set_property(GLOBAL PROPERTY PYTHON_EXTENSIONS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# --- Detect SIMD --------------------------------------------------------------

include("src/scripts/cmake/FindAVX2.cmake")
include("src/scripts/cmake/FindNEON.cmake")
include("src/scripts/cmake/FindSSE2.cmake")
include("src/scripts/cmake/CythonExtension.cmake")

# --- Compile Prodigal code ----------------------------------------------------

include_directories(BEFORE PUBLIC ${CMAKE_CURRENT_LIST_DIR}/pyrodigal/prodigal)
include_directories(AFTER  PUBLIC include/prodigal)

# --- Compile Python extension -------------------------------------------------

add_subdirectory("src")
