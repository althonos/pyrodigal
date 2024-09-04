find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
get_property(PYTHON_EXTENSIONS_SOURCE_DIR GLOBAL PROPERTY PYTHON_EXTENSIONS_SOURCE_DIR)

# --- Detect PyInterpreterState_GetID ------------------------------------------

include(CheckCSourceCompiles)

set(SAFE_CMAKE_REQUIRED_INCLUDES "${CMAKE_REQUIRED_INCLUDES}")
set(CMAKE_REQUIRED_INCLUDES "${Python_INCLUDE_DIRS}")
set(PYINTERPRETER_STATE_SOURCE
"
#include <stdint.h>
    #include <stdlib.h>
    #include <Python.h>

    int main(int argc, char *argv[]) {{
      PyInterpreterState_GetID(NULL);
      return 0;
    }
")
check_c_source_compiles("${PYINTERPRETER_STATE_SOURCE}" HAVE_PYINTERPRETERSTATE_GETID)
set(CMAKE_REQUIRED_INCLUDES "${SAFE_CMAKE_REQUIRED_INCLUDES}")

# --- Prepare Cython directives and constants ----------------------------------

set(CYTHON_DIRECTIVES
    -X cdivision=True
    -X nonecheck=False
    -E SSE2_BUILD_SUPPORT=$<IF:$<BOOL:${HAVE_SSE2}>,True,False>
    -E AVX2_BUILD_SUPPORT=$<IF:$<BOOL:${HAVE_AVX2}>,True,False>
    -E NEON_BUILD_SUPPORT=$<IF:$<BOOL:${HAVE_NEON}>,True,False>
    -E MMX_BUILD_SUPPORT=False
    -E AVX512_BUILD_SUPPORT=False
    -E SYS_IMPLEMENTATION_NAME="cpython"
    -E SYS_VERSION_INFO_MAJOR=${Python_VERSION_MAJOR}
    -E SYS_VERSION_INFO_MINOR=${Python_VERSION_MINOR}
    -E TARGET_CPU=${CMAKE_SYSTEM_PROCESSOR}
    -E TARGET_SYSTEM="linux"
    -E PYPY=$<IF:$<STREQUAL:${Python_INTERPRETER_ID},PyPy>,True,False>
    -E PROJECT_VERSION=${CMAKE_PROJECT_VERSION}
    -E HAVE_PYINTERPRETERSTATE_GETID=$<IF:$<BOOL:${HAVE_PYINTERPRETERSTATE_GETID}>,True,False>
)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set(CYTHON_DIRECTIVES
    ${CYTHON_DIRECTIVES}
    -X cdivision_warnings=True
    -X warn.undeclared=True
    -X warn.unreachable=True
    -X warn.maybe_uninitialized=True
    -X warn.unused=True
    -X warn.unused_arg=True
    -X warn.unused_result=True
    -X warn.multiple_declarators=True
  )
else()
  set(CYTHON_DIRECTIVES
    ${CYTHON_DIRECTIVES}
    -X boundscheck=True
    -X wraparound=True
  )
endif()

# --- Declare Cython extension -------------------------------------------------

macro(cython_extension _name)
  cmake_parse_arguments(CYTHON_EXTENSION "" "" "LINKS" ${ARGN} )

  # Make sure that the source directory is known
  if(NOT DEFINED PYTHON_EXTENSIONS_SOURCE_DIR)
    message(FATAL_ERROR "The PYTHON_EXTENSIONS_SOURCE_DIR variable has not been set.")
  endif()

  # Generate C file from Cython file
  add_custom_command(
    OUTPUT ${_name}.c
    COMMENT
      "Cythonizing ${_name}.pyx to ${_name}.c"
    COMMAND
      Python::Interpreter -m cython
        "${CMAKE_CURRENT_SOURCE_DIR}/${_name}.pyx"
        --output-file ${_name}.c
        -I "${CYTHON_HEADERS_DIR}"
        ${CYTHON_DIRECTIVES}
    MAIN_DEPENDENCY
      ${_name}.pyx
    VERBATIM)

  # Build fully-qualified module name as the target name
  string(REGEX REPLACE "^${PYTHON_EXTENSIONS_SOURCE_DIR}/?" "" _dest_folder ${CMAKE_CURRENT_SOURCE_DIR})
  string(REPLACE "/" "." _target ${_dest_folder}.${_name})

  # Add Python module
  python_add_library(${_target} MODULE WITH_SOABI ${_name}.pyx ${_name}.pxd ${_name}.c)
  set_target_properties(${_target} PROPERTIES OUTPUT_NAME ${_name} )
  target_include_directories(${_target} AFTER PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})  
  target_link_libraries(${_target} PUBLIC ${CYTHON_EXTENSION_LINKS})

  # Preserve the relative project structure in the install directory
  string(REGEX REPLACE "^${PYTHON_EXTENSIONS_SOURCE_DIR}/?" "" _dest_folder ${CMAKE_CURRENT_SOURCE_DIR})
  install(TARGETS ${_target} DESTINATION ${_dest_folder} )
  message(DEBUG "Install folder for extension ${_name}: ${_dest_folder}")

  # Add the targets to the list of Cython extensions
  get_property(_ext GLOBAL PROPERTY PYRODIGAL_CYTHON_EXTENSIONS)
  list(APPEND _ext ${_target})
  set_property(GLOBAL PROPERTY PYRODIGAL_CYTHON_EXTENSIONS ${_ext})
endmacro()
