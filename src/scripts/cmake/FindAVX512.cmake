#.rst:
# FindAVX512
# ----------
#
# Finds AVX512 support
#
# This module can be used to detect AVX512 support in a C compiler.  If
# the compiler supports AVX512, the flags required to compile with
# AVX512 support are returned in variables for the different languages.
# The variables may be empty if the compiler does not need a special
# flag to support AVX512.
#
# The following variables are set:
#
# ::
#
#    AVX512_C_FLAGS - flags to add to the C compiler for AVX512 support
#    AVX512_FOUND - true if AVX512 is detected
#
#=============================================================================

set(_AVX512_REQUIRED_VARS)
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set(CMAKE_REQUIRED_QUIET ${AVX512_FIND_QUIETLY})

# sample AVX512 source code to test
set(AVX512_C_TEST_SOURCE
"
#include <immintrin.h>
int foo() {
    __m512i vOne     = _mm512_set1_epi8(1);
    __mmask64 result = _mm512_cmpeq_epi8_mask(vOne,vOne);
    return result;
}
int main(void) { return (int)foo(); }
")

# if these are set then do not try to find them again,
# by avoiding any try_compiles for the flags
if((DEFINED AVX512_C_FLAGS) OR (DEFINED HAVE_AVX512))
else()
  if(CMAKE_C_COMPILER_ID STREQUAL "MSVC")
    # MSVC can compile AVX intrinsics without the arch flag, however it
    # will detect that AVX code is found and "consider using /arch:AVX".
    set(AVX512_C_FLAG_CANDIDATES
      "/arch:AVX512")
  else()
    set(AVX512_C_FLAG_CANDIDATES
      #Empty, if compiler automatically accepts AVX512
      " "
      #clang, gcc
      "-mavx512f -mavx512bw"
    )
  endif()

  include(CheckCSourceCompiles)

  foreach(FLAG IN LISTS AVX512_C_FLAG_CANDIDATES)
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(HAVE_AVX512 CACHE)
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Try AVX512 C flag = [${FLAG}]")
    endif()
    check_c_source_compiles("${AVX512_C_TEST_SOURCE}" HAVE_AVX512)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    if(HAVE_AVX512)
      set(AVX512_C_FLAGS_INTERNAL "${FLAG}")
      break()
    endif()
  endforeach()

  unset(AVX512_C_FLAG_CANDIDATES)
  
  set(AVX512_C_FLAGS "${AVX512_C_FLAGS_INTERNAL}"
    CACHE STRING "C compiler flags for AVX512 intrinsics")
endif()

list(APPEND _AVX512_REQUIRED_VARS AVX512_C_FLAGS)

set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if(_AVX512_REQUIRED_VARS)
  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args(AVX512
                                    REQUIRED_VARS ${_AVX512_REQUIRED_VARS})

  mark_as_advanced(${_AVX512_REQUIRED_VARS})

  unset(_AVX512_REQUIRED_VARS)
else()
  message(SEND_ERROR "FindAVX512 requires C or CXX language to be enabled")
endif()
