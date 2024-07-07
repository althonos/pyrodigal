project(
  'pyrodigal',
  'c', 'cython',
  version : '3.4.1',
  # default_options : ['warning_level=0']
)

# load filesystem module
fs = import('fs')

# Load C compiler and check if `-lm` is needed
cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)

# load Python dependencies
py_mod = import('python')
python = py_mod.find_installation(pure: false)
py_dep = python.dependency()


# --- SIMD detection -----------------------------------------------------------

# ---- MMX ----

mmx_code = (
  '''
  #include <mmintrin.h>
  int main(int argc, char *argv[]) {
      __m64 a = _mm_set1_pi16(1);
      short x = (short) _m_to_int(a);
      return (x == 1) ? 0 : 1;
  }
  '''
)
mmx_supported = cc.compiles(mmx_code, name : 'MMX', args : [])
mmx_flags = ['-D__MMX__']

# ---- SSE2 ----

sse2_code = (
  '''
  #include <emmintrin.h>
  int main(int argc, char *argv[]) {
      __m128i a = _mm_set1_epi16(-1);
              a = _mm_and_si128(a, a);
      short   x = _mm_extract_epi16(a, 1);
      return (x == -1) ? 0 : 1;
  }
  '''
)
sse2_supported = cc.compiles(sse2_code, name : 'SSE2', args : [])
sse2_flags = ['-D__SSE2__']
if not sse2_supported
  sse2_supported = cc.compiles(sse2_code, name : 'SSE with -msse2', args : ['-msse2'])
  if sse2_supported
    sse2_flags += '-msse2'
  endif
endif

# ---- AVX2 ----

avx2_code = (
  '''
  #include <immintrin.h>
  int main(int argc, char *argv[]) {
      __m256i a = _mm256_set1_epi16(-1);
              a = _mm256_abs_epi16(a);
      short   x = _mm256_extract_epi16(a, 1);
      return (x == 1) ? 0 : 1;
  }
  '''
)
avx2_supported = cc.compiles(avx2_code, name : 'AVX2', args : [])
avx2_flags = ['-D__AVX2__']
if not avx2_supported
  avx2_supported = cc.compiles(avx2_code, name : 'AVX2 with -mavx2', args : ['-mavx2'])
  if avx2_supported
    avx2_flags += '-mavx2'
  endif
endif

# ---- AVX512 ----

avx512_code = (
  '''
  #include <x86intrin.h>
  int main(int argc, char *argv[]) {
      __m512i   a = _mm512_set1_epi16(-1);
      __m512i   b = _mm512_set1_epi16(0);
                a = _mm512_abs_epi16(a);
      __mmask32 c = _mm512_cmpgt_epi16_mask(a, b);
      return (c == 0xFFFFFFFF) ? 0 : 1;
  }
  '''
)
avx512_supported = cc.compiles(avx512_code, name : 'AVX512', args : [])
avx512_flags = ['-D__AVX512F__']
if not avx512_supported
  avx512_supported = cc.compiles(avx512_code, name : 'AVX512 with -mavx512bw -mavx512f', args : ['-mavx512f', '-mavx512bw'])
  if avx512_supported
    avx512_flags += ['-mavx512f', '-mavx512bw']
  endif
endif

# ---- NEON ----

neon_code = (
  '''
  #include <arm_neon.h>
  int main(int argc, char *argv[]) {
      int16x8_t a = vdupq_n_s16(-1);
                a = vabsq_s16(a);
      short     x = vgetq_lane_s16(a, 1);
      return (x == 1) ? 0 : 1;
  }
  '''
)
neon_supported = cc.compiles(neon_code, name : 'NEON', args : [])
neon_flags = ['-D__ARM_NEON__']
if not neon_supported
  neon_supported = cc.compiles(neon_code, name : 'NEON with -mfpu=neon', args : ['-mfpu=neon'])
  if neon_supported
    neon_flags += '-mfpu=neon'
  endif
endif


# --- Prodigal -----------------------------------------------------------------

training_splits = ['training_header.c']
foreach i : range(50)
  training_splits += 'training@0@.c'.format(i)
endforeach

# split training.c
training = custom_target(
    'training',
    input : 'vendor/Prodigal/training.c',
    output : training_splits,
    command : [
      python,
      meson.current_source_dir() / 'tools' / 'split_training.py',
      '--input',
      '@INPUT@',
      '--output-dir',
      '@OUTDIR@'
    ],
)

# build Prodigal as a library
prodigal_sources = [
    'vendor/Prodigal/bitmap.c',
    'vendor/Prodigal/dprog.c',
    'vendor/Prodigal/gene.c',
    'vendor/Prodigal/metagenomic.c',
    'vendor/Prodigal/node.c',
    'vendor/Prodigal/sequence.c',
]

# copy files out of source so that includes don't have locality precedence
prodigal_files = []
foreach p : prodigal_sources
  prodigal_files += configure_file(
    input : p, 
    output : fs.name(p), 
    copy : true,
  )
endforeach

prodigal_inc = include_directories(
  'pyrodigal/prodigal',
  'vendor/Prodigal',
)

prodigal = static_library('prodigal',
    prodigal_files,
    dependencies: m_dep,
    include_directories : prodigal_inc,
)

prodigal_training = static_library('training',
    training,
    include_directories : prodigal_inc
)

# --- Pyrodigal ----------------------------------------------------------------

# check if `PyInterpreterState_GetID` is available
getid_code = (
  '''
  #include <stdint.h>
  #include <stdlib.h>
  #include <Python.h>

  int main(int argc, char *argv[]) {
      PyInterpreterState_GetID(NULL);
      return 0;
  }
  '''
)
has_getid = cc.compiles(
  getid_code,
  dependencies : py_dep,
  name :  'PyInterpreterState_GetID',
  args : ['-Werror=implicit-function-declaration'],
)

subdir('pyrodigal')

