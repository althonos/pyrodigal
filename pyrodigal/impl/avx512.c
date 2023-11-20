#include "sequence.h"
#include "avx512.h"
#include "generic.h"

#ifdef __AVX2__

#include <immintrin.h>

#include "template.h"

#define simd_t            __m512i
#define simd_load(m)      _mm512_load_si512((__m512i*) (m))
#define simd_store(x, m)  _mm512_store_si512((__m512i*) (m), x)
#define simd_set1(x)      _mm512_set1_epi8(x)
#define simd_or(x, y)     _mm512_or_si512(x, y)
#define simd_eq(x, y)     _mm512_movm_epi8(_mm512_cmpeq_epi8_mask(x, y))
#define simd_and(x, y)    _mm512_and_si512(x, y)
#define simd_andnot(x, y) _mm512_andnot_si512(y, x)

#define SIMD_LANES 64
#define SIMD_MASK  0x3F

void skippable_avx512(
    const int8_t*  restrict strands,
    const uint8_t* restrict types,
    const uint8_t* restrict frames,
    const int min,
    const int i,
          uint8_t* restrict skip
) {
    skippable_simd(strands, types, frames, min, i, skip);
}

#endif
