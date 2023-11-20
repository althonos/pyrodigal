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
#define simd_eq(x, y)     _mm512_cmpeq_epi8_mask(x, y)

#define mask_t            __mmask64
#define mask_or(x, y)     _kor_mask64(x, y)
#define mask_and(x, y)    _kand_mask64(x, y)
#define mask_andnot(x, y) _kandn_mask64(y, x)
#define mask_convert(x)   _mm512_movm_epi8(x)

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
