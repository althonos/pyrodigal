#include "sequence.h"
#include "mmx.h"
#include "generic.h"

#ifdef __MMX__

#include <mmintrin.h>

#include "template.h"

#define simd_t            __m64
#define simd_load(m)      (((__m64*) (m))[0])
#define simd_store(x, m)  (((__m64*) (m))[0] = x)
#define simd_set1(x)      _mm_set1_pi8(x)
#define simd_eq(x, y)     _mm_cmpeq_pi8(x, y)

#define mask_t            __m64
#define mask_or(x, y)     _mm_or_si64(x, y)
#define mask_and(x, y)    _mm_and_si64(x, y)
#define mask_andnot(x, y) _mm_andnot_si64(y, x)
#define mask_convert(x)   (x)

#define SIMD_LANES 8
#define SIMD_MASK  0x7

void skippable_mmx(
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
