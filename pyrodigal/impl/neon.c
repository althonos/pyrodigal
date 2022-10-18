#include "sequence.h"
#include "neon.h"
#include "generic.h"

#ifdef __ARM_NEON__

#include <arm_neon.h>

#include "template.h"

#define simd_t            uint8x16_t
#define simd_load(m)      vld1q_u8((uint8_t*) (m))
#define simd_store(x, m)  vst1q_u8((uint8_t*) (m), x)
#define simd_set1(x)      vdupq_n_u8(x)
#define simd_or(x, y)     vorrq_u8(x, y)
#define simd_eq(x, y)     vceqq_u8(x, y)
#define simd_and(x, y)    vandq_u8(x, y)
#define simd_andnot(x, y) vbicq_u8(x, y)

#define SIMD_LANES 16
#define SIMD_MASK  0xF

void skippable_neon(
    const int8_t* strands,
    const uint8_t* types,
    const uint8_t* frames,
    const int min,
    const int i,
    uint8_t* skip
) {
    skippable_simd(strands, types, frames, min, i, skip);
}

#endif
