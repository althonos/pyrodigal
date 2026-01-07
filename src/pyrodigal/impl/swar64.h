#ifndef _PYRODIGAL_IMPL_SWAR64_H
#define _PYRODIGAL_IMPL_SWAR64_H

#include <stdint.h>

#include "sequence.h"
#include "generic.h"
#include "template.h"

#define SWAR_HIBIT        0x8080808080808080ULL
#define SWAR_LOBIT        0x0101010101010101ULL

#define simd_t            uint64_t
#define simd_load(m)      swar_load(m)
#define simd_store(x, m)  swar_store(x, m)
#define simd_set1(x)      (SWAR_LOBIT * ((uint64_t) x))
#define simd_eq(x, y)     ((~(((x ^ y) | SWAR_HIBIT) - SWAR_LOBIT) & SWAR_HIBIT) >> 7)

#define mask_t            simd_t
#define mask_eq(x, y)     simd_eq(x, y)
#define mask_or(x, y)     (x | y)
#define mask_and(x, y)    (x & y)
#define mask_andnot(x, y) (x & (~y))
#define mask_convert(x)   (x & SWAR_LOBIT)

#define SIMD_LANES 8
#define SIMD_MASK  0x7

static inline simd_t swar_load(const uint8_t* data) {
    return *((const simd_t*) data);
}

static inline void swar_store(simd_t x, uint8_t* data) {
    *((simd_t*) data) = x;
}

void skippable_swar64(
    const uint8_t* restrict strands,
    const uint8_t* restrict types,
    const uint8_t* restrict frames,
    const int min,
    const int i,
          uint8_t* restrict skip
) {
    skippable_simd(strands, types, frames, min, i, skip);
}

#endif
