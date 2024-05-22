#define skippable_simd(strands, types, frames, min, i, skip) do {         \
        const simd_t ALL_STOPS = simd_set1(STOP);                         \
        const simd_t ALL_FWD   = simd_set1(1);                            \
                                                                          \
        int j;                                                            \
        mask_t x0, x1, x2, x3, x4, x5;                                    \
        mask_t s;                                                         \
        mask_t n1_fwd;                                                    \
        simd_t n1_types;                                                  \
        simd_t n1_frames;                                                 \
        const mask_t n2_fwd    = simd_eq(simd_set1(strands[i]), ALL_FWD); \
        const simd_t n2_types  = simd_set1(types[i]);                     \
        const simd_t n2_frames = simd_set1(frames[i]);                    \
                                                                          \
        for (j = min; j < ((min + SIMD_MASK) & (~SIMD_MASK)); j++)        \
            skippable_generic_single(strands, types, frames, j, i, skip); \
        for (; j + SIMD_LANES < i; j += SIMD_LANES) {                     \
            n1_fwd     = simd_eq(simd_load(&strands[j]), ALL_FWD);        \
            n1_types   = simd_load(&types[j]);                            \
            n1_frames  = simd_load(&frames[j]);                           \
                                                                          \
            x0 =                 mask_eq(n1_fwd, n2_fwd);                 \
            x0 = mask_andnot(x0, simd_eq(n2_types, ALL_STOPS));           \
            x0 = mask_andnot(x0, simd_eq(n1_types, ALL_STOPS));           \
            s  =             x0;                                          \
                                                                          \
            x1 =                 n1_fwd;                                  \
            x1 = mask_andnot(x1, n2_fwd);                                 \
            x1 = mask_andnot(x1, simd_eq(n1_types, ALL_STOPS));           \
            s  = mask_or(    x1, s);                                      \
                                                                          \
            x2 =                 simd_eq(n1_types, ALL_STOPS);            \
            x2 = mask_andnot(x2, n1_fwd);                                 \
            x2 = mask_and(   x2, n2_fwd);                                 \
            s  = mask_or(    x2, s);                                      \
                                                                          \
            x3 =                 simd_eq(n2_types, ALL_STOPS);            \
            x3 = mask_andnot(x3, simd_eq(n1_types, ALL_STOPS));           \
            x3 = mask_andnot(x3, n1_fwd);                                 \
            x3 = mask_and(   x3, n2_fwd);                                 \
            s  = mask_or(    x3, s);                                      \
                                                                          \
            x4 =                 mask_eq(n1_fwd,     n2_fwd);             \
            x4 = mask_andnot(x4, simd_eq(n1_frames,  n2_frames));         \
            x4 = mask_and(   x4, n1_fwd)                 ;                \
            x4 = mask_andnot(x4, simd_eq(n1_types,   ALL_STOPS));         \
            x4 = mask_and(   x4, simd_eq(n2_types,   ALL_STOPS));         \
            s  = mask_or(    x4, s);                                      \
                                                                          \
            x5 =                 mask_eq(n1_fwd,     n2_fwd);             \
            x5 = mask_andnot(x5, simd_eq(n1_frames,  n2_frames));         \
            x5 = mask_andnot(x5, n1_fwd);                                 \
            x5 = mask_and(   x5, simd_eq(n1_types,   ALL_STOPS));         \
            x5 = mask_andnot(x5, simd_eq(n2_types,   ALL_STOPS));         \
            s  = mask_or(    x5, s);                                      \
                                                                          \
            simd_store(mask_convert(s), &skip[j]);                        \
        }                                                                 \
        for (; j < i; j++)                                                \
            skippable_generic_single(strands, types, frames, j, i, skip); \
                                                                          \
    } while(0)
