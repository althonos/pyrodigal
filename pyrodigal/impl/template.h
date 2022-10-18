#define skippable_simd(strands, types, frames, min, i, skip) do {         \
        const simd_t ALL_STOPS = simd_set1(STOP);                         \
        const simd_t ALL_FWD   = simd_set1(1);                            \
        const simd_t ALL_BWD   = simd_set1(-1);                           \
                                                                          \
        int j;                                                            \
        simd_t x0, x1, x2, x3, x4, x5;                                    \
        simd_t s;                                                         \
        simd_t n1_strands;                                                \
        simd_t n1_types;                                                  \
        simd_t n1_frames;                                                 \
        const simd_t n2_strands = simd_set1(strands[i]);                  \
        const simd_t n2_types   = simd_set1(types[i]);                    \
        const simd_t n2_frames  = simd_set1(frames[i]);                   \
                                                                          \
        for (j = min; j < ((min + SIMD_MASK) & (~SIMD_MASK)); j++)        \
            skippable_generic_single(strands, types, frames, j, i, skip); \
        for (; j + SIMD_LANES < i; j += SIMD_LANES) {                     \
            n1_strands = simd_load(&strands[j]);                          \
            n1_types   = simd_load(&types[j]);                            \
            n1_frames  = simd_load(&frames[j]);                           \
                                                                          \
            x0 =                 simd_eq(n1_strands, n2_strands);         \
            x0 = simd_andnot(x0, simd_eq(n2_types, ALL_STOPS));           \
            x0 = simd_andnot(x0, simd_eq(n1_types, ALL_STOPS));           \
            s  =             x0;                                          \
                                                                          \
            x1 =                 simd_eq(n2_strands, ALL_BWD);            \
            x1 = simd_and(   x1, simd_eq(n1_strands, ALL_FWD));           \
            x1 = simd_andnot(x1, simd_eq(n1_types, ALL_STOPS));           \
            s  = simd_or(    x1, s);                                      \
                                                                          \
            x2 =              simd_eq(n1_types, ALL_STOPS);               \
            x2 = simd_and(x2, simd_eq(n1_strands, ALL_BWD));              \
            x2 = simd_and(x2, simd_eq(n2_strands, ALL_FWD));              \
            s  = simd_or( x2, s);                                         \
                                                                          \
            x3 =                 simd_eq(n2_types, ALL_STOPS);            \
            x3 = simd_and(   x3, simd_eq(n1_strands, ALL_BWD));           \
            x3 = simd_and(   x3, simd_eq(n2_strands, ALL_FWD));           \
            x3 = simd_andnot(x3, simd_eq(n1_types, ALL_STOPS));           \
            s  = simd_or(    x3, s);                                      \
                                                                          \
            x4 =                 simd_eq(n1_strands, n2_strands);         \
            x4 = simd_and(   x4, simd_eq(n1_strands, ALL_FWD));           \
            x4 = simd_andnot(x4, simd_eq(n1_types,   ALL_STOPS));         \
            x4 = simd_and(   x4, simd_eq(n2_types,   ALL_STOPS));         \
            x4 = simd_andnot(x4, simd_eq(n1_frames,  n2_frames));         \
            s = simd_or(     x4, s);                                      \
                                                                          \
            x5 =                 simd_eq(n1_strands, n2_strands);         \
            x5 = simd_and(   x5, simd_eq(n1_strands, ALL_BWD));           \
            x5 = simd_and(   x5, simd_eq(n1_types,   ALL_STOPS));         \
            x5 = simd_andnot(x5, simd_eq(n2_types,   ALL_STOPS));         \
            x5 = simd_andnot(x5, simd_eq(n1_frames,  n2_frames));         \
            s  = simd_or(    x5, s);                                      \
                                                                          \
            simd_store(s, &skip[j]);                                      \
        }                                                                 \
        for (; j < i; j++)                                                \
            skippable_generic_single(strands, types, frames, j, i, skip); \
                                                                          \
    } while(0)
