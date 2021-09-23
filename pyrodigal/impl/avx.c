#include "training.h"
#include "node.h"
#include "dprog.h"
#include "avx.h"

#ifdef __AVX2__

#include <immintrin.h>

void skippable_avx(
    const int8_t* strands,
    const uint8_t* types,
    const uint8_t* frames,

    const int min,
    const int i,
    uint8_t* skip
) {

    const __m256i all_stops  = _mm256_set1_epi8(STOP);
    const __m256i all_fwd    = _mm256_set1_epi8(1);
    const __m256i all_bwd    = _mm256_set1_epi8(-1);

    int j;
    __m256i x;
    __m256i s;
    __m256i n1_strands;
    __m256i n1_types;
    __m256i n1_frames;
    __m256i n2_strands = _mm256_set1_epi8(strands[i]);
    __m256i n2_types   = _mm256_set1_epi8(types[i]);
    __m256i n2_frames  = _mm256_set1_epi8(frames[i]);

    for (j = (min + 0x1F) & (~0x1F); j + 31 < i; j += 32) {
        n1_strands = _mm256_load_si256((__m256i*) &strands[j]);
        n1_types =   _mm256_load_si256((__m256i*) &types[j]);
        n1_frames =  _mm256_load_si256((__m256i*) &frames[j]);
        s = _mm256_set1_epi8(0);
        // 5'fwd->5'fwd
        // n1->strand == n2->strand && n2->type != STOP && n1->type != STOP
        x =                     _mm256_cmpeq_epi8(n1_strands, n2_strands);
        x = _mm256_andnot_si256(_mm256_cmpeq_epi8(n2_types, all_stops),    x);
        x = _mm256_andnot_si256(_mm256_cmpeq_epi8(n1_types, all_stops),    x);
        s = _mm256_or_si256(                                         s, x);
        // 5'fwd->5'ref, 5'fwd->3'rev
        // n2->strand == -1 && n1->strand == 1 && n1->type != STOP
        x =                     _mm256_cmpeq_epi8(n2_strands, all_bwd);
        x = _mm256_and_si256(   _mm256_cmpeq_epi8(n1_strands, all_fwd), x);
        x = _mm256_andnot_si256(_mm256_cmpeq_epi8(n1_types, all_stops), x);
        s = _mm256_or_si256(                                         s, x);
        // 5'fwd
        // n1->type == STOP && n1->strand == -1 && n2->strand == -1
        x =                  _mm256_cmpeq_epi8(n1_types, all_stops);
        x = _mm256_and_si256(_mm256_cmpeq_epi8(n1_strands, all_bwd), x);
        x = _mm256_and_si256(_mm256_cmpeq_epi8(n2_strands, all_fwd), x);
        s = _mm256_or_si256(                                      s, x);
        // 5'rev->3'fwd
        // n2->type == STOP && n1->strand == -1 && n2->strand == 1 && n1->type != STOP
        x =                     _mm256_cmpeq_epi8(n2_types, all_stops);
        x = _mm256_and_si256(   _mm256_cmpeq_epi8(n1_strands, all_bwd), x);
        x = _mm256_and_si256(   _mm256_cmpeq_epi8(n2_strands, all_fwd), x);
        x = _mm256_andnot_si256(_mm256_cmpeq_epi8(n1_types, all_stops), x);
        s = _mm256_or_si256(                                         s, x);
        // 5'fwd->3'fwd
        // n1->strand == n2->strand && n1->strand == 1 && n1->type != STOP && n2->type == STOP && n1->ndx%3 != n2->ndx%3
        x =                     _mm256_cmpeq_epi8(n1_strands, n2_strands);
        x = _mm256_and_si256(   _mm256_cmpeq_epi8(n1_strands, all_fwd),    x);
        x = _mm256_andnot_si256(_mm256_cmpeq_epi8(n1_types,   all_stops),  x);
        x = _mm256_and_si256(   _mm256_cmpeq_epi8(n2_types,   all_stops),  x);
        x = _mm256_andnot_si256(_mm256_cmpeq_epi8(n1_frames,  n2_frames),  x);
        s = _mm256_or_si256(                                            s, x);
        // 3'rev->5'rev
        // n1->strand == n2->strand && n1->strand == -1 && n1->type == STOP && n2->type != STOP && n1->ndx%3 != n2->ndx%3
        x =                     _mm256_cmpeq_epi8(n1_strands, n2_strands);
        x = _mm256_and_si256(   _mm256_cmpeq_epi8(n1_strands, all_bwd),    x);
        x = _mm256_and_si256(   _mm256_cmpeq_epi8(n1_types,   all_stops),  x);
        x = _mm256_andnot_si256(_mm256_cmpeq_epi8(n2_types,   all_stops),  x);
        x = _mm256_andnot_si256(_mm256_cmpeq_epi8(n1_frames,  n2_frames),  x);
        s = _mm256_or_si256(                                            s, x);
        // store result mask
        _mm256_store_si256((__m256i*) &skip[j], s);
    }
}
#endif
