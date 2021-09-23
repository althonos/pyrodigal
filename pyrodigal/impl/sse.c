#include "training.h"
#include "node.h"
#include "dprog.h"
#include "sse.h"

#ifdef __SSE2__

#include <emmintrin.h>

void skippable_sse(
    const int8_t* strands,
    const uint8_t* types,
    const uint8_t* frames,

    const int min,
    const int i,
    uint8_t* skip
) {

  const __m128i all_stops  = _mm_set1_epi8(STOP);
  const __m128i all_fwd    = _mm_set1_epi8(1);
  const __m128i all_bwd    = _mm_set1_epi8(-1);

  int j;
  __m128i x;
  __m128i s;
  __m128i n1_strands;
  __m128i n1_types;
  __m128i n1_frames;
  __m128i n2_strands = _mm_set1_epi8(strands[i]);
  __m128i n2_types   = _mm_set1_epi8(types[i]);
  __m128i n2_frames  = _mm_set1_epi8(frames[i]);

  for (j = (min + 0xF) & (~0xF); j + 15 < i; j += 16) {
      n1_strands = _mm_load_si128((__m128i*) &strands[j]);
      n1_types =   _mm_load_si128((__m128i*) &types[j]);
      n1_frames =  _mm_load_si128((__m128i*) &frames[j]);
      s = _mm_xor_si128(s, s);
      // 5'fwd->5'fwd
      // n1->strand == n2->strand && n2->type != STOP && n1->type != STOP
      x =                  _mm_cmpeq_epi8(n1_strands, n2_strands);
      x = _mm_andnot_si128(_mm_cmpeq_epi8(n2_types, all_stops),    x);
      x = _mm_andnot_si128(_mm_cmpeq_epi8(n1_types, all_stops),    x);
      s = _mm_or_si128(                                         s, x);
      // 5'fwd->5'ref, 5'fwd->3'rev
      // n2->strand == -1 && n1->strand == 1 && n1->type != STOP
      x =                  _mm_cmpeq_epi8(n2_strands, all_bwd);
      x = _mm_and_si128(   _mm_cmpeq_epi8(n1_strands, all_fwd), x);
      x = _mm_andnot_si128(_mm_cmpeq_epi8(n1_types, all_stops), x);
      s = _mm_or_si128(                                      s, x);
      // 5'fwd
      // n1->type == STOP && n1->strand == -1 && n2->strand == -1
      x =               _mm_cmpeq_epi8(n1_types, all_stops);
      x = _mm_and_si128(_mm_cmpeq_epi8(n1_strands, all_bwd), x);
      x = _mm_and_si128(_mm_cmpeq_epi8(n2_strands, all_fwd), x);
      s = _mm_or_si128(                                   s, x);
      // 5'rev->3'fwd
      // n2->type == STOP && n1->strand == -1 && n2->strand == 1 && n1->type != STOP
      x =                  _mm_cmpeq_epi8(n2_types, all_stops);
      x = _mm_and_si128(   _mm_cmpeq_epi8(n1_strands, all_bwd), x);
      x = _mm_and_si128(   _mm_cmpeq_epi8(n2_strands, all_fwd), x);
      x = _mm_andnot_si128(_mm_cmpeq_epi8(n1_types, all_stops), x);
      s = _mm_or_si128(                                      s, x);
      // 5'fwd->3'fwd
      // n1->strand == n2->strand && n1->strand == 1 && n1->type != STOP && n2->type == STOP && n1->ndx%3 != n2->ndx%3
      x =                  _mm_cmpeq_epi8(n1_strands, n2_strands);
      x = _mm_and_si128(   _mm_cmpeq_epi8(n1_strands, all_fwd),    x);
      x = _mm_andnot_si128(_mm_cmpeq_epi8(n1_types,   all_stops),  x);
      x = _mm_and_si128(   _mm_cmpeq_epi8(n2_types,   all_stops),  x);
      x = _mm_andnot_si128(_mm_cmpeq_epi8(n1_frames,  n2_frames),  x);
      s = _mm_or_si128(                                        s,  x);
      // 3'rev->5'rev
      // n1->strand == n2->strand && n1->strand == -1 && n1->type == STOP && n2->type != STOP && n1->ndx%3 != n2->ndx%3
      x =                  _mm_cmpeq_epi8(n1_strands, n2_strands);
      x = _mm_and_si128(   _mm_cmpeq_epi8(n1_strands, all_bwd),    x);
      x = _mm_and_si128(   _mm_cmpeq_epi8(n1_types,   all_stops),  x);
      x = _mm_andnot_si128(_mm_cmpeq_epi8(n2_types,   all_stops),  x);
      x = _mm_andnot_si128(_mm_cmpeq_epi8(n1_frames,  n2_frames),  x);
      s = _mm_or_si128(                                        s,  x);

      // store result mask
      _mm_store_si128((__m128i*) &skip[j], s);
  }
}
#endif
