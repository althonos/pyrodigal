#include "training.h"
#include "node.h"
#include "dprog.h"
#include "neon.h"

#ifdef __ARM_NEON__

#include <arm_neon.h>

void skippable_neon(
    const int8_t* strands,
    const uint8_t* types,
    const uint8_t* frames,

    const int min,
    const int i,
    uint8_t* skip
) {

  const uint8x16_t all_stops  = vdupq_n_u8(STOP);
  const uint8x16_t all_fwd    = vdupq_n_u8(1);
  const uint8x16_t all_bwd    = vdupq_n_u8(-1);

  int j;
  uint8x16_t x;
  uint8x16_t s;
  uint8x16_t n1_strands;
  uint8x16_t n1_types;
  uint8x16_t n1_frames;
  uint8x16_t n2_strands = vdupq_n_u8(strands[i]);
  uint8x16_t n2_types   = vdupq_n_u8(types[i]);
  uint8x16_t n2_frames  = vdupq_n_u8(frames[i]);

  for (j = (min + 0xF) & (~0xF); j + 15 < i; j += 16) {
      n1_strands = vld1q_u8((uint8_t*) &strands[j]);
      n1_types   = vld1q_u8(&types[j]);
      n1_frames  = vld1q_u8(&frames[j]);
      s          = vdupq_n_u8(0);
      // 5'fwd->5'fwd
      // n1->strand == n2->strand && n2->type != STOP && n1->type != STOP
      x =                  vceqq_u8(n1_strands, n2_strands);
      x = vandq_u8(vmvnq_u8(vceqq_u8(n2_types, all_stops)),    x);
      x = vandq_u8(vmvnq_u8(vceqq_u8(n1_types, all_stops)),    x);
      s = vorrq_u8(                                         s, x);
      // 5'fwd->5'ref, 5'fwd->3'rev
      // n2->strand == -1 && n1->strand == 1 && n1->type != STOP
      x =                  vceqq_u8(n2_strands, all_bwd);
      x = vandq_u8(   vceqq_u8(n1_strands, all_fwd), x);
      x = vandq_u8(vmvnq_u8(vceqq_u8(n1_types, all_stops)), x);
      s = vorrq_u8(                                      s, x);
      // 5'fwd
      // n1->type == STOP && n1->strand == -1 && n2->strand == -1
      x =               vceqq_u8(n1_types, all_stops);
      x = vandq_u8(vceqq_u8(n1_strands, all_bwd), x);
      x = vandq_u8(vceqq_u8(n2_strands, all_fwd), x);
      s = vorrq_u8(                                   s, x);
      // 5'rev->3'fwd
      // n2->type == STOP && n1->strand == -1 && n2->strand == 1 && n1->type != STOP
      x =                  vceqq_u8(n2_types, all_stops);
      x = vandq_u8(   vceqq_u8(n1_strands, all_bwd), x);
      x = vandq_u8(   vceqq_u8(n2_strands, all_fwd), x);
      x = vandq_u8(vmvnq_u8(vceqq_u8(n1_types, all_stops)), x);
      s = vorrq_u8(                                      s, x);
      // 5'fwd->3'fwd
      // n1->strand == n2->strand && n1->strand == 1 && n1->type != STOP && n2->type == STOP && n1->ndx%3 != n2->ndx%3
      x =                  vceqq_u8(n1_strands, n2_strands);
      x = vandq_u8(   vceqq_u8(n1_strands, all_fwd),    x);
      x = vandq_u8(vmvnq_u8(vceqq_u8(n1_types,   all_stops)),  x);
      x = vandq_u8(   vceqq_u8(n2_types,   all_stops),  x);
      x = vandq_u8(vmvnq_u8(vceqq_u8(n1_frames,  n2_frames)),  x);
      s = vorrq_u8(                                        s,  x);
      // 3'rev->5'rev
      // n1->strand == n2->strand && n1->strand == -1 && n1->type == STOP && n2->type != STOP && n1->ndx%3 != n2->ndx%3
      x =                   vceqq_u8(n1_strands, n2_strands);
      x = vandq_u8(         vceqq_u8(n1_strands, all_bwd),    x);
      x = vandq_u8(         vceqq_u8(n1_types,   all_stops),  x);
      x = vandq_u8(vmvnq_u8(vceqq_u8(n2_types,   all_stops)), x);
      x = vandq_u8(vmvnq_u8(vceqq_u8(n1_frames,  n2_frames)), x);
      s = vorrq_u8(                                        s, x);

      // store result mask
      vst1q_u8(&skip[j], s);
  }
}
#endif
