#ifndef _PYRODIGAL_IMPL_GENERIC_H
#define _PYRODIGAL_IMPL_GENERIC_H

#include <stdint.h>

void skippable_generic(const int8_t*, const uint8_t*, const uint8_t*, const int, const int, uint8_t*);
static inline void skippable_generic_single(
    const int8_t* strands,
    const uint8_t* types,
    const uint8_t* frames,
    const int j,
    const int i,
    uint8_t* skip
) {
    int_fast8_t  n1_strand = strands[j];
    uint_fast8_t n1_frame  = frames[j];
    uint_fast8_t n1_type   = types[j];

    int_fast8_t  n2_strand = strands[i];
    uint_fast8_t n2_frame  = frames[i];
    uint_fast8_t n2_type   = types[i];

    skip[j]  = n1_type   != STOP && n2_type != STOP && n1_strand == n2_strand;
    skip[j] |= n1_strand ==  1   && n1_type != STOP && n2_strand == -1;
    skip[j] |= n1_strand == -1   && n1_type == STOP && n2_strand ==  1;
    skip[j] |= n1_strand == -1   && n1_type != STOP && n2_strand ==  1 && n2_type == STOP;
    skip[j] |= n1_strand == n2_strand && n1_strand == 1  && n1_type != STOP && n2_type == STOP && n1_frame != n2_frame;
    skip[j] |= n1_strand == n2_strand && n1_strand == -1 && n1_type == STOP && n2_type != STOP && n1_frame != n2_frame;
}

#endif
