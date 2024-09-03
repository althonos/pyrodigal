#include "sequence.h"
#include "dprog.h"
#include "generic.h"

void skippable_generic(
    const int8_t*  restrict strands,
    const uint8_t* restrict types,
    const uint8_t* restrict frames,
    const int min,
    const int i,
          uint8_t* restrict skip
) {
    for (int j = min; j < i; j++)
      skippable_generic_single(strands, types, frames, j, i, skip);
}
