#ifndef _PYRODIGAL_SEQUENCE_H
#define _PYRODIGAL_SEQUENCE_H
#include <stdint.h>

enum nucleotide {
    A = 0b000,
    G = 0b001,
    C = 0b010,
    T = 0b011,
    N = 0b110
};

const uint8_t _complement [N+1] = {   T,   C,   G,   A,   N,   N,   N };
const char    _letters    [N+1] = { 'A', 'G', 'C', 'T', 'N', 'N', 'N' };

static inline int _is_a(const uint8_t* digits, const int slen, const int i, const int strand) {
    return (strand == 1) ? digits[i] == A : digits[slen - 1 - i] == T;
}

static inline int _is_g(const uint8_t* digits, const int slen, const int i, const int strand) {
    return (strand == 1) ? digits[i] == G : digits[slen - 1 - i] == C;
}

static inline int _is_c(const uint8_t* digits, const int slen, const int i, const int strand) {
    return (strand == 1) ? digits[i] == C : digits[slen - 1 - i] == G;
}

static inline int _is_t(const uint8_t* digits, const int slen, const int i, const int strand) {
    return (strand == 1) ? digits[i] == T : digits[slen - 1 - i] == A;
}

static inline int _is_gc(const uint8_t* digits, const int slen, const int i, const int strand) {
    // NB(@althonos): In the original Prodigal implementation, any unknown
    //                character gets encoded as a C, so it gets counted
    //                when computing the most GC frame. We reproduce this
    //                behaviour here, but a better solution would be to
    //                count only known letters.
    uint_fast8_t nuc = (strand == 1) ? digits[i] : digits[slen - 1 - i];
    return (nuc != A) && (nuc != T);
}

static inline int _is_start(const uint8_t* digits, const int slen, const int i, const int tt, const int strand) {
    uint_fast8_t x0, x1, x2;

    if (strand == 1) {
        x0 = digits[i];
        x1 = digits[i+1];
        x2 = digits[i+2];
    } else {
        x0 = digits[slen - 1 - i] ^ 0b11;
        x1 = digits[slen - 2 - i] ^ 0b11;
        x2 = digits[slen - 3 - i] ^ 0b11;
    }

    // ATG
    if ((x0 == A) && (x1 == T) && (x2 == G))
        return 1;
    // Codes that only use ATG
    if ((tt == 6) || (tt == 10) || (tt == 14) || (tt == 15) || (tt == 16) || (tt == 2))
        return 0;
    // GTG
    if ((x0 == G) && (x1 == T) && (x2 == G))
        return !((tt == 1) || (tt == 3) || (tt == 12) || (tt == 2));
    // TTG
    if ((x0 == T) && (x1 == T) && (x2 == G))
        return !((tt < 4) || (tt == 9) || ((tt >= 21) && (tt < 25)));

    // other codons
    return 0;
}

static inline int _is_stop(const uint8_t* digits, const int slen, const int i, const int tt, const int strand) {
    uint_fast8_t x0, x1, x2;

    if (strand == 1) {
        x0 = digits[i];
        x1 = digits[i+1];
        x2 = digits[i+2];
    } else {
        x0 = digits[slen - 1 - i] ^ 0b11;
        x1 = digits[slen - 2 - i] ^ 0b11;
        x2 = digits[slen - 3 - i] ^ 0b11;
    }

    // TAG
    if ((x0 == T) && (x1 == A) && (x2 == G))
        return !((tt == 6) || (tt == 15) || (tt == 16) || (tt == 22));
    // TGA
    if ((x0 == T) && (x1 == G) && (x2 == A))
        return !(
                (tt ==  2) || (tt ==  3) || (tt ==  4) || (tt ==  5)
             || (tt ==  9) || (tt == 10) || (tt == 13) || (tt == 14)
             || (tt == 21) || (tt == 25)
        );
    // TAA
    if ((x0 == T) && (x1 == A) && (x2 == A))
        return !((tt == 6) || (tt == 14));

    // Code 2: AGA / AGG
    if (tt == 2)
        return (x0 == A) && (x1 == G) && ((x2 == A) || (x2 == G));
    // Code 22: TCA
    if (tt == 22)
        return (x0 == T) && (x1 == C) && (x2 == A);
    // Code 23: TTA
    if (tt == 23)
        return (x0 == T) && (x1 == T) && (x2 == A);

    // other codons
    return 0;
}

static inline int _is_atg(const uint8_t* digits, const int slen, const int i, const int strand) {
    uint_fast8_t x0, x1, x2;

    if (strand == 1) {
        x0 = digits[i];
        x1 = digits[i+1];
        x2 = digits[i+2];
    } else {
        x0 = digits[slen - 1 - i] ^ 0b11;
        x1 = digits[slen - 2 - i] ^ 0b11;
        x2 = digits[slen - 3 - i] ^ 0b11;
    }

    return (x0 == A) && (x1 == T) && (x2 == G);
}

static inline int _is_gtg(const uint8_t* digits, const int slen, const int i, const int strand) {
    uint_fast8_t x0, x1, x2;

    if (strand == 1) {
        x0 = digits[i];
        x1 = digits[i+1];
        x2 = digits[i+2];
    } else {
        x0 = digits[slen - 1 - i] ^ 0b11;
        x1 = digits[slen - 2 - i] ^ 0b11;
        x2 = digits[slen - 3 - i] ^ 0b11;
    }

    return (x0 == G) && (x1 == T) && (x2 == G);
}

static inline int _is_ttg(const uint8_t* digits, const int slen, const int i, const int strand) {
    uint_fast8_t x0, x1, x2;

    if (strand == 1) {
        x0 = digits[i];
        x1 = digits[i+1];
        x2 = digits[i+2];
    } else {
        x0 = digits[slen - 1 - i] ^ 0b11;
        x1 = digits[slen - 2 - i] ^ 0b11;
        x2 = digits[slen - 3 - i] ^ 0b11;
    }

    return (x0 == T) && (x1 == T) && (x2 == G);
}

static int _mer_ndx(const uint8_t* digits, const int slen, const int i, const int length, const int strand) {
    int j, k;
    int ndx = 0;
    if (strand == 1) {
        for (j = 0, k = i; j < length; k++, j++) {
            ndx |= (digits[k] & 0b11) << 2*j;
        }
    } else {
        for (j = 0, k = slen - 1 - i; j < length; k--, j++) {
            ndx |= (_complement[digits[k]] & 0b11) << 2*j;
        }
    }
  return ndx;
}
#endif
