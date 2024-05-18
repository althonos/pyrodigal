from libc.stdint cimport uint8_t

cdef extern from "_sequence.h" nogil:

    cdef enum nucleotide:
        A = 0b000
        G = 0b001
        C = 0b010
        T = 0b011
        N = 0b110

    const uint8_t _complement[N+1]
    const char    _letters   [N+1]

    cdef bint _is_a(const uint8_t* digits, const int slen, const int i, const int strand) noexcept
    cdef bint _is_g(const uint8_t* digits, const int slen, const int i, const int strand) noexcept
    cdef bint _is_c(const uint8_t* digits, const int slen, const int i, const int strand) noexcept
    cdef bint _is_t(const uint8_t* digits, const int slen, const int i, const int strand) noexcept
    cdef bint _is_gc(const uint8_t* digits, const int slen, const int i, const int strand) noexcept

    cdef bint _is_start(const uint8_t* digits, const int slen, const int i, const int tt, const int strand) noexcept
    cdef bint _is_stop(const uint8_t* digits, const int slen, const int i, const int tt, const int strand) noexcept

    cdef bint _is_atg(const uint8_t* digits, const int slen, const int i, const int strand) noexcept
    cdef bint _is_gtg(const uint8_t* digits, const int slen, const int i, const int strand) noexcept
    cdef bint _is_ttg(const uint8_t* digits, const int slen, const int i, const int strand) noexcept

    cdef int _mer_ndx(const uint8_t* digits, const int slen, const int i, const int length, const int strand) noexcept

    cdef char _amino(const uint8_t* digits, const int slen, const int i, const int tt, const int strand, const bint strict) noexcept