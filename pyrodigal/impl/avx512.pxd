from libc.stdint cimport int8_t, uint8_t

cdef extern from "impl/avx512.h" nogil:
    void skippable_avx512(const int8_t*, const uint8_t*, const uint8_t*, const int, const int, uint8_t*) noexcept
