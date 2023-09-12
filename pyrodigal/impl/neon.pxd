from libc.stdint cimport int8_t, uint8_t

cdef extern from "impl/neon.h" nogil:
    void skippable_neon(const int8_t*, const uint8_t*, const uint8_t*, const int, const int, uint8_t*) noexcept
