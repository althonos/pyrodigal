# coding: utf-8
# cython: language_level=3, linetrace=True, binding=True

from libc.stdint cimport int8_t, uint8_t

from ..lib cimport BaseConnectionScorer

cdef extern from "neon.h" nogil:
    void skippable_neon(const int8_t*, const uint8_t*, const uint8_t*, const int, const int, uint8_t*) noexcept

cdef class NEONConnectionScorer(BaseConnectionScorer):
    def __cinit__(self):
        self.skippable = skippable_neon
        self.enabled = True
