# coding: utf-8
# cython: language_level=3, linetrace=True, binding=True

from libc.stdint cimport int8_t, uint8_t

from ..lib cimport BaseConnectionScorer

cdef extern from "swar64.h" nogil:
    void skippable_swar64(const uint8_t*, const uint8_t*, const uint8_t*, const int, const int, uint8_t*) noexcept

cdef class SWAR64ConnectionScorer(BaseConnectionScorer):
    def __cinit__(self):
        self.skippable = skippable_swar64
        self.enabled = True
