from pyrodigal.prodigal.bitmap cimport bitmap_t
from pyrodigal.prodigal.training cimport _training


cdef extern from "sequence.h":

    cdef size_t MAX_MASKS = 5000
    cdef size_t MASK_SIZE = 50

    cdef struct _mask:
        int begin
        int end

    # given a bitmap_t sequence `seq` of length `len`, write the reverse
    # complement of `seq` into `rseq`, ignoring chars set a ``
    cdef void rcom_seq(const bitmap_t seq, bitmap_t rseq, bitmap_t useq, int len)

    cdef char is_a(bitmap_t seq, int n) nogil
    cdef char is_c(bitmap_t seq, int n) nogil
    cdef char is_g(bitmap_t seq, int n) nogil
    cdef char is_t(bitmap_t seq, int n) nogil
    cdef char is_gc(bitmap_t seq, int n) nogil

    cdef char amino(bitmap_t seq, int n, _training* tinf, bint is_init) nogil
    cdef int amino_num(char) nogil
    cdef char amino_letter(int) nogil
