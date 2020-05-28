from pyrodigal.prodigal.bitmap cimport bitmap_t
from pyrodigal.prodigal.training cimport _training


cdef extern from "sequence.h" nogil:

    size_t MAX_MASKS = 5000
    size_t MASK_SIZE = 50

    struct _mask:
        int begin
        int end

    enum:
        ATG = 0
        GTG = 1
        TTG = 2
        STOP = 3

    # given a bitmap_t sequence `seq` of length `len`, write the reverse
    # complement of `seq` into `rseq`, ignoring chars set a ``
    void rcom_seq(const bitmap_t seq, bitmap_t rseq, bitmap_t useq, int slen)

    char is_a(bitmap_t seq, int n)
    char is_c(bitmap_t seq, int n)
    char is_g(bitmap_t seq, int n)
    char is_t(bitmap_t seq, int n)
    char is_gc(bitmap_t seq, int n)

    double gc_content(bitmap_t seq, int a, int b)

    char amino(bitmap_t seq, int n, _training* tinf, bint is_init)
    int amino_num(char)
    char amino_letter(int)

    int* calc_most_gc_frame(bitmap_t seq, int slen)
