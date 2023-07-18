from pyrodigal.prodigal.training cimport _training


cdef extern from "metagenomic.h" nogil:

    cdef const ssize_t NUM_BIN
    cdef const ssize_t NUM_META
    cdef const ssize_t SAMPLE_LEN
    cdef const ssize_t MAX_SAMPLE

    cdef struct _metagenomic_bin:
        int index
        int clusnum
        char desc[500]
        double weight
        double gc
        _training* tinf

    void initialize_metagenomic_bins(_metagenomic_bin*)
