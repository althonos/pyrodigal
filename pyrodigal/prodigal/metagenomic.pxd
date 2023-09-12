from pyrodigal.prodigal.training cimport _training


cdef extern from "metagenomic.h" nogil:

    cdef const ssize_t NUM_BIN
    cdef const ssize_t NUM_META
    cdef const ssize_t SAMPLE_LEN
    cdef const ssize_t MAX_SAMPLE

    cdef struct _metagenomic_bin:
        char desc[500]
        _training* tinf
        
        # NOTE(@althonos): The attributes below are unused so let's just mask
        #                  them to prevent our code from fiddling with them.
        # int index
        # int clusnum
        # double weight
        # double gc

    void initialize_metagenomic_bins(_metagenomic_bin*) noexcept
