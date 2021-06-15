cdef extern from "_utils.h" nogil:

    cdef struct _mini_training:
        double gc
        int trans_table
