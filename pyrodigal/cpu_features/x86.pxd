cdef extern from "cpuinfo_x86.h" nogil:

    ctypedef struct X86Features:
        int mmx

        int sse
        int sse2
        int sse3
        int ssse3
        int sse4_1
        int sse4_2
        int sse4a

        int avx
        int avx2

        int avx512f
        int avx512cd

    ctypedef struct X86Info:
        X86Features features
        int         family
        int         model
        int         stepping
        char        vendor[13]

    cdef X86Info GetX86Info()
