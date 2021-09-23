cdef extern from "cpuinfo_arm.h" nogil:

    ctypedef struct ArmFeatures:
        int swp
        int half
        int thumb
        int _26bit
        int fastmult
        int fpa
        int vfp
        int edsp
        int java
        int iwmmxt
        int crunch
        int thumbee
        int neon
        int vfpv3
        int vfpv3d16
        int tls
        int vfpv4
        int idiva
        int idivt
        int vfpd32
        int lpae
        int evtstrm
        int aes
        int pmull
        int sha1
        int sha2
        int crc32

    ctypedef struct ArmInfo:
        ArmFeatures features
        int         implementer
        int         architecture
        int         variant
        int         part
        int         revision

    cdef ArmInfo GetArmInfo()
