cdef extern from "cpuinfo_aarch64.h" nogil:

    ctypedef struct Aarch64Features:
        int fp
        int asimd
        int evtstrm
        int aes
        int pmull
        int sha1
        int sha2
        int crc32
        int atomics
        int fphp
        int asimdhp
        int cpuid
        int asimdrdm
        int jscvt
        int fcma
        int lrcpc
        int dcpop
        int sha3
        int sm3
        int sm4
        int asimddp
        int sha512
        int sve
        int asimdfhm
        int dit
        int uscat
        int ilrcpc
        int flagm
        int ssbs
        int sb
        int paca
        int pacg
        int dcpodp
        int sve2
        int sveaes
        int svepmull
        int svebitperm
        int svesha3
        int svesm4
        int flagm2
        int frint
        int svei8mm
        int svef32mm
        int svef64mm
        int svebf16
        int i8mm
        int bf16
        int dgh
        int rng
        int bti

    ctypedef struct Aarch64Info:
        Aarch64Features features
        int             implementer
        int             variant
        int             part
        int             revision

    cdef Aarch64Info GetAarch64Info()
