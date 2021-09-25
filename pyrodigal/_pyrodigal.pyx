# coding: utf-8
# cython: language_level=3, linetrace=True

"""Bindings to Prodigal, an ORF finder for genomes and metagenomes.

Example:
    Pyrodigal can work on any DNA sequence stored in either a text or a byte
    array. To load a sequence from one of the common sequence formats, you can
    for instance use `Biopython <https://github.com/biopython/biopython>`_::

        >>> import gzip
        >>> import Bio.SeqIO
        >>> with gzip.open("pyrodigal/tests/data/KK037166.fna.gz") as f:
        ...     record = Bio.SeqIO.read(f, "fasta")

    Then use Pyrodigal to find the genes in *metagenomic* mode (without
    training first), and then build a map of codon frequencies for each
    gene::

        >>> from collections import Counter
        >>> import pyrodigal
        >>> p = pyrodigal.Pyrodigal(meta=True)
        >>> for prediction in p.find_genes(record.seq.encode()):
        ...     gene_seq = prediction.sequence()
        ...     codon_counter = Counter()
        ...     for i in range(len(gene_seq), 3):
        ...         codon_counter[gene_seq[i:i+3]] += 1
        ...     codon_frequencies = {
        ...         codon:count/(len(gene_seq)//3)
        ...         for codon, count in codon_counter.items()
        ...     }

Caution:
    In Pyrodigal, sequences are assumed to contain only the usual nucleotides
    (A/T/G/C) as lowercase or uppercase letters; any other symbol will be
    treated as an unknown nucleotide. Be careful to remove the gap characters
    if loading sequences from a multiple alignment file.

See Also:
    The `academic paper for Prodigal <https://doi.org/10.1186/1471-2105-11-119>_`
    which describes the algorithm in use.

"""

# ----------------------------------------------------------------------------

from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AsString
from cpython.exc cimport PyErr_CheckSignals
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from cpython.ref cimport Py_INCREF
from cpython.tuple cimport PyTuple_New, PyTuple_SET_ITEM
from libc.math cimport sqrt, log, pow, fmax, fmin
from libc.stdint cimport int8_t, uint8_t, uintptr_t
from libc.stdio cimport printf
from libc.stdlib cimport abs, malloc, calloc, free, qsort
from libc.string cimport memcpy, memchr, memset, strstr

from pyrodigal.prodigal cimport bitmap, dprog, gene, node, sequence
from pyrodigal.prodigal.bitmap cimport bitmap_t
from pyrodigal.prodigal.gene cimport _gene
from pyrodigal.prodigal.metagenomic cimport NUM_META, _metagenomic_bin, initialize_metagenomic_bins
from pyrodigal.prodigal.node cimport _motif, _node, MIN_EDGE_GENE, MIN_GENE, cross_mask
from pyrodigal.prodigal.sequence cimport _mask, node_type, rcom_seq
from pyrodigal.prodigal.training cimport _training
from pyrodigal._unicode cimport *

IF TARGET_CPU == "x86":
    from pyrodigal.cpu_features.x86 cimport GetX86Info, X86Info
    IF SSE2_BUILD_SUPPORT:
        from pyrodigal.impl.sse cimport skippable_sse
    IF AVX2_BUILD_SUPPORT:
        from pyrodigal.impl.avx cimport skippable_avx
ELIF TARGET_CPU == "arm" or TARGET_CPU == "aarch64":
    IF TARGET_CPU == "arm":
        from pyrodigal.cpu_features.arm cimport GetArmInfo, ArmInfo
    IF NEON_BUILD_SUPPORT:
        from pyrodigal.impl.neon cimport skippable_neon

# ----------------------------------------------------------------------------

import warnings
import threading

# --- Module-level constants -------------------------------------------------

cdef int    IDEAL_SINGLE_GENOME = 100000
cdef int    MIN_SINGLE_GENOME   = 20000
cdef int    WINDOW              = 120
cdef size_t MIN_GENES_ALLOC     = 8
cdef size_t MIN_NODES_ALLOC     = 8 * MIN_GENES_ALLOC
cdef set    TRANSLATION_TABLES  = set(range(1, 7)) | set(range(9, 17)) | set(range(21, 26))

# --- Input sequence ---------------------------------------------------------

cdef enum:
    A = 0b000
    G = 0b001
    C = 0b010
    T = 0b011
    N = 0b110  # use 6 so that N & 0x3 == C

cdef uint8_t _translation[N+1]
for i in range(N+1):
    _translation[i] = N
_translation[A] = T
_translation[T] = A
_translation[C] = G
_translation[G] = C

cdef Py_UCS4 _letters[N+1]
for i in range(N+1):
    _letters[i] = "N"
_letters[A] = "A"
_letters[T] = "T"
_letters[C] = "C"
_letters[G] = "G"


cdef class Sequence:
    """A digitized input sequence.
    """

    def __cinit__(self):
        self.slen = 0
        self.gc = 0.0
        self.digits = NULL

    def __dealloc__(self):
        PyMem_Free(self.digits)

    def __len__(self):
        return self.slen

    def __sizeof__(self):
        return self.slen * sizeof(uint8_t) + sizeof(self)

    cdef int _allocate(self, int slen) except 1:
        self.slen = slen
        self.digits = <uint8_t*> PyMem_Malloc(slen * sizeof(uint8_t))
        if self.digits == NULL:
            raise MemoryError()
        with nogil:
            memset(self.digits, 0, slen * sizeof(uint8_t))
        return 0

    @classmethod
    def from_bytes(cls, const unsigned char[:] sequence):
        """from_bytes(cls, sequence)\n--

        Create a new `Sequence` object from an ASCII-encoded sequence.

        Arguments:
            sequence (`bytes`): The ASCII-encoded sequence to use. Any object
                implementing the *buffer protocol* is supported.

        """
        cdef int           i
        cdef int           j
        cdef unsigned char letter
        cdef Sequence      seq
        cdef int           gc_count = 0

        seq = Sequence.__new__(Sequence)
        seq._allocate(sequence.shape[0])

        with nogil:
            for i, j in enumerate(range(0, seq.slen * 2, 2)):
                letter = sequence[i]
                if letter == b'A' or letter == b'a':
                    seq.digits[i] = A
                elif letter == b'T' or letter == b't':
                    seq.digits[i] = T
                elif letter == b'G' or letter == b'g':
                    seq.digits[i] = G
                    gc_count += 1
                elif letter == b'C' or letter == b'c':
                    seq.digits[i] = C
                    gc_count += 1
                else:
                    seq.digits[i] = N
            if seq.slen > 0:
                seq.gc = (<double> gc_count) / (<double> seq.slen)

        return seq

    @classmethod
    def from_string(cls, str sequence):
        """from_string(cls, sequence)\n--

        Create a new `Sequence` object from a Unicode sequence.

        Arguments:
            sequence (`str`): The Unicode sequence to use.

        """
        cdef int      i
        cdef int      j
        cdef Py_UCS4  letter
        cdef Sequence seq
        cdef int      kind
        cdef void*    data
        cdef int      gc_count = 0

        # make sure the unicode string is in canonical form,
        # --> won't be needed anymore in Python 3.12
        IF SYS_VERSION_INFO_MAJOR <= 3 and SYS_VERSION_INFO_MINOR < 12:
            PyUnicode_READY(sequence)

        seq = Sequence.__new__(Sequence)
        seq._allocate(PyUnicode_GET_LENGTH(sequence))

        kind = PyUnicode_KIND(sequence)
        data = PyUnicode_DATA(sequence)

        with nogil:
            for i, j in enumerate(range(0, seq.slen * 2, 2)):
                letter = PyUnicode_READ(kind, data, i)
                if letter == u'A' or letter == u'a':
                    seq.digits[i] = A
                elif letter == u'T' or letter == u't':
                    seq.digits[i] = T
                elif letter == u'G' or letter == u'g':
                    seq.digits[i] = G
                    gc_count += 1
                elif letter == u'C' or letter == u'c':
                    seq.digits[i] = C
                    gc_count += 1
                else:
                    seq.digits[i] = N
            if seq.slen > 0:
                seq.gc = (<double> gc_count) / (<double> seq.slen)

        return seq

    cdef inline bint _is_a(self, int i, int strand = 1) nogil:
        cdef uint8_t x

        if strand == 1:
            return self.digits[i] == A
        else:
            return self.digits[self.slen - 1 - i] == T

    cdef inline bint _is_g(self, int i, int strand = 1) nogil:
        cdef uint8_t x

        if strand == 1:
            return self.digits[i] == G
        else:
            return self.digits[self.slen - 1 - i] == C

    cdef inline bint _is_gc(self, int i, int strand = 1) nogil:
        cdef uint8_t x

        if strand == 1:
            x = self.digits[i]
        else:
            x = self.digits[self.slen - 1 - i]

        # NB(@althonos): In the original Prodigal implementation, any unknown
        #                character gets encoded as a C, so it gets counted
        #                when computing the GC percent. We reproduce this
        #                behaviour here, but a better solution would be to
        #                count only known letters.
        return x == C or x == G or x == N

    cdef inline bint _is_start(self, int i, int tt, int strand = 1) nogil:
        cdef uint8_t x0
        cdef uint8_t x1
        cdef uint8_t x2

        if strand == 1:
            x0 = self.digits[i]
            x1 = self.digits[i+1]
            x2 = self.digits[i+2]
        else:
            x0 = _translation[self.digits[self.slen - 1 - i]]
            x1 = _translation[self.digits[self.slen - 2 - i]]
            x2 = _translation[self.digits[self.slen - 3 - i]]

        # ATG
        if x0 == A and x1 == T and x2 == G:
            return True
        # Codes that only use ATG
        if tt == 6 or tt == 10 or tt == 14 or tt == 15 or tt == 16 or tt == 2:
            return False
        # GTG
        if x0 == G and x1 == T and x2 == G:
            return not (tt == 1 or tt == 3 or tt == 12 or tt == 2)
        # TTG
        if x0 == T and x1 == T and x2 == G:
            return not (tt < 4 or tt == 9 or (tt >= 21 and tt < 25))
        # other codons
        return False

    cdef inline bint _is_stop(self, int i, int tt, int strand = 1) nogil:
        cdef uint8_t x0
        cdef uint8_t x1
        cdef uint8_t x2

        if strand == 1:
            x0 = self.digits[i]
            x1 = self.digits[i+1]
            x2 = self.digits[i+2]
        else:
            x0 = _translation[self.digits[self.slen - 1 - i]]
            x1 = _translation[self.digits[self.slen - 2 - i]]
            x2 = _translation[self.digits[self.slen - 3 - i]]

        # TAG
        if x0 == T and x1 == A and x2 == G:
            return not (tt == 6 or tt == 15 or tt == 16 or tt == 22)
        # TGA
        if x0 == T and x1 == G and x2 == A:
            return not (
                    tt == 2 or tt == 3 or tt == 4 or tt == 5
                 or tt == 9 or tt == 10 or tt == 13 or tt == 14
                 or tt == 21 or tt == 25
            )
        # TAA
        if x0 == T and x1 == A and x2 == A:
            return not (tt == 6 or tt == 14)
        # Code 2
        if tt == 2:
            return x0 == A and x1 == G and (x2 == A or x2 == G)
        # Code 22
        elif tt == 22:
            return x0 == T and x1 == C and x2 == A
        # Code 23
        elif tt == 23:
            return x0 == T and x1 == T and x2 == A
        # other codons
        return False

    cdef inline bint _is_atg(self, int i, int strand = 1) nogil:
        cdef uint8_t x0
        cdef uint8_t x1
        cdef uint8_t x2
        if strand == 1:
            x0 = self.digits[i]
            x1 = self.digits[i+1]
            x2 = self.digits[i+2]
        else:
            x0 = _translation[self.digits[self.slen - 1 - i]]
            x1 = _translation[self.digits[self.slen - 2 - i]]
            x2 = _translation[self.digits[self.slen - 3 - i]]
        return x0 == A and x1 == T and x2 == G

    cdef inline bint _is_gtg(self, int i, int strand = 1) nogil:
        cdef uint8_t x0
        cdef uint8_t x1
        cdef uint8_t x2
        if strand == 1:
            x0 = self.digits[i]
            x1 = self.digits[i+1]
            x2 = self.digits[i+2]
        else:
            x0 = _translation[self.digits[self.slen - 1 - i]]
            x1 = _translation[self.digits[self.slen - 2 - i]]
            x2 = _translation[self.digits[self.slen - 3 - i]]
        return x0 == G and x1 == T and x2 == G

    cdef inline bint _is_ttg(self, int i, int strand = 1) nogil:
        cdef uint8_t x0
        cdef uint8_t x1
        cdef uint8_t x2
        if strand == 1:
            x0 = self.digits[i]
            x1 = self.digits[i+1]
            x2 = self.digits[i+2]
        else:
            x0 = _translation[self.digits[self.slen - 1 - i]]
            x1 = _translation[self.digits[self.slen - 2 - i]]
            x2 = _translation[self.digits[self.slen - 3 - i]]
        return x0 == T and x1 == T and x2 == G

    cdef inline int _mer_ndx(self, int i, int length, int strand = 1) nogil:

        cdef int     j
        cdef int     k
        cdef uint8_t x
        cdef int     ndx = 0

        if strand == 1:
            for j, k in enumerate(range(i, i+length)):
                x = self.digits[k]
                ndx |= (x & 0b11) << 2*j
        else:
            for j, k in enumerate(range(self.slen - 1 - i, self.slen - 1 - i - length, -1)):
                x = self.digits[k]
                ndx |= (_translation[x] & 0b11) << 2*j

        return ndx

    cdef char _amino(self, int i, int tt, int strand = 1, bint is_init = False) nogil:

        cdef uint8_t x0
        cdef uint8_t x1
        cdef uint8_t x2

        if strand == 1:
            x0 = self.digits[i]
            x1 = self.digits[i+1]
            x2 = self.digits[i+2]
        else:
            x0 = _translation[self.digits[self.slen - 1 - i]]
            x1 = _translation[self.digits[self.slen - 2 - i]]
            x2 = _translation[self.digits[self.slen - 3 - i]]

        if self._is_stop(i, tt, strand=strand):
            return b"*"
        if self._is_start(i, tt, strand=strand) and is_init:
            return b"M"
        if x0 == T and x1 == T and (x2 == T or x2 == C):
            return b"F"
        if x0 == T and x1 == T and (x2 == A or x2 == G):
            return b"L"
        if x0 == T and x1 == C and x2 != N:
            return b"S"
        if x0 == T and x1 == A and x2 == T:
            return b"Y"
        if x0 == T and x1 == A and x2 == C:
            return b"Y"
        if x0 == T and x1 == A and x2 == A:
            if tt == 6:
                return b"Q"
            elif tt == 14:
                return b"Y"
        if x0 == T and x1 == A and x2 == G:
            if tt == 6 or tt == 15:
                return b"Q"
            elif tt == 22:
                return b"L"
        if x0 == T and x1 == G and (x2 == T or x2 == C):
            return b"C"
        if x0 == T and x1 == G and x2 == A:
            return b"G" if tt == 25 else b"W"
        if x0 == T and x1 == G and x2 == G:
            return b"W"
        if x0 == C and x1 == T and (x2 == T or x2 == C or x2 == A):
            return b"T" if tt == 3 else b"L"
        if x0 == C and x1 == T and x2 == G:
            return b"T" if tt == 3 else b"S" if tt == 12 else b"L"
        if x0 == C and x1 == C and x2 != N:
            return b"P"
        if x0 == C and x1 == A and (x2 == T or x2 == C):
            return b"H"
        if x0 == C and x1 == A and (x2 == A or x2 == G):
            return b"Q"
        if x0 == C and x1 == G and x2 != N:
            return b"R"
        if x0 == A and x1 == T and (x2 == T or x2 == C):
            return b"I"
        if x0 == A and x1 == T and x2 == A:
            return b"M" if tt == 2 or tt == 3 or tt == 5 or tt == 13 or tt == 22 else b"I"
        if x0 == A and x1 == T and x2 == G:
            return b"M"
        if x0 == A and x1 == C and x2 != N:
            return b"T"
        if x0 == A and x1 == A and (x2 == T or x2 == C):
            return b"N"
        if x0 == A and x1 == A and x2 == A:
            return b"N" if tt == 9 or tt == 14 or tt == 21 else b"K"
        if x0 == A and x1 == A and x2 == G:
            return b"K"
        if x0 == A and x1 == G and (x2 == T or x2 == C):
            return b"S"
        if x0 == A and x1 == G and (x2 == A or x2 == G):
            return b"G" if tt == 13 else b"S" if tt == 5 or tt == 9 or tt == 14 or tt == 21 else b"R"
        if x0 == G and x1 == T and x2 != N:
            return b"V"
        if x0 == G and x1 == C and x2 != N:
            return b"A"
        if x0 == G and x1 == A and (x2 == T or x2 == C):
            return b"D"
        if x0 == G and x1 == A and (x2 == A or x2 == G):
            return b"E"
        if x0 == G and x1 == G and x2 != N:
            return b"G"

        return b'X'

# --- Connection Scorer ------------------------------------------------------

_TARGET_CPU           = TARGET_CPU
_AVX2_RUNTIME_SUPPORT = False
_NEON_RUNTIME_SUPPORT = False
_SSE2_RUNTIME_SUPPORT = False
_AVX2_BUILD_SUPPORT   = False
_NEON_BUILD_SUPPORT   = False
_SSE2_BUILD_SUPPORT   = False

IF TARGET_CPU == "x86":
    cdef X86Info cpu_info = GetX86Info()
    _SSE2_RUNTIME_SUPPORT = cpu_info.features.sse2 != 0
    _AVX2_RUNTIME_SUPPORT = cpu_info.features.avx2 != 0
    _SSE2_BUILD_SUPPORT   = SSE2_BUILD_SUPPORT
    _AVX2_BUILD_SUPPORT   = AVX2_BUILD_SUPPORT
ELIF TARGET_CPU == "arm":
    cdef ArmInfo cpu_info = GetArmInfo()
    _NEON_RUNTIME_SUPPORT = cpu_info.features.neon != 0
    _NEON_BUILD_SUPPORT   = NEON_BUILD_SUPPORT
ELIF TARGET_CPU == "aarch64":
    _NEON_RUNTIME_SUPPORT = True
    _NEON_BUILD_SUPPORT   = NEON_BUILD_SUPPORT

cdef enum simd_backend:
    NONE = 0
    SSE2 = 1
    AVX2 = 2
    NEON = 3

cdef class ConnectionScorer:

    def __cinit__(self):
        self.capacity = 0
        self.skip_connection = self.skip_connection_raw = NULL
        self.node_types      = self.node_types_raw      = NULL
        self.node_strands    = self.node_strands_raw    = NULL
        self.node_frames     = self.node_frames_raw     = NULL

    def __init__(self, unicode backend="detect"):
        IF TARGET_CPU == "x86":
            if backend =="detect":
                self.backend = simd_backend.NONE
                IF SSE2_BUILD_SUPPORT:
                    if _SSE2_RUNTIME_SUPPORT:
                        self.backend = simd_backend.SSE2
                IF AVX2_BUILD_SUPPORT:
                    if _AVX2_RUNTIME_SUPPORT:
                        self.backend = simd_backend.AVX2
            elif backend == "sse":
                IF not SSE2_BUILD_SUPPORT:
                    raise RuntimeError("Extension was compiled without SSE2 support")
                ELSE:
                    if not _SSE2_RUNTIME_SUPPORT:
                        raise RuntimeError("Cannot run SSE2 instructions on this machine")
                    self.backend = simd_backend.SSE2
            elif backend == "avx":
                IF not AVX2_BUILD_SUPPORT:
                    raise RuntimeError("Extension was compiled without AVX2 support")
                ELSE:
                    if not _AVX2_RUNTIME_SUPPORT:
                        raise RuntimeError("Cannot run AVX2 instructions on this machine")
                    self.backend = simd_backend.AVX2
            elif backend is None:
                self.backend = simd_backend.NONE
            else:
                raise ValueError(f"Unsupported backend on this architecture: {backend}")
        ELIF TARGET_CPU == "arm" or TARGET_CPU == "aarch64":
            if backend =="detect":
                self.backend = simd_backend.NONE
                IF NEON_BUILD_SUPPORT:
                    if _NEON_RUNTIME_SUPPORT:
                        self.backend = simd_backend.NEON
            elif backend == "neon":
                IF not NEON_BUILD_SUPPORT:
                    raise RuntimeError("Extension was compiled without NEON support")
                ELSE:
                    if not _NEON_RUNTIME_SUPPORT:
                        raise RuntimeError("Cannot run NEON instructions on this machine")
                    self.backend = simd_backend.NEON
            elif backend is None:
                self.backend = simd_backend.NONE
            else:
                raise ValueError(f"Unsupported backend on this architecture: {backend}")
        ELSE:
            self.backend = simd_backend.NONE

    def __dealloc__(self):
        PyMem_Free(self.node_types_raw)
        PyMem_Free(self.node_strands_raw)
        PyMem_Free(self.node_frames_raw)
        PyMem_Free(self.skip_connection_raw)

    cdef int _index(self, Nodes nodes) nogil except 1:
        cdef size_t i
        # reallocate if needed
        if self.capacity < nodes.length:
            with gil:
                # reallocate new memory
                self.skip_connection_raw = <uint8_t*> PyMem_Realloc(self.skip_connection_raw, nodes.length * sizeof(uint8_t) + 0x1F)
                self.node_types_raw      = <uint8_t*> PyMem_Realloc(self.node_types_raw, nodes.length      * sizeof(uint8_t) + 0x1F)
                self.node_strands_raw    = <int8_t*>  PyMem_Realloc(self.node_strands_raw, nodes.length    * sizeof(int8_t)  + 0x1F)
                self.node_frames_raw     = <uint8_t*> PyMem_Realloc(self.node_frames_raw, nodes.length     * sizeof(uint8_t) + 0x1F)
                # check that allocations were successful
                if self.skip_connection_raw == NULL:
                    raise MemoryError("Failed to allocate memory for scoring bypass index")
                if self.node_types_raw == NULL:
                    raise MemoryError("Failed to allocate memory for node type array")
                if self.node_strands_raw == NULL:
                    raise MemoryError("Failed to allocate memory for node strand array")
                if self.node_frames_raw == NULL:
                    raise MemoryError("Failed to allocate memory for node frame array")
            # record new capacity
            self.capacity = nodes.length
            # compute pointers to aligned memory
            self.skip_connection = <uint8_t*> ((<uintptr_t> self.skip_connection_raw + 0x1F) & (~0x1F))
            self.node_types      = <uint8_t*> ((<uintptr_t> self.node_types_raw      + 0x1F) & (~0x1F))
            self.node_strands    = <int8_t*>  ((<uintptr_t> self.node_strands_raw    + 0x1F) & (~0x1F))
            self.node_frames     = <uint8_t*> ((<uintptr_t> self.node_frames_raw     + 0x1F) & (~0x1F))
        # copy data from the array of nodes
        for i in range(nodes.length):
            self.node_types[i]      = nodes.nodes[i].type
            self.node_strands[i]    = nodes.nodes[i].strand
            self.node_frames[i]     = nodes.nodes[i].ndx % 3
            self.skip_connection[i] = False
        # return 0 if no exceptions were raised
        return 0

    def index(self, Nodes nodes not None):
        with nogil:
            self._index(nodes)

    cdef int _compute_skippable(self, int min, int i) nogil except 1:
        if self.backend != simd_backend.NONE:
            memset(&self.skip_connection[min], 0, sizeof(uint8_t) * (i - min))
        IF AVX2_BUILD_SUPPORT:
            if self.backend == simd_backend.AVX2:
                skippable_avx(self.node_strands, self.node_types, self.node_frames, min, i, self.skip_connection)
        IF SSE2_BUILD_SUPPORT:
            if self.backend == simd_backend.SSE2:
                skippable_sse(self.node_strands, self.node_types, self.node_frames, min, i, self.skip_connection)
        IF NEON_BUILD_SUPPORT:
            if self.backend == simd_backend.NEON:
                skippable_neon(self.node_strands, self.node_types, self.node_frames, min, i, self.skip_connection)
        return 0

    def compute_skippable(self, int min, int i):
        assert self.skip_connection != NULL
        assert i < <int> self.capacity
        assert min <= i
        with nogil:
            self._compute_skippable(min, i)

    cdef int _score_connections(self, Nodes nodes, int min, int i, TrainingInfo tinf, bint final=False) nogil except 1:
        cdef int j
        cdef _node*     raw_nodes = nodes.nodes
        cdef _training* raw_tinf  = tinf.tinf

        for j in range(min, i):
            if self.skip_connection[j] == 0:
                dprog.score_connection(raw_nodes, j, i, raw_tinf, final)

        return 0

    def score_connections(self, Nodes nodes not None, int min, int i, TrainingInfo tinf not None, bint final=False):
        assert self.skip_connection != NULL
        assert i < <int> nodes.length
        assert min <= i
        with nogil:
            self._score_connections(nodes, min, i, tinf, final)

# --- Nodes ------------------------------------------------------------------

cdef class Node:
    """A dynamic programming node used by Prodigal to score ORFs.

    .. versionadded:: 0.5.4

    """

    @property
    def type(self):
        """`str`: The node type (ATG, GTG, TTG, or Stop).
        """
        assert self.node != NULL
        return ["ATG", "GTG", "TTG" , "Stop"][self.node.type]

    @property
    def edge(self):
        """`bool`: `True` if the node runs off the edge.
        """
        assert self.node != NULL
        return self.node.edge

    @property
    def gc_bias(self):
        """`int`: The frame of highest GC content within this node.
        """
        assert self.node != NULL
        return self.node.gc_bias

    @property
    def cscore(self):
        """`float`: The coding score for this node, based on 6-mer usage.
        """
        assert self.node != NULL
        return self.node.cscore

    @property
    def gc_cont(self):
        """`float`: The GC content for the node.
        """
        assert self.node != NULL
        return self.node.gc_cont

    @property
    def score(self):
        """`float`: The score of the total solution up to this node.
        """
        assert self.node != NULL
        return self.node.score

cdef class Nodes:
    """A list of dynamic programming nodes used by Prodigal to score ORFs.

    .. versionadded:: 0.5.4

    """

    def __cinit__(self):
        self.nodes = NULL
        self.capacity = 0
        self.length = 0

    def __init__(self):
        self._clear()

    def __dealloc__(self):
        PyMem_Free(self.nodes)

    def __copy__(self):
        return self.copy()

    def __len__(self):
        return self.length

    def __getitem__(self, ssize_t index):
        cdef Node node
        if index < 0:
            index += <ssize_t> self.length
        if index >= <ssize_t> self.length or index < 0:
            raise IndexError("list index out of range")
        node = Node.__new__(Node)
        node.owner = self
        node.node = &self.nodes[index]
        return node

    def __sizeof__(self):
        return self.capacity * sizeof(_node) + sizeof(self)

    cdef inline _node* _add_node(
        self,
        const int  ndx,
        const int  type,
        const int  strand,
        const int  stop_val,
        const bint edge,
    ) nogil except NULL:
        """Add a single node to the vector, and return a pointer to that node.
        """

        cdef size_t old_capacity = self.capacity
        cdef _node* node

        if self.length >= self.capacity:
            self.capacity = MIN_NODES_ALLOC if self.capacity == 0 else self.capacity*2
            with gil:
                self.nodes = <_node*> PyMem_Realloc(self.nodes, self.capacity * sizeof(_node))
                if self.nodes == NULL:
                    raise MemoryError("Failed to reallocate node array")
            memset(&self.nodes[old_capacity], 0, (self.capacity - old_capacity) * sizeof(_node))

        self.length += 1
        node = &self.nodes[self.length - 1]
        node.ndx = ndx
        node.type = type
        node.strand = strand
        node.stop_val = stop_val
        node.edge = edge
        return node

    cdef int _clear(self) nogil except 1:
        """Remove all nodes from the vector.
        """
        cdef size_t old_length
        old_length, self.length = self.length, 0
        memset(self.nodes, 0, old_length * sizeof(_node))

    def clear(self):
        """Remove all nodes from the vector.
        """
        with nogil:
            self._clear()

    cpdef Nodes copy(self):
        cdef Nodes new = Nodes.__new__(Nodes)
        new.capacity = self.capacity
        new.length = self.length
        new.nodes = <_node*> PyMem_Malloc(new.capacity * sizeof(_node))
        if new.nodes == NULL:
            raise MemoryError("Failed to reallocate node array")
        memcpy(new.nodes, self.nodes, new.capacity * sizeof(_node))
        return new

    cdef int _sort(self) nogil except 1:
        """Sort all nodes in the vector by their index and strand.
        """
        qsort(self.nodes, self.length, sizeof(_node), node.compare_nodes)

    def sort(self):
        """Sort all nodes in the vector by their index and strand.
        """
        with nogil:
            self._sort()

# --- Genes ------------------------------------------------------------------

cdef class Gene:
    """A single raw gene found by Prodigal within a DNA sequence.

    .. versionadded:: 0.5.4

    """

    @property
    def begin(self):
        """`int`: The coordinate at which the gene begins.
        """
        return self.gene.begin

    @property
    def end(self):
        """`int`: The coordinate at which the gene ends.
        """
        return self.gene.end

    @property
    def start_ndx(self):
        """`int`: The index of the start node in the `Nodes` list.
        """
        return self.gene.start_ndx

    @property
    def stop_ndx(self):
        """`int`: The index of the stop node in the `Nodes` list.
        """
        return self.gene.stop_ndx

cdef class Genes:
    """A list of raw genes found by Prodigal in a single sequence.

    .. versionadded:: 0.5.4

    """

    def __cinit__(self):
        self.genes = NULL
        self.capacity = 0
        self.length = 0

    def __init__(self):
        self._clear()

    def __dealloc__(self):
        PyMem_Free(self.genes)

    def __len__(self):
        return self.length

    def __getitem__(self, ssize_t index):
        cdef Gene gene
        if index < 0:
            index += <ssize_t> self.length
        if index >= <ssize_t> self.length or index < 0:
            raise IndexError("list index out of range")
        gene = Gene.__new__(Gene)
        gene.owner = self
        gene.gene = &self.genes[index]
        return gene

    def __sizeof__(self):
        return self.capacity * sizeof(_gene) + sizeof(self)

    cdef inline _gene* _add_gene(
        self,
        const int begin,
        const int end,
        const int start_ndx,
        const int stop_ndx,
    ) nogil except NULL:
        """Add a single gene to the vector, and return a pointer to that gene.
        """

        cdef size_t old_capacity = self.capacity
        cdef _gene* gene

        if self.length >= self.capacity:
            self.capacity = MIN_GENES_ALLOC if self.capacity == 0 else self.capacity*2
            with gil:
                self.genes = <_gene*> PyMem_Realloc(self.genes, self.capacity * sizeof(_gene))
                if self.genes == NULL:
                    raise MemoryError("Failed to reallocate gene array")
            memset(&self.genes[old_capacity], 0, (self.capacity - old_capacity) * sizeof(_gene))

        self.length += 1
        gene = &self.genes[self.length - 1]
        gene.begin = begin
        gene.end = end
        gene.start_ndx = start_ndx
        gene.stop_ndx = stop_ndx
        return gene

    cdef int _clear(self) nogil except 1:
        """Remove all genes from the vector.
        """
        cdef size_t old_length
        old_length, self.length = self.length, 0
        memset(self.genes, 0, old_length * sizeof(_gene))

    def clear(self):
        """Remove all genes from the vector.
        """
        with nogil:
            self._clear()

# --- Training Info ----------------------------------------------------------

cdef class TrainingInfo:
    """A collection of parameters obtained after training.

    .. versionadded:: 0.5.0

    """

    def __cinit__(self):
        self.owned = False
        self.tinf = NULL

    def __init__(self, double gc, double start_weight=4.35, int translation_table=11):
        if self.tinf != NULL:
            raise RuntimeError("TrainingInfo.__init__ called more than once")
        # reallocate the training info memory, if needed
        self.owned = True
        self.tinf = <_training*> PyMem_Malloc(sizeof(_training))
        if self.tinf == NULL:
            raise MemoryError("Failed to allocate training info")
        # clear memory
        memset(self.tinf, 0, sizeof(_training))
        # set the variables
        self.tinf.gc = gc
        self.tinf.st_wt = start_weight
        self.tinf.trans_table = translation_table

    def __dealloc__(self):
        if self.owned:
            PyMem_Free(self.tinf)

    @property
    def translation_table(self):
        """`int`: The translation table used during training.
        """
        assert self.tinf != NULL
        return self.tinf.trans_table

    @translation_table.setter
    def translation_table(self, int table):
        assert self.tinf != NULL
        if table not in TRANSLATION_TABLES:
            raise ValueError(f"{table} is not a valid translation table index")
        self.raw.trans_table = table

    @property
    def gc(self):
        """`float`: The GC content of the training sequence.
        """
        assert self.tinf != NULL
        return self.tinf.gc

    @gc.setter
    def gc(self, double gc):
        assert self.tinf != NULL
        self.tinf.gc = gc

    @property
    def bias(self):
        """`tuple` of `float`: The GC frame bias for each of the 3 positions.
        """
        assert self.tinf != NULL
        return tuple(self.tinf.bias)

    @bias.setter
    def bias(self, object bias):
        assert self.tinf != NULL
        self.tinf.bias = bias

    @property
    def type_weights(self):
        """`tuple` of `float`: The weights for the ATG, GTG and TTG codons.
        """
        assert self.tinf != NULL
        return tuple(self.tinf.type_wt)

    @type_weights.setter
    def type_weights(self, object type_weights):
        assert self.tinf != NULL
        self.tinf.type_wt = type_weights

    @property
    def uses_sd(self):
        """`bool`: `True` if the sequence uses a Shine/Dalgarno motif.
        """
        return self.tinf.uses_sd

    @uses_sd.setter
    def uses_sd(self, bint uses_sd):
        assert self.tinf != NULL
        self.tinf.uses_sd = uses_sd

    @property
    def start_weight(self):
        """`float`: The start score weight to use for the training sequence.
        """
        assert self.tinf != NULL
        return self.tinf.st_wt

    @start_weight.setter
    def start_weight(self, double st_wt):
        self.tinf.st_wt = st_wt

# --- Metagenomic Bins -------------------------------------------------------

cdef class MetagenomicBin:
    """A pre-trained collection used to find genes in metagenomic mode.
    """

    @property
    def index(self):
        assert self.bin != NULL
        return self.bin.index

    @property
    def description(self):
        assert self.bin != NULL
        return self.bin.desc.decode('ascii')

# Allocate raw C memory for the C structs
cdef _metagenomic_bin _METAGENOMIC_BINS[NUM_META]
cdef ssize_t _i
for _i in range(NUM_META):
    memset(&_METAGENOMIC_BINS[_i], 0, sizeof(_metagenomic_bin))
    _METAGENOMIC_BINS[_i].tinf = <_training*> PyMem_Malloc(sizeof(_training))
    if not _METAGENOMIC_BINS[_i].tinf:
        raise MemoryError()
    memset(_METAGENOMIC_BINS[_i].tinf, 0, sizeof(_training))
initialize_metagenomic_bins(_METAGENOMIC_BINS)

# Create a tuple of objects exposing the C metagenomic bins
cdef MetagenomicBin _bin
cdef tuple _m = PyTuple_New(NUM_META)
for _i in range(NUM_META):
    _bin = MetagenomicBin.__new__(MetagenomicBin, )
    _bin.bin = &_METAGENOMIC_BINS[_i]
    _bin.training_info = TrainingInfo.__new__(TrainingInfo)
    _bin.training_info.owned = False
    _bin.training_info.tinf = _bin.bin.tinf
    PyTuple_SET_ITEM(_m, _i, _bin)
    Py_INCREF(_bin)
METAGENOMIC_BINS = _m

# --- Predictions ------------------------------------------------------------

cdef class Prediction:
    """A single predicted gene found by Prodigal.
    """

    @property
    def _gene_data(self):
        return self.gene.gene.gene_data.decode('ascii')

    @property
    def _score_data(self):
        return self.gene.gene.score_data.decode('ascii')

    @property
    def begin(self):
        """`int`: The coordinate at which the gene begins.
        """
        return self.gene.begin

    @property
    def end(self):
        """`int`: The coordinate at which the gene ends.
        """
        return self.gene.end

    @property
    def strand(self):
        """`int`: *-1* if the gene is on the reverse strand, *+1* otherwise.
        """
        return self.owner.nodes.nodes[self.gene.gene.start_ndx].strand

    @property
    def partial_begin(self):
        """`bool`: whether the gene overlaps with the start of the sequence.
        """
        if self.strand == 1:
            return self.owner.nodes.nodes[self.gene.gene.start_ndx].edge == 1
        else:
            return self.owner.nodes.nodes[self.gene.gene.stop_ndx].edge == 1

    @property
    def partial_end(self):
        """`bool`: whether the gene overlaps with the end of the sequence.
        """
        if self.strand == 1:
            return self.owner.nodes.nodes[self.gene.gene.stop_ndx].edge == 1
        else:
            return self.owner.nodes.nodes[self.gene.gene.start_ndx].edge == 1

    @property
    def start_type(self):
        """`str`: The start codon of this gene.

        Can be one of ``ATG``, ``GTG`` or ``TTG``, or ``Edge`` if `Pyrodigal`
        has been initialized in open ends mode and the gene starts right at the
        beginning of the input sequence.
        """
        node = self.owner.nodes.nodes[self.gene.gene.start_ndx]
        start_type = 3 if node.edge else node.type
        return ["ATG", "GTG", "TTG" , "Edge"][start_type]

    @property
    def rbs_motif(self):
        """``str``, optional: The motif of the Ribosome Binding Site.

        Possible non-`None` values are ``GGA/GAG/AGG``, ``3Base/5BMM``,
        ``4Base/6BMM``, ``AGxAG``, ``GGxGG``, ``AGGAG(G)/GGAGG``, ``AGGA``,
        ``AGGA/GGAG/GAGG``, ``GGAG/GAGG``, ``AGGAG/GGAGG``, ``AGGAG``,
        ``GGAGG`` or ``AGGAGG``.

        """
        cdef char* data = self.gene.gene.gene_data
        cdef char* i = strstr(data, "rbs_motif")
        cdef char* j = <char*> memchr(i, b';', 30)
        cdef size_t length = j - i
        if i[10:length] == b"None":
            return None
        return i[10:length].decode("ascii")

    @property
    def rbs_spacer(self):
        """`str`, optional: The number of base pair between the RBS and the CDS.

        Possible non-`None` values are ``3-4bp``, ``5-10bp``, ``11-12bp`` or
        ``13-15bp``.

        """
        cdef char* data = self.gene.gene.gene_data
        cdef char* i = strstr(data, "rbs_spacer")
        cdef char* j = <char*> memchr(i, b';', 30)
        cdef size_t length = j - i
        if i[11:length] == b"None":
            return None
        return i[11:length].decode("ascii")

    @property
    def gc_cont(self):
        """`float`: The GC content of the gene (between *0* and *1*).
        """
        return self.owner.nodes.nodes[self.gene.gene.start_ndx].gc_cont

    @property
    def translation_table(self):
        """`int`: The translation table used to find the gene.
        """
        return self.owner.training_info.translation_table

    @property
    def cscore(self):
        """`float`: The coding score for the start node, based on 6-mer usage.

        .. versionadded:: 0.5.1

        """
        return self.owner.nodes.nodes[self.gene.gene.start_ndx].cscore

    @property
    def rscore(self):
        """`float`: The score for the RBS motif.

        .. versionadded:: 0.5.1

        """
        return self.owner.nodes.nodes[self.gene.gene.start_ndx].rscore

    @property
    def sscore(self):
        """`float`: The score for the strength of the start codon.

        .. versionadded:: 0.5.1

        """
        return self.owner.nodes.nodes[self.gene.gene.start_ndx].sscore

    @property
    def tscore(self):
        """`float`: The score for the codon kind (ATG/GTG/TTG).

        .. versionadded:: 0.5.1

        """
        return self.owner.nodes.nodes[self.gene.gene.start_ndx].tscore

    @property
    def uscore(self):
        """`float`: The score for the upstream regions.

        .. versionadded:: 0.5.1

        """
        return self.owner.nodes.nodes[self.gene.gene.start_ndx].uscore

    cpdef double confidence(self):
        """confidence(self)\n--

        Estimate the confidence of the prediction.

        Returns:
            `float`: A confidence percentage between *0* and *100*.

        """
        cdef int ndx = self.gene.gene.start_ndx
        return gene.calculate_confidence(
            self.nodes.nodes[ndx].cscore + self.nodes.nodes[ndx].sscore,
            self.owner.training_info.tinf.st_wt
        )

    cpdef unicode sequence(self):
        """sequence(self)\n--

        Build the nucleotide sequence of this predicted gene.

        .. versionadded:: 0.5.4

        """

        cdef size_t   i
        cdef size_t   j
        cdef size_t   begin
        cdef size_t   end
        cdef size_t   unk
        cdef size_t   length
        cdef Py_UCS4  nuc
        cdef _gene*   gene   = self.gene.gene
        cdef int      slen   = self.owner.sequence.slen
        cdef int      strand = self.owner.nodes.nodes[gene.start_ndx].strand
        cdef uint8_t* digits   = self.owner.sequence.digits

        # compute the right length to hold the nucleotides
        length = (<size_t> gene.end) - (<size_t> gene.begin) + 1

        # compute the offsets in the sequence bitmap
        if strand == 1:
            begin = gene.begin - 1
            end = gene.end
            unk = gene.begin - 1
        else:
            begin = slen - gene.end
            end = slen + 1 - gene.begin
            unk = gene.end - 1

        # NB(@althonos): For some reason, PyPy3.6 (v7.3.3) is not happy with
        #                the use of the PyUnicode API here, and will just not
        #                write any letter with PyUnicode_WRITE. The bug
        #                doesn't seem to affect `Prediction.translate`, so
        #                I'm not sure what's going on, but in that case we
        #                can build an ASCII string and decode afterwards.
        IF SYS_VERSION_INFO_MAJOR <= 3 and SYS_VERSION_INFO_MINOR < 7 and SYS_IMPLEMENTATION_NAME == "pypy":
            cdef bytes dna
            cdef int   kind
            cdef void* data
            # create an empty byte buffer that we can write to
            dna = PyBytes_FromStringAndSize(NULL, length)
            data = <void*> PyBytes_AsString(dna)
        ELSE:
            cdef unicode  dna
            cdef int      kind
            cdef void*    data
            # create an empty string that we can write to
            dna  = PyUnicode_New(length, 0x7F)
            kind = PyUnicode_KIND(dna)
            data = PyUnicode_DATA(dna)

        with nogil:
            for i, j in enumerate(range(begin, end)):
                if strand == 1:
                    nuc = _letters[digits[j]]
                else:
                    nuc = _letters[_translation[digits[slen - 1 - j]]]
                IF SYS_VERSION_INFO_MAJOR <= 3 and SYS_VERSION_INFO_MINOR < 7 and SYS_IMPLEMENTATION_NAME == "pypy":
                    (<char*> data)[i] = nuc
                ELSE:
                    PyUnicode_WRITE(kind, data, i, nuc)
                unk += strand

        IF SYS_VERSION_INFO_MAJOR <= 3 and SYS_VERSION_INFO_MINOR < 7 and SYS_IMPLEMENTATION_NAME == "pypy":
            return dna.decode("ascii")
        ELSE:
            return dna

    cpdef unicode translate(
        self,
        object translation_table=None,
        Py_UCS4 unknown_residue=u"X",
    ):
        """translate(self, translation_table=None, unknown_residue="X")\n--

        Translate the predicted gene into a protein sequence.

        Arguments:
            translation_table (`int`, optional): An alternative translation table
                to use to translate the gene. Use ``None`` (the default) to
                translate using the translation table this gene was found with.
            unknown_residue (`str`): A single character to use for residues
                translated from codons with unknown nucleotides.

        Returns:
            `str`: The proteins sequence as a string using the right translation
            table and the standard single letter alphabet for proteins.

        Raises:
            `ValueError`: when ``translation_table`` is not a valid
                genetic code number.

        """

        cdef size_t         nucl_length
        cdef size_t         prot_length
        cdef size_t         i
        cdef size_t         j
        cdef int            tt
        cdef object         protein
        cdef int            kind
        cdef void*          data
        cdef Py_UCS4        aa
        cdef _gene*         gene        = self.gene.gene
        cdef int            slen        = self.owner.sequence.slen
        cdef int            edge        = self.owner.nodes.nodes[gene.start_ndx].edge
        cdef int            strand      = self.owner.nodes.nodes[gene.start_ndx].strand
        cdef bitmap_t       seq
        cdef size_t         begin
        cdef size_t         end
        cdef size_t         unk

        # HACK: support changing the translation table (without allocating a
        #       new a training info structure) by manipulating where the
        #       table would be read from in the fields of the struct
        if translation_table is None:
            tt = self.owner.training_info.tinf.trans_table
        elif translation_table not in TRANSLATION_TABLES:
            raise ValueError(f"{translation_table} is not a valid translation table index")
        else:
            tt = translation_table

        # compute the right length to hold the protein
        nucl_length = (<size_t> gene.end) - (<size_t> gene.begin) + 1
        prot_length = nucl_length//3
        # create an empty protein string that we can write to
        # with the appropriate functions
        protein = PyUnicode_New(prot_length, 0x7F)
        kind    = PyUnicode_KIND(protein)
        data    = PyUnicode_DATA(protein)

        # compute the offsets in the sequence bitmaps:
        # - begin is the coordinates of the first nucleotide in the gene
        # - unk is the coordinate of the first nucleotide in the useq bitmap
        if strand == 1:
            begin = gene.begin - 1
            end = gene.end - 1
        else:
            begin = slen - gene.end
            end = slen - gene.begin

        with nogil:
            for i, j in enumerate(range(begin, end, 3)):
                aa = self.owner.sequence._amino(j, tt, strand=strand, is_init=i==0 and not edge)
                PyUnicode_WRITE(kind, data, i, aa)

        # return the string containing the protein sequence
        return protein

cdef class Predictions:
    """A list of predictions found by Prodigal in a single sequence.

    Attributes:
        sequence (`pyrodigal.Sequence`): The compressed input sequence for
            which the predictions were obtained.
        training_info (`pyrodigal.TrainingInfo`): A reference to the training
            info these predictions were obtained with.
        genes (`pyrodigal.Genes`): A list to the predicted genes in the
            input sequence.
        nodes (`pyrodigal.Nodes`): A list to the raw nodes found in the
            input sequence.

    """

    def __init__(self, Genes genes, Nodes nodes, Sequence sequence, TrainingInfo training_info):
        self.genes = genes
        self.nodes = nodes
        self.sequence = sequence
        self.training_info = training_info

    def __bool__(self):
        return self.genes.length > 0

    def __len__(self):
        return self.genes.length

    def __getitem__(self, ssize_t index):
        cdef Prediction pred
        if index < 0:
            index += <ssize_t> self.genes.length
        if index >= <ssize_t> self.genes.length or index < 0:
            raise IndexError("list index out of range")
        pred = Prediction.__new__(Prediction)
        pred.owner = self
        pred.gene = self.genes[index]
        return pred

# --- Pyrodigal --------------------------------------------------------------

cdef class Pyrodigal:
    """An efficient ORF finder for genomes, progenomes and metagenomes.

    Attributes:
        meta (`bool`): Whether or not this object is configured to
           find genes using the metagenomic bins or manually created
           training infos.
        closed (`bool`): Whether or not proteins can run off edges when
           finding genes in a sequence.
        training_info (`~pyrodigal.TrainingInfo`): The object storing the
           training information, or `None` if the object is in metagenomic
           mode or hasn't been trained yet.

    """

    def __cinit__(self):
        self._num_seq = 1

    def __init__(self, meta=False, closed=False):
        """__init__(self, meta=False, closed=False)\n--

        Instantiate and configure a new ORF finder.

        Arguments:
            meta (`bool`): Set to `True` to run in metagenomic mode, using a
                pre-trained profiles for better results with metagenomic or
                progenomic inputs. Defaults to `False`.
            closed (`bool`): Set to `True` to consider sequences ends 'closed',
                which prevents proteins from running off edges. Defaults to
                `False`.

        """
        self.meta = meta
        self.closed = closed
        self.lock = threading.Lock()

    cpdef Predictions find_genes(self, object sequence):
        """find_genes(self, sequence)\n--

        Find all the genes in the input DNA sequence.

        Arguments:
            sequence (`str` or buffer): The nucleotide sequence to use,
                either as a string of nucleotides, or as an object implementing
                the buffer protocol. Letters not corresponding to an usual
                nucleotide (not any of "ATGC") will be ignored.

        Returns:
            `Predictions`: A collection of all the genes found in the input.

        Raises:
            `MemoryError`: When allocation of an internal buffers fails.
            `RuntimeError`: On calling this method without `train` in *single* mode.
            `TypeError`: When ``sequence`` does not implement the buffer protocol.

        """
        cdef int         n
        cdef Sequence    seq
        cdef Predictions preds

        if not self.meta and self.training_info is None:
            raise RuntimeError("cannot find genes without having trained in single mode")

        if isinstance(sequence, str):
            seq = Sequence.from_string(sequence)
        else:
            seq = Sequence.from_bytes(sequence)

        with self.lock:
            n = self._num_seq
            self._num_seq += 1
        if self.meta:
            return find_genes_meta(seq, self.closed, n)
        else:
            return find_genes_single(seq, self.training_info, self.closed, n)

    cpdef TrainingInfo train(self, object sequence, bint force_nonsd=False, double st_wt=4.35, int translation_table=11):
        """train(self, sequence, force_nonsd=False, st_wt=4.35, translation_table=11)\n--

        Search optimal parameters for the ORF finder using a training sequence.

        Arguments:
            sequence (`str` or buffer): The nucleotide sequence to use, either
                as a string of nucleotides, or as an object implementing the
                buffer protocol.
            force_nonsd (`bool`, optional): Set to ``True`` to bypass the heuristic
                algorithm that tries to determine if the organism the training
                sequence belongs to uses a Shine-Dalgarno motif or not.
            st_wt (`float`, optional): The start score weight to use. The default
                value has been manually selected by the Prodigal authors as an
                appropriate value for 99% of genomes.
            translation_table (`int`, optional): The translation table to use.
                Check the `list of genetic codes <https://w.wiki/47wo>`_ for
                the available values.

        Raises:
            `MemoryError`: When allocation of an internal buffers fails.
            `RuntimeError`: When calling this method while in *metagenomic* mode.
            `TypeError`: When ``sequence`` does not implement the buffer protocol.
            `ValueError`: When ``translation_table`` is not a valid number.

        """
        cdef Sequence     seq
        cdef int          slen
        cdef TrainingInfo tinf

        if self.meta:
            raise RuntimeError("cannot use training sequence in metagenomic mode")
        if translation_table not in TRANSLATION_TABLES:
            raise ValueError(f"{translation_table} is not a valid translation table index")

        if isinstance(sequence, str):
            seq = Sequence.from_string(sequence)
        else:
            seq = Sequence.from_bytes(sequence)

        if seq.slen < MIN_SINGLE_GENOME:
            raise ValueError(
                f"sequence must be at least {MIN_SINGLE_GENOME} characters ({seq.slen} found)"
            )
        elif seq.slen < IDEAL_SINGLE_GENOME:
            warnings.warn(
                f"sequence should be at least {IDEAL_SINGLE_GENOME} characters ({seq.slen} found)"
            )

        tinf = train(seq, self.closed, force_nonsd, st_wt, translation_table)
        with self.lock:
            self.training_info = tinf

        return tinf

# --- C-level API reimplementation -------------------------------------------

cpdef int add_nodes(Nodes nodes, Sequence seq, TrainingInfo tinf, bint closed=False) nogil except -1:
    """Adds nodes to the node list, based on the sequence.

    Genes must be larger than 90bp in length, unless they run off the edge,
    in which case they only have to be 50bp.

    """
    cdef int    i
    cdef int    last[3]
    cdef int    min_dist[3]
    cdef bint   saw_start[3]
    cdef int    slmod        = seq.slen % 3
    cdef int    nn           = 0
    cdef int    tt           = tinf.tinf.trans_table
    # TODO(@althonos): handle region masking
    cdef _mask* mlist = NULL
    cdef int    nm    = 0

    # If sequence is smaller than a codon, there are no nodes to add
    if seq.slen < 3:
        return nn

    # Forward strand nodes
    for i in range(3):
        last[(i+slmod)%3] = seq.slen + i
        saw_start[i%3] = False
        min_dist[i%3] = MIN_EDGE_GENE
        if not closed:
            while last[(i+slmod)%3] + 3 > seq.slen:
                last[(i+slmod)%3] -= 3
    for i in reversed(range(seq.slen-2)):
        if seq._is_stop(i, tt):
            if saw_start[i%3]:
                nodes._add_node(
                    ndx = last[i%3],
                    type = node_type.STOP,
                    strand = 1,
                    stop_val = i,
                    edge = not seq._is_stop(last[i%3], tt),
                )
                nn += 1
            min_dist[i%3] = MIN_GENE
            last[i%3] = i
            saw_start[i%3] = False
            continue
        if last[i%3] >= seq.slen:
            continue
        if not cross_mask(i, last[i%3], mlist, nm):
            if last[i%3] - i + 3 >= min_dist[i%3] and seq._is_start(i, tt):
                if seq._is_atg(i):
                    saw_start[i%3] = True
                    nodes._add_node(
                        ndx = i,
                        type = node_type.ATG,
                        stop_val = last[i%3],
                        strand = 1,
                        edge = False
                    )
                    nn += 1
                elif seq._is_ttg(i):
                    saw_start[i%3] = True
                    nodes._add_node(
                        ndx = i,
                        type = node_type.TTG,
                        stop_val = last[i%3],
                        strand = 1,
                        edge = False
                    )
                    nn += 1
                elif seq._is_gtg(i):
                    saw_start[i%3] = True
                    nodes._add_node(
                        ndx = i,
                        type = node_type.GTG,
                        stop_val = last[i%3],
                        strand = 1,
                        edge = False
                    )
                    nn += 1
            if i <= 2 and not closed and last[i%3] - i > MIN_EDGE_GENE:
                saw_start[i%3] = True
                nodes._add_node(
                    ndx = i,
                    type = node_type.ATG,
                    stop_val = last[i%3],
                    strand = 1,
                    edge = True,
                )
                nn += 1
    for i in range(3):
        if saw_start[i%3]:
            nodes._add_node(
                ndx = last[i%3],
                type = node_type.STOP,
                strand = 1,
                stop_val = i - 6,
                edge = not seq._is_stop(last[i%3], tt)
            )
            nn += 1
    # Reverse strand nodes
    for i in range(3):
        last[(i + slmod) % 3] = seq.slen + i
        saw_start[i%3] = False
        min_dist[i%3] = MIN_EDGE_GENE
        if not closed:
            while last[(i+slmod) % 3] + 3 > seq.slen:
                last[(i+slmod)%3] -= 3
    for i in reversed(range(seq.slen-2)):
        if seq._is_stop(i, tt, strand=-1):
            if saw_start[i%3]:
                nodes._add_node(
                    ndx = seq.slen - last[i%3] - 1,
                    type = node_type.STOP,
                    strand = -1,
                    stop_val = seq.slen - i - 1,
                    edge = not seq._is_stop(last[i%3], tt, strand=-1)
                )
                nn += 1
            min_dist[i%3] = MIN_GENE
            last[i%3] = i
            saw_start[i%3] = False
            continue
        if last[i%3] >= seq.slen:
            continue
        if not cross_mask(i, last[i%3], mlist, nm):
            if last[i%3] - i + 3 >= min_dist[i%3] and seq._is_start(i, tt, strand=-1):
                if seq._is_atg(i, strand=-1):
                    saw_start[i%3] = True
                    nodes._add_node(
                        ndx = seq.slen - i - 1,
                        type = node_type.ATG,
                        strand = -1,
                        stop_val = seq.slen - last[i%3] - 1,
                        edge = False
                    )
                    nn += 1
                elif seq._is_gtg(i, strand=-1):
                    saw_start[i%3] = True
                    nodes._add_node(
                        ndx = seq.slen - i - 1,
                        type = node_type.GTG,
                        strand = -1,
                        stop_val = seq.slen - last[i%3] - 1,
                        edge = False
                    )
                    nn += 1
                elif seq._is_ttg(i, strand=-1):
                    saw_start[i%3] = 1
                    nodes._add_node(
                        ndx = seq.slen - i - 1,
                        type = node_type.TTG,
                        strand = -1,
                        stop_val = seq.slen - last[i%3] - 1,
                        edge = False,
                    )
                    nn += 1
        if i <= 2 and not closed and last[i%3] - i > MIN_EDGE_GENE:
            saw_start[i%3] = 1
            node = nodes._add_node(
                ndx = seq.slen - i - 1,
                type = node_type.ATG,
                strand = -1,
                stop_val = seq.slen - last[i%3] - 1,
                edge = True,
            )
            nn += 1
    for i in range(3):
        if saw_start[i%3]:
            node = nodes._add_node(
                ndx = seq.slen - last[i%3] - 1,
                type = node_type.STOP,
                strand = -1,
                stop_val = seq.slen - i + 5,
                edge = not seq._is_stop(last[i%3], tt, strand=-1),
            )
            nn += 1

    return nn

cpdef int add_genes(Genes genes, Nodes nodes, int ipath) nogil except -1:
    """Adds genes to the gene list, based on the nodes.
    """
    cdef int  path      = ipath
    cdef int  ng        = 0
    cdef int  begin     = 0
    cdef int  end       = 0
    cdef int  start_ndx = 0
    cdef int  stop_ndx  = 0

    if path == -1:
        return 0
    while nodes.nodes[path].traceb != -1:
        path = nodes.nodes[path].traceb
    while path != -1:
        if nodes.nodes[path].elim == 1:
            pass
        elif nodes.nodes[path].strand == 1:
            if nodes.nodes[path].type != node_type.STOP:
                begin = nodes.nodes[path].ndx + 1
                start_ndx = path
            else:
                end = nodes.nodes[path].ndx + 3
                stop_ndx = path
                genes._add_gene(begin, end, start_ndx, stop_ndx)
                ng += 1
        else:
            if nodes.nodes[path].type != node_type.STOP:
                end = nodes.nodes[path].ndx + 1
                start_ndx = path
                genes._add_gene(begin, end, start_ndx, stop_ndx)
                ng += 1
            else:
                begin = nodes.nodes[path].ndx - 1
                stop_ndx = path
        path = nodes.nodes[path].tracef

    return ng

cpdef void calc_dicodon_gene(TrainingInfo tinf, Sequence seq, Nodes nodes, int ipath) nogil:
    """Compute the dicodon frequency in genes and in the background.

    Stores the log-likelihood of each 6-mer relative to the background.

    """
    cdef int    i
    cdef int    right
    cdef int    counts[4096]
    cdef double prob[4096]
    cdef double bg[4096]
    cdef int    in_gene      = 0
    cdef int    path         = ipath
    cdef int    left         = -1
    cdef int    glob

    for i in range(4096):
        prob[i] = 0.0
        bg[i] = 0.0

    # get background counts
    # (shortened code from `calc_mer_bg` without malloc)
    glob = 0
    memset(counts, 0, 4096*sizeof(int))
    for i in range(seq.slen - 5):
        counts[seq._mer_ndx(i, length=6, strand=+1)] += 1
        counts[seq._mer_ndx(i, length=6, strand=-1)] += 1
        glob += 2
    for i in range(4096):
        bg[i] = (<double> counts[i]) / (<double> glob)

    # get counts in genes
    glob = 0
    memset(counts, 0, 4096*sizeof(int))
    while path != -1:
        if nodes.nodes[path].strand == 1:
            if nodes.nodes[path].type == node_type.STOP:
                in_gene = 1
                right = nodes.nodes[path].ndx+2
            elif in_gene == 1:
                left = nodes.nodes[path].ndx
                for i in range(left, right-5, 3):
                    counts[seq._mer_ndx(i, length=6, strand=1)] += 1
                    glob += 1
                in_gene = 0
        else:
            if nodes.nodes[path].type != node_type.STOP:
                in_gene = -1
                left = seq.slen - nodes.nodes[path].ndx - 1
            elif in_gene == -1:
                right = seq.slen - nodes.nodes[path].ndx + 1
                for i in range(left, right-5, 3):
                    counts[seq._mer_ndx(i, length=6, strand=-1)] += 1
                    glob += 1
                in_gene = 0
        path = nodes.nodes[path].traceb

    # compute log likelihood
    for i in range(4096):
        prob[i] = (<double> counts[i])/(<double> glob)
        if prob[i] == 0 and bg[i] != 0:
            tinf.tinf.gene_dc[i] = -5.0;
        elif bg[i] == 0:
            tinf.tinf.gene_dc[i] = 0.0
        else:
            tinf.tinf.gene_dc[i] = log(prob[i]/bg[i])
        if tinf.tinf.gene_dc[i] > 5.0:
            tinf.tinf.gene_dc[i] = 5.0
        elif tinf.tinf.gene_dc[i] < -5.0:
            tinf.tinf.gene_dc[i] = -5.0

cdef int* calc_most_gc_frame(Sequence seq) nogil except NULL:
    cdef int  i
    cdef int  j
    cdef int  win
    cdef int* fwd
    cdef int* bwd
    cdef int* tot
    cdef int* gp
    cdef int  slen = seq.slen

    gp = <int*> malloc(seq.slen*sizeof(int));
    fwd = <int*> malloc(seq.slen*sizeof(int));
    bwd = <int*> malloc(seq.slen*sizeof(int));
    tot = <int*> malloc(seq.slen*sizeof(int));
    if fwd == NULL or bwd == NULL or gp == NULL or tot == NULL:
        free(gp)
        free(bwd)
        free(tot)
        with gil:
            raise MemoryError("could not allocate GC frame buffers")

    for i in range(slen):
        fwd[i] = 0
        bwd[i] = 0
        tot[i] = 0
        gp[i] = -1

    for j in range(0, slen):
        if j < 3:
            fwd[j] = seq._is_gc(j)
            bwd[seq.slen-j-1] = seq._is_gc(seq.slen-j-1)
        else:
            fwd[j] = fwd[j-3] + seq._is_gc(j)
            bwd[slen-j-1] = bwd[slen-j+2] + seq._is_gc(slen-j-1)
    for i in range(slen):
        tot[i] = fwd[i] + bwd[i] - seq._is_gc(i);
        if i >= WINDOW//2:
            tot[i] -= fwd[i-WINDOW//2]
        if i + WINDOW//2 < seq.slen:
            tot[i] -= bwd[i+WINDOW//2]
    free(bwd);
    for i in range(0, slen-2, 3):
        win = sequence.max_fr(tot[i], tot[i+1], tot[i+2])
        for j in range(3):
            gp[i+j] = win
    free(tot);

    return gp;

cpdef int calc_orf_gc(Nodes nodes, Sequence seq, TrainingInfo tinf) nogil except -1:
    cdef int i
    cdef int j
    cdef int last[3]
    cdef int phase
    cdef double gc[3]
    cdef double gsize = 0.0

    # direct strand
    gc[0] = gc[1] = gc[2] = 0.0
    for i in reversed(range(<int> nodes.length)):
        if nodes.nodes[i].strand == 1:
            phase = nodes.nodes[i].ndx %3
            if nodes.nodes[i].type == node_type.STOP:
                last[phase] = j = nodes.nodes[i].ndx
                gc[phase] = seq._is_gc(j) + seq._is_gc(j+1) + seq._is_gc(j+2)
            else:
                for j in range(last[phase] - 3, nodes.nodes[i].ndx - 1, -3):
                    gc[phase] += seq._is_gc(j) + seq._is_gc(j+1) + seq._is_gc(j+2)
                gsize = abs(nodes.nodes[i].stop_val - nodes.nodes[i].ndx) + 3.0
                nodes.nodes[i].gc_cont = gc[phase] / gsize
                last[phase] = nodes.nodes[i].ndx

    # reverse strand
    gc[0] = gc[1] = gc[2] = 0.0
    for i in range(<int> nodes.length):
        if nodes.nodes[i].strand == -1:
            phase = nodes.nodes[i].ndx % 3
            if nodes.nodes[i].type == node_type.STOP:
                last[phase] = j = nodes.nodes[i].ndx
                gc[phase] = seq._is_gc(j) + seq._is_gc(j-1) + seq._is_gc(j-2)
            else:
                for j in range(last[phase] + 3, nodes.nodes[i].ndx + 1, 3):
                    gc[phase] += seq._is_gc(j) + seq._is_gc(j+1) + seq._is_gc(j+2)
                gsize = abs(nodes.nodes[i].stop_val - nodes.nodes[i].ndx) + 3.0
                nodes.nodes[i].gc_cont = gc[phase] / gsize
                last[phase] = nodes.nodes[i].ndx

cpdef void count_upstream_composition(Sequence seq, TrainingInfo tinf, int pos, int strand=1) nogil:
    cdef int start
    cdef int j
    cdef int i     = 0

    if strand == 1:
        for j in range(1, 3):
            if pos - j >= 0:
                tinf.tinf.ups_comp[i][seq.digits[pos-j] & 0b11] += 1
            i += 1
        for j in range(15, 45):
            if pos - j >= 0:
                tinf.tinf.ups_comp[i][seq.digits[pos-j] & 0b11] += 1
            i += 1
    else:
        start = seq.slen - 1 - pos
        for j in range(1, 3):
            if pos + j < seq.slen:
                tinf.tinf.ups_comp[i][_translation[seq.digits[pos+j]] & 0b11] += 1
            i += 1
        for j in range(15, 45):
            if pos + j < seq.slen:
                tinf.tinf.ups_comp[i][_translation[seq.digits[pos+j]] & 0b11] += 1
            i += 1

cpdef int dynamic_programming(Nodes nodes, TrainingInfo tinf, ConnectionScorer scorer, bint final=False) nogil:
    cdef int    i
    cdef int    j
    cdef int    min
    cdef int    path
    cdef int    nxt
    cdef int    tmp
    cdef int    max_ndx = -1
    cdef double max_sc  = -1.0

    if nodes.length == 0:
        return -1

    for i in range(<int> nodes.length):
        nodes.nodes[i].score = 0
        nodes.nodes[i].traceb = -1
        nodes.nodes[i].tracef = -1

    for i in range(<int> nodes.length):
        # Set up distance constraints for making connections,
        # but make exceptions for giant ORFS.
        min = 0 if i < dprog.MAX_NODE_DIST else i - dprog.MAX_NODE_DIST
        if nodes.nodes[i].strand == -1 and nodes.nodes[i].type != node_type.STOP and nodes.nodes[min].ndx >= nodes.nodes[i].stop_val:
            if nodes.nodes[i].ndx != nodes.nodes[i].stop_val:
                min = 0
        elif nodes.nodes[i].strand == 1 and nodes.nodes[i].type == node_type.STOP and nodes.nodes[min].ndx >= nodes.nodes[i].stop_val:
            if nodes.nodes[i].ndx != nodes.nodes[i].stop_val:
                min = 0
        min = 0 if min < dprog.MAX_NODE_DIST else min - dprog.MAX_NODE_DIST
        # Check which nodes can be skipped
        scorer._compute_skippable(min, i)
        # Score connections
        scorer._score_connections(nodes, min, i, tinf, final)

    for i in reversed(range(<int> nodes.length)):
        if nodes.nodes[i].strand == 1 and nodes.nodes[i].type != node_type.STOP:
            continue
        if nodes.nodes[i].strand == -1 and nodes.nodes[i].type == node_type.STOP:
            continue
        if nodes.nodes[i].score > max_sc:
            max_sc = nodes.nodes[i].score
            max_ndx = i

    # First Pass: untangle the triple overlaps
    path = max_ndx
    while nodes.nodes[path].traceb != -1:
        nxt = nodes.nodes[path].traceb
        if nodes.nodes[path].strand == -1 and nodes.nodes[path].type == node_type.STOP and nodes.nodes[nxt].strand == 1 and nodes.nodes[nxt].type == node_type.STOP and nodes.nodes[path].ov_mark != -1 and nodes.nodes[path].ndx > nodes.nodes[nxt].ndx:
            tmp = nodes.nodes[path].star_ptr[nodes.nodes[path].ov_mark];
            i = tmp
            while nodes.nodes[i].ndx != nodes.nodes[tmp].stop_val:
                i -= 1
            nodes.nodes[path].traceb = tmp
            nodes.nodes[tmp].traceb = i
            nodes.nodes[i].ov_mark = -1
            nodes.nodes[i].traceb = nxt
        path = nodes.nodes[path].traceb

    # Second Pass: Untangle the simple overlaps
    path = max_ndx
    while nodes.nodes[path].traceb != -1:
        nxt = nodes.nodes[path].traceb
        if nodes.nodes[path].strand == -1 and nodes.nodes[path].type != node_type.STOP and nodes.nodes[nxt].strand == 1 and nodes.nodes[nxt].type == node_type.STOP:
            i = path
            while nodes.nodes[i].ndx != nodes.nodes[path].stop_val:
                i -= 1
            nodes.nodes[path].traceb = i
            nodes.nodes[i].traceb = nxt
        if nodes.nodes[path].strand == 1 and nodes.nodes[path].type == node_type.STOP and nodes.nodes[nxt].strand == 1 and nodes.nodes[nxt].type == node_type.STOP:
            nodes.nodes[path].traceb = nodes.nodes[nxt].star_ptr[nodes.nodes[path].ndx%3]
            nodes.nodes[nodes.nodes[path].traceb].traceb = nxt
        if nodes.nodes[path].strand == -1 and nodes.nodes[path].type == node_type.STOP and nodes.nodes[nxt].strand == -1 and nodes.nodes[nxt].type == node_type.STOP:
            nodes.nodes[path].traceb = nodes.nodes[path].star_ptr[nodes.nodes[nxt].ndx%3]
            nodes.nodes[nodes.nodes[path].traceb].traceb = nxt
        path = nodes.nodes[path].traceb

    # Mark forward pointers
    path = max_ndx
    while nodes.nodes[path].traceb != -1:
        nodes.nodes[nodes.nodes[path].traceb].tracef = path
        path = nodes.nodes[path].traceb

    return -1 if nodes.nodes[max_ndx].traceb == -1 else max_ndx

cpdef int find_best_upstream_motif(Nodes nodes, int ni, Sequence seq, TrainingInfo tinf, int stage) nogil except -1:
    cdef int i
    cdef int j
    cdef int start
    cdef int spacer
    cdef int spacendx
    cdef int index
    cdef int max_spacer   = 0
    cdef int max_spacendx = 0
    cdef int max_len      = 0
    cdef int max_ndx      = 0
    cdef double max_sc    = -100.0
    cdef double score     = 0.0

    if nodes.nodes[ni].type == node_type.STOP or nodes.nodes[ni].edge:
        return 0

    if nodes.nodes[ni].strand == 1:
        start = nodes.nodes[ni].ndx
    else:
        start = seq.slen - 1 - nodes.nodes[ni].ndx

    for i in reversed(range(4)):
        for j in range(start-18-i, start-5-i):
            if j < 0:
                continue
            spacer = start - j - i - 3

            if j <= start - 16 - i:
                spacendx = 3
            elif j <= start - 14 - i:
                spacendx = 2
            elif j >= start - 7 - i:
                spacendx = 1
            else:
                spacendx = 0

            index = seq._mer_ndx(j, length=i+3, strand=nodes.nodes[ni].strand)
            score = tinf.tinf.mot_wt[i][spacendx][index]
            if score > max_sc:
                max_sc = score
                max_spacendx = spacendx
                max_spacer = spacer
                max_ndx = index
                max_len = i+3

    if stage == 2 and (max_sc == -4.0 or max_sc < tinf.tinf.no_mot + 0.69):
        nodes.nodes[ni].mot.ndx = 0
        nodes.nodes[ni].mot.len = 0
        nodes.nodes[ni].mot.spacendx = 0
        nodes.nodes[ni].mot.spacer = 0
        nodes.nodes[ni].mot.score = tinf.tinf.no_mot
    else:
        nodes.nodes[ni].mot.ndx = max_ndx
        nodes.nodes[ni].mot.len = max_len
        nodes.nodes[ni].mot.spacendx = max_spacendx
        nodes.nodes[ni].mot.spacer = max_spacer
        nodes.nodes[ni].mot.score = max_sc

cpdef void raw_coding_score(Nodes nodes, Sequence seq, TrainingInfo tinf) nogil:
    cdef double  score[3]
    cdef double  lfac
    cdef double  no_stop
    cdef double  gsize
    cdef ssize_t last[3]
    cdef size_t  phase
    cdef ssize_t j
    cdef ssize_t i
    cdef ssize_t nn = nodes.length

    if tinf.tinf.trans_table != 11:
        no_stop =  ((1-tinf.tinf.gc)*(1-tinf.tinf.gc)*tinf.tinf.gc)     / 8.0
        no_stop += ((1-tinf.tinf.gc)*(1-tinf.tinf.gc)*(1-tinf.tinf.gc)) / 8.0
        no_stop = 1 - no_stop
    else:
        no_stop =  ((1-tinf.tinf.gc)*(1-tinf.tinf.gc)*tinf.tinf.gc)     / 4.0
        no_stop += ((1-tinf.tinf.gc)*(1-tinf.tinf.gc)*(1-tinf.tinf.gc)) / 8.0
        no_stop = 1 - no_stop

    # Initial Pass: Score coding potential (start->stop)
    score[0] = score[1] = score[2] = 0.0
    for i in reversed(range(nn)):
        if nodes.nodes[i].strand == 1:
            phase = nodes.nodes[i].ndx%3
            if nodes.nodes[i].type == node_type.STOP:
                last[phase] = nodes.nodes[i].ndx
                score[phase] = 0.0
            else:
                for j in range(last[phase] - 3, nodes.nodes[i].ndx - 1, -3):
                    score[phase] += tinf.tinf.gene_dc[seq._mer_ndx(j, length=6, strand=1)];
                nodes.nodes[i].cscore = score[phase]
                last[phase] = nodes.nodes[i].ndx
    score[0] = score[1] = score[2] = 0.0
    for i in range(nn):
        if nodes.nodes[i].strand == -1:
            phase = nodes.nodes[i].ndx%3
            if nodes.nodes[i].type == node_type.STOP:
                last[phase] = nodes.nodes[i].ndx
                score[phase] = 0.0
            else:
                for j in range(last[phase] + 3, nodes.nodes[i].ndx + 1, 3):
                    score[phase] += tinf.tinf.gene_dc[seq._mer_ndx(seq.slen-1-j, length=6, strand=-1)]
                nodes.nodes[i].cscore = score[phase]
                last[phase] = nodes.nodes[i].ndx

    # Second Pass: Penalize start nodes with ascending coding to their left
    score[0] = score[1] = score[2] = -10000.0;
    for i in range(nn):
        if nodes.nodes[i].strand == 1:
            phase = nodes.nodes[i].ndx%3
            if nodes.nodes[i].type == node_type.STOP:
                score[phase] = -10000.0
            elif nodes.nodes[i].cscore > score[phase]:
                score[phase] = nodes.nodes[i].cscore
            else:
                nodes.nodes[i].cscore -= score[phase] - nodes.nodes[i].cscore
    score[0] = score[1] = score[2] = -10000.0
    for i in reversed(range(nn)):
        if nodes.nodes[i].strand == -1:
            phase = nodes.nodes[i].ndx%3
            if nodes.nodes[i].type == node_type.STOP:
                score[phase] = -10000.0
            elif nodes.nodes[i].cscore > score[phase]:
                score[phase] = nodes.nodes[i].cscore
            else:
                nodes.nodes[i].cscore -= (score[phase] - nodes.nodes[i].cscore);

    # Third Pass: Add length-based factor to the score
    # Penalize start nodes based on length to their left
    for i in range(nn):
        if nodes.nodes[i].strand == 1:
            phase = nodes.nodes[i].ndx%3
            if nodes.nodes[i].type == node_type.STOP:
                score[phase] = -10000.0;
            else:
                gsize = ((<double> nodes.nodes[i].stop_val - nodes.nodes[i].ndx)+3.0)/3.0
                if gsize > 1000.0:
                    lfac = log((1-pow(no_stop, 1000.0))/pow(no_stop, 1000.0))
                    lfac -= log((1-pow(no_stop, 80))/pow(no_stop, 80))
                    lfac *= (gsize - 80) / 920.0
                else:
                    lfac = log((1-pow(no_stop, gsize))/pow(no_stop, gsize))
                    lfac -= log((1-pow(no_stop, 80))/pow(no_stop, 80))
                if lfac > score[phase]:
                    score[phase] = lfac
                else:
                    lfac -= fmax(fmin(score[phase] - lfac, lfac), 0);
                if lfac > 3.0 and nodes.nodes[i].cscore < 0.5*lfac:
                    nodes.nodes[i].cscore = 0.5*lfac
                nodes.nodes[i].cscore += lfac
    for i in reversed(range(nn)):
        if nodes.nodes[i].strand == -1:
            phase = nodes.nodes[i].ndx%3;
            if nodes.nodes[i].type == node_type.STOP:
                score[phase] = -10000.0;
            else:
                gsize = ((<double> nodes.nodes[i].ndx - nodes.nodes[i].stop_val)+3.0)/3.0
                if(gsize > 1000.0):
                    lfac = log((1-pow(no_stop, 1000.0))/pow(no_stop, 1000.0));
                    lfac -= log((1-pow(no_stop, 80))/pow(no_stop, 80));
                    lfac *= (gsize - 80) / 920.0;
                else:
                    lfac = log((1-pow(no_stop, gsize))/pow(no_stop, gsize));
                    lfac -= log((1-pow(no_stop, 80))/pow(no_stop, 80));
                if lfac > score[phase]:
                    score[phase] = lfac
                else:
                    lfac -= fmax(fmin(score[phase] - lfac, lfac), 0);
                if lfac > 3.0 and nodes.nodes[i].cscore < 0.5*lfac:
                    nodes.nodes[i].cscore = 0.5*lfac
                nodes.nodes[i].cscore += lfac

cpdef void rbs_score(Nodes nodes, Sequence seq, TrainingInfo tinf) nogil:
    cdef int i
    cdef int j
    cdef int cur_sc[2]
    cdef int slen      = seq.slen

    for i in range(<int> nodes.length):
        if nodes.nodes[i].type == node_type.STOP or nodes.nodes[i].edge:
            continue
        nodes.nodes[i].rbs[0] = nodes.nodes[i].rbs[1] = 0

        if nodes.nodes[i].strand == 1:
            for j in range(nodes.nodes[i].ndx - 20, nodes.nodes[i].ndx - 5):
                if j < 0:
                    continue
                cur_sc[0] = shine_dalgarno_exact(seq, j, nodes.nodes[i].ndx, tinf, strand=1)
                cur_sc[1] = shine_dalgarno_mm(seq, j, nodes.nodes[i].ndx, tinf, strand=1)
                if cur_sc[0] > nodes.nodes[i].rbs[0]:
                    nodes.nodes[i].rbs[0] = cur_sc[0]
                if cur_sc[1] > nodes.nodes[i].rbs[1]:
                    nodes.nodes[i].rbs[1] = cur_sc[1]
        else:
            for j in range(slen - nodes.nodes[i].ndx - 21, slen - nodes.nodes[i].ndx - 6):
                if j + 1 > slen:
                    continue
                cur_sc[0] = shine_dalgarno_exact(seq, j, slen-1-nodes.nodes[i].ndx, tinf, strand=-1)
                cur_sc[1] = shine_dalgarno_mm(seq, j, slen-1-nodes.nodes[i].ndx, tinf, strand=-1)
                if cur_sc[0] > nodes.nodes[i].rbs[0]:
                    nodes.nodes[i].rbs[0] = cur_sc[0]
                if cur_sc[1] > nodes.nodes[i].rbs[1]:
                    nodes.nodes[i].rbs[1] = cur_sc[1]

cpdef void score_nodes(Nodes nodes, Sequence seq, TrainingInfo tinf, bint closed=False, bint is_meta=False) nogil:
    """score_nodes(nodes, seq, tinf, closed=False, is_meta=False)\n--

    Score the start nodes currently scored in ``nodes``.

    Note:
        This function is reimplemented from the ``score_nodes`` function of
        ``node.c`` from the Prodigal source, and already contains the patch
        from `hyattpd/Prodigal#88 <https://github.com/hyattpd/Prodigal/pull/88>`_.

    """
    cdef size_t i
    cdef size_t j
    cdef size_t orf_length
    cdef double negf
    cdef double posf
    cdef double rbs1
    cdef double rbs2
    cdef double sd_score
    cdef double edge_gene
    cdef double min_meta_len

    # Calculate raw coding potential for every start-stop pair
    calc_orf_gc(nodes, seq, tinf)
    raw_coding_score(nodes, seq, tinf)

    # Calculate raw RBS Scores for every start node.
    if tinf.tinf.uses_sd:
        rbs_score(nodes, seq, tinf)
    else:
        for i in range(nodes.length):
            if nodes.nodes[i].type == node_type.STOP or nodes.nodes[i].edge:
                continue
            find_best_upstream_motif(nodes, i, seq, tinf, stage=2)

    # Score the start nodes
    for i in range(nodes.length):
        if nodes.nodes[i].type == node_type.STOP:
            continue

        # compute ORF length, we'll need it later several times
        if nodes.nodes[i].ndx > nodes.nodes[i].stop_val:
            orf_length = nodes.nodes[i].ndx - nodes.nodes[i].stop_val
        else:
            orf_length = nodes.nodes[i].stop_val - nodes.nodes[i].ndx

        # Does this gene run off the edge?
        edge_gene = 0
        if nodes.nodes[i].edge:
            edge_gene += 1
        if (
                nodes.nodes[i].strand == 1 and not seq._is_stop(nodes.nodes[i].stop_val, tinf.tinf.trans_table)
            or  nodes.nodes[i].strand == -1 and not seq._is_stop(seq.slen - 1 - nodes.nodes[i].stop_val, tinf.tinf.trans_table, strand=-1)
        ):
            edge_gene += 1

        # Edge Nodes : stops with no starts, give a small bonus
        if nodes.nodes[i].edge:
            nodes.nodes[i].tscore = node.EDGE_BONUS*tinf.tinf.st_wt / edge_gene
            nodes.nodes[i].uscore = 0.0
            nodes.nodes[i].rscore = 0.0
        else:
            # Type Score
            nodes.nodes[i].tscore = tinf.tinf.type_wt[nodes.nodes[i].type] * tinf.tinf.st_wt

            # RBS Motif Score
            rbs1 = tinf.tinf.rbs_wt[nodes.nodes[i].rbs[0]];
            rbs2 = tinf.tinf.rbs_wt[nodes.nodes[i].rbs[1]];
            sd_score = node.dmax(rbs1, rbs2) * tinf.tinf.st_wt
            if tinf.tinf.uses_sd:
                nodes.nodes[i].rscore = sd_score
            else:
                nodes.nodes[i].rscore = tinf.tinf.st_wt * nodes.nodes[i].mot.score
                if nodes.nodes[i].rscore < sd_score and tinf.tinf.no_mot > -0.5:
                    nodes.nodes[i].rscore = sd_score

            # Upstream Score
            score_upstream_composition(nodes, i, seq, tinf)

            # Penalize upstream score if choosing this start would stop the gene
            # from running off the edge
            if not closed and nodes.nodes[i].ndx <= 2 and nodes.nodes[i].strand == 1:
                nodes.nodes[i].uscore += node.EDGE_UPS*tinf.tinf.st_wt
            elif not closed and nodes.nodes[i].ndx >= seq.slen - 3 and nodes.nodes[i].strand == -1:
                nodes.nodes[i].uscore += node.EDGE_UPS*tinf.tinf.st_wt
            elif i < 500 and nodes.nodes[i].strand == 1:
                for j in reversed(range(i)):
                    if nodes.nodes[j].edge and nodes.nodes[i].stop_val == nodes.nodes[j].stop_val:
                        nodes.nodes[i].uscore += node.EDGE_UPS*tinf.tinf.st_wt
                        break
            elif i + 500>= nodes.length and nodes.nodes[i].strand == -1:
                for j in range(i+1, nodes.length):
                    if nodes.nodes[j].edge and nodes.nodes[i].stop_val == nodes.nodes[j].stop_val:
                        nodes.nodes[i].uscore += node.EDGE_UPS*tinf.tinf.st_wt
                        break

        # Convert starts at base 1 and slen to edge genes if closed = 0
        if (
                not closed
            and not nodes.nodes[i].edge
            and ((nodes.nodes[i].ndx <= 2 and nodes.nodes[i].strand == 1) or (nodes.nodes[i].ndx >= seq.slen - 3 and nodes.nodes[i].strand == -1))
        ):
            edge_gene += 1
            nodes.nodes[i].edge = True
            nodes.nodes[i].tscore = 0.0
            nodes.nodes[i].uscore = node.EDGE_BONUS*tinf.tinf.st_wt / edge_gene
            nodes.nodes[i].rscore = 0.0

        # Penalize starts with no stop codon
        if not nodes.nodes[i].edge and edge_gene == 1:
            nodes.nodes[i].uscore -= 0.5*node.EDGE_BONUS*tinf.tinf.st_wt

        # Penalize non-edge genes < 250bp
        if edge_gene == 0 and orf_length < 250:
            negf = 250.0 / <float> orf_length
            posf = <float> orf_length / 250.0
            nodes.nodes[i].rscore *= negf if nodes.nodes[i].rscore < 0 else posf
            nodes.nodes[i].uscore *= negf if nodes.nodes[i].uscore < 0 else posf
            nodes.nodes[i].tscore *= negf if nodes.nodes[i].tscore < 0 else posf

        # Coding Penalization in Metagenomic Fragments:  Internal genes must
        # have a score of 5.0 and be >= 120bp.  High GC genes are also penalized.                                  */
        if (
                is_meta
            and seq.slen < 3000
            and edge_gene == 0
            and (nodes.nodes[i].cscore < 5.0 or orf_length < 120)
        ):
            nodes.nodes[i].cscore -= node.META_PEN * node.dmax(0, (3000.0 - seq.slen) / 2700.0)

        # Base Start Score
        nodes.nodes[i].sscore = nodes.nodes[i].tscore + nodes.nodes[i].rscore + nodes.nodes[i].uscore

        # Penalize starts if coding is negative.  Larger penalty for edge genes,
        # since the start is offset by a smaller amount of coding than normal.
        if nodes.nodes[i].cscore < 0.0:
            if edge_gene > 0 and not nodes.nodes[i].edge:
                if not is_meta or seq.slen > 1500:
                    nodes.nodes[i].sscore -= tinf.tinf.st_wt
                else:
                    nodes.nodes[i].sscore -= 10.31 - 0.004*seq.slen
            elif is_meta and seq.slen < 3000 and nodes.nodes[i].edge:
                min_meta_len = sqrt(<double> seq.slen)*5.0
                if orf_length >= min_meta_len:
                    if nodes.nodes[i].cscore >= 0:
                        nodes.nodes[i].cscore = -1.0
                    nodes.nodes[i].sscore = 0.0
                    nodes.nodes[i].uscore = 0.0
            else:
                nodes.nodes[i].sscore -= 0.5
        elif (
                is_meta
            and nodes.nodes[i].cscore < 5.0
            and orf_length < 120
            and nodes.nodes[i].sscore < 0.0
        ):
            nodes.nodes[i].sscore -= tinf.tinf.st_wt

cpdef void score_upstream_composition(Nodes nodes, int ni, Sequence seq, TrainingInfo tinf) nogil:
    cdef int i
    cdef int start
    cdef int count = 0
    cdef int mer
    cdef int strand

    if nodes.nodes[ni].strand == 1:
        start = nodes.nodes[ni].ndx
        strand = 1
    else:
        start = seq.slen - 1 - nodes.nodes[ni].ndx
        strand = -1

    nodes.nodes[ni].uscore = 0.0;
    for i in range(1, 45):
        if i > 2 and i < 15:
            continue
        if start < i:
            continue

        if strand == 1:
            mer = seq.digits[start - i] & 0b11
        else:
            mer = _translation[seq.digits[seq.slen - 1 - start + i]] & 0b11

        nodes.nodes[ni].uscore += 0.4 * tinf.tinf.st_wt * tinf.tinf.ups_comp[count][mer]
        count += 1

cpdef int shine_dalgarno_exact(Sequence seq, int pos, int start, TrainingInfo tinf, int strand=1) nogil:
    cdef int i
    cdef int j
    cdef int k
    cdef int mism
    cdef int rdis
    cdef int limit
    cdef int max_val
    cdef int cmp_val
    cdef int cur_val = 0
    cdef int match[6]
    cdef int cur_ctr
    cdef int dis_flag

    limit = start - 4 - pos
    if limit > 6:
        limit = 6
    for i in range(limit, 6):
        match[i] = -10

    # Compare the 6-base region to AGGAGG
    for i in range(limit):
        if pos + i < 0:
            continue
        if i%3 == 0 and seq._is_a(pos+i, strand=strand):
            match[i] = 2
        elif i%3 != 0 and seq._is_g(pos+i, strand=strand):
            match[i] = 3
        else:
            match[i] = -10

    # Find the maximally scoring motif
    max_val = 0
    for i in range(limit, 2, -1):
        for j in range(limit+1-i):
            cur_ctr = -2
            mism = 0;
            for k in range(j, j+i):
                cur_ctr += match[k]
                if match[k] < 0:
                    mism += 1
            if mism > 0:
                continue
            rdis = start - (pos + j + i)
            if rdis < 5 and i < 5:
                dis_flag = 2
            elif rdis < 5 and i >= 5:
                dis_flag = 1
            elif rdis > 10 and rdis <= 12 and i < 5:
                dis_flag = 1
            elif rdis > 10 and rdis <= 12 and i >= 5:
                dis_flag = 2
            elif rdis >= 13:
                dis_flag = 3
            else:
                dis_flag = 0
            if rdis > 15 or cur_ctr < 6.0:
                continue

            # Exact-Matching RBS Motifs
            if cur_ctr < 6:
                cur_val = 0
            elif cur_ctr == 6 and dis_flag == 2:
                cur_val = 1
            elif cur_ctr == 6 and dis_flag == 3:
                cur_val = 2
            elif cur_ctr == 8 and dis_flag == 3:
                cur_val = 3
            elif cur_ctr == 9 and dis_flag == 3:
                cur_val = 3
            elif cur_ctr == 6 and dis_flag == 1:
                cur_val = 6
            elif cur_ctr == 11 and dis_flag == 3:
                cur_val = 10
            elif cur_ctr == 12 and dis_flag == 3:
                cur_val = 10
            elif cur_ctr == 14 and dis_flag == 3:
                cur_val = 10
            elif cur_ctr == 8 and dis_flag == 2:
                cur_val = 11
            elif cur_ctr == 9 and dis_flag == 2:
                cur_val = 11
            elif cur_ctr == 8 and dis_flag == 1:
                cur_val = 12
            elif cur_ctr == 9 and dis_flag == 1:
                cur_val = 12
            elif cur_ctr == 6 and dis_flag == 0:
                cur_val = 13
            elif cur_ctr == 8 and dis_flag == 0:
                cur_val = 15
            elif cur_ctr == 9 and dis_flag == 0:
                cur_val = 16
            elif cur_ctr == 11 and dis_flag == 2:
                cur_val = 20
            elif cur_ctr == 11 and dis_flag == 1:
                cur_val = 21
            elif cur_ctr == 11 and dis_flag == 0:
                cur_val = 22
            elif cur_ctr == 12 and dis_flag == 2:
                cur_val = 20
            elif cur_ctr == 12 and dis_flag == 1:
                cur_val = 23
            elif cur_ctr == 12 and dis_flag == 0:
                cur_val = 24
            elif cur_ctr == 14 and dis_flag == 2:
                cur_val = 25
            elif cur_ctr == 14 and dis_flag == 1:
                cur_val = 26
            elif cur_ctr == 14 and dis_flag == 0:
                cur_val = 27

            if tinf.tinf.rbs_wt[cur_val] < tinf.tinf.rbs_wt[max_val]:
                continue
            if tinf.tinf.rbs_wt[cur_val] == tinf.tinf.rbs_wt[max_val] and cur_val < max_val:
                continue
            max_val = cur_val

    return max_val

cpdef int shine_dalgarno_mm(Sequence seq, int pos, int start, TrainingInfo tinf, int strand=1) nogil:
    cdef int i
    cdef int j
    cdef int k
    cdef int mism
    cdef int rdis
    cdef int limit
    cdef int max_val
    cdef int cmp_val
    cdef int cur_val = 0
    cdef int match[6]
    cdef int cur_ctr
    cdef int dis_flag

    limit = start - 4 - pos
    if limit > 6:
        limit = 6
    for i in range(limit, 6):
        match[i] = -10

    # Compare the 6-base region to AGGAGG
    for i in range(limit):
        if pos + i < 0:
            continue
        if i%3 == 0:
            match[i] = 2 if seq._is_a(pos+i, strand=strand) else -3
        else:
            match[i] = 3 if seq._is_g(pos+i, strand=strand) else -2

    # Find the maximally scoring motif
    max_val = 0
    for i in range(limit, 4, -1):
        for j in range(limit+1-i):
            cur_ctr = -2
            mism = 0;
            for k in range(j, j+i):
                cur_ctr += match[k]
                if match[k] < 0.0:
                    mism += 1
                    if k <= j+1 or k >= j+i-2:
                      cur_ctr -= 10
            if mism != 1:
                continue
            rdis = start - (pos + j + i)
            if rdis < 5:
                dis_flag = 1
            elif rdis > 10 and rdis <= 12:
                dis_flag = 2
            elif rdis >= 13:
                dis_flag = 3
            else:
                dis_flag = 0
            if rdis > 15 or cur_ctr < 6:
                continue

            # Single-Matching RBS Motifs
            if cur_ctr < 6:
                cur_val = 0
            elif cur_ctr == 6 and dis_flag == 3:
                cur_val = 2
            elif cur_ctr == 7 and dis_flag == 3:
                cur_val = 2
            elif cur_ctr == 9 and dis_flag == 3:
                cur_val = 3
            elif cur_ctr == 6 and dis_flag == 2:
                cur_val = 4
            elif cur_ctr == 6 and dis_flag == 1:
                cur_val = 5
            elif cur_ctr == 6 and dis_flag == 0:
                cur_val = 9
            elif cur_ctr == 7 and dis_flag == 2:
                cur_val = 7
            elif cur_ctr == 7 and dis_flag == 1:
                cur_val = 8
            elif cur_ctr == 7 and dis_flag == 0:
                cur_val = 14
            elif cur_ctr == 9 and dis_flag == 2:
                cur_val = 17
            elif cur_ctr == 9 and dis_flag == 1:
                cur_val = 18
            elif cur_ctr == 9 and dis_flag == 0:
                cur_val = 19

            if tinf.tinf.rbs_wt[cur_val] < tinf.tinf.rbs_wt[max_val]:
                continue
            if tinf.tinf.rbs_wt[cur_val] == tinf.tinf.rbs_wt[max_val] and cur_val < max_val:
                continue
            max_val = cur_val

    return max_val

cpdef void train_starts_sd(Nodes nodes, Sequence seq, TrainingInfo tinf) nogil:
    cdef int phase
    cdef int rbs[3]
    cdef int type[3]
    cdef int bndx[3]
    cdef int max_rb
    cdef double sum
    cdef double rbg[28]
    cdef double rreal[28]
    cdef double best[3]
    cdef double tbg[3]
    cdef double treal[3]
    cdef double sthresh   = 35.0
    cdef double wt        = tinf.tinf.st_wt

    cdef ssize_t i
    cdef ssize_t j
    cdef ssize_t nn = nodes.length

    # reset training info
    for j in range(3):
        tinf.tinf.type_wt[j] = 0.0
    for j in range(28):
        tinf.tinf.rbs_wt[j] = 0.0
    for i in range(32):
        for j in range(4):
            tinf.tinf.ups_comp[i][j] = 0.0

    # Build the background of random types
    for i in range(3):
        tbg[i] = 0.0
    for i in range(nn):
        if nodes.nodes[i].type == node_type.STOP:
            continue
        tbg[nodes.nodes[i].type] += 1.0
    sum = 0.0
    for i in range(3):
        sum += tbg[i]
    for i in range(3):
        tbg[i] /= sum

    # Iterate 10 times through the list of nodes
    # Converge upon optimal weights for ATG vs GTG vs TTG and RBS motifs
    # (convergence typically takes 4-5 iterations, but we run a few
    # extra to be safe)
    for i in range(10):
        # Recalculate the RBS motif background */
        for j in range(28):
            rbg[j] = 0.0
        for j in range(nn):
            if nodes.nodes[j].type == node_type.STOP or nodes.nodes[j].edge:
                continue
            if tinf.tinf.rbs_wt[nodes.nodes[j].rbs[0]] > tinf.tinf.rbs_wt[nodes.nodes[j].rbs[1]]+1.0 or nodes.nodes[j].rbs[1] == 0:
                max_rb = nodes.nodes[j].rbs[0]
            elif tinf.tinf.rbs_wt[nodes.nodes[j].rbs[0]] < tinf.tinf.rbs_wt[nodes.nodes[j].rbs[1]]-1.0 or nodes.nodes[j].rbs[0] == 0:
                max_rb = nodes.nodes[j].rbs[1]
            elif nodes.nodes[j].rbs[0] > nodes.nodes[j].rbs[1]:
                max_rb = nodes.nodes[j].rbs[0]
            else:
                max_rb = nodes.nodes[j].rbs[1]
            rbg[max_rb] += 1.0;

        sum = 0.0
        for j in range(28):
            sum += rbg[j]
        for j in range(28):
            rbg[j] /= sum

        for j in range(28):
            rreal[j] = 0.0
        for j in range(3):
            treal[j] = 0.0

        # Forward strand pass
        for j in range(3):
            best[j] = 0.0; bndx[j] = -1; rbs[j] = 0; type[j] = 0;
        for j in range(nn):
            if nodes.nodes[j].type != node_type.STOP and nodes.nodes[j].edge:
                continue
            if nodes.nodes[j].strand == 1:
                phase = nodes.nodes[j].ndx % 3
                if nodes.nodes[j].type == node_type.STOP:
                    if best[phase] >= sthresh and nodes.nodes[bndx[phase]].ndx%3 == phase:
                        rreal[rbs[phase]] += 1.0
                        treal[type[phase]] += 1.0
                        if i == 9:
                            count_upstream_composition(seq, tinf, nodes.nodes[bndx[phase]].ndx, strand=1)
                    best[phase] = 0.0; bndx[phase] = -1; rbs[phase] = 0; type[phase] = 0;
                else:
                    if tinf.tinf.rbs_wt[nodes.nodes[j].rbs[0]] > tinf.tinf.rbs_wt[nodes.nodes[j].rbs[1]]+1.0 or nodes.nodes[j].rbs[1] == 0:
                        max_rb = nodes.nodes[j].rbs[0]
                    elif tinf.tinf.rbs_wt[nodes.nodes[j].rbs[0]] < tinf.tinf.rbs_wt[nodes.nodes[j].rbs[1]]-1.0 or nodes.nodes[j].rbs[0] == 0:
                        max_rb = nodes.nodes[j].rbs[1]
                    elif nodes.nodes[j].rbs[0] > nodes.nodes[j].rbs[1]:
                        max_rb = nodes.nodes[j].rbs[0]
                    else:
                        max_rb = nodes.nodes[j].rbs[1]
                    if nodes.nodes[j].cscore + wt*tinf.tinf.rbs_wt[max_rb] + wt*tinf.tinf.type_wt[nodes.nodes[j].type] >= best[phase]:
                        best[phase] = nodes.nodes[j].cscore + wt*tinf.tinf.rbs_wt[max_rb] + wt*tinf.tinf.type_wt[nodes.nodes[j].type]
                        bndx[phase] = j
                        type[phase] = nodes.nodes[j].type
                        rbs[phase] = max_rb

        # Reverse strand pass
        for j in range(3):
            best[j] = 0.0
            bndx[j] = -1
            rbs[j] = 0
            type[j] = 0
        for j in reversed(range(nn)):
            if nodes.nodes[j].type != node_type.STOP and nodes.nodes[j].edge:
                continue
            if nodes.nodes[j].strand == -1:
                phase =  nodes.nodes[j].ndx % 3
                if nodes.nodes[j].type == node_type.STOP:
                    if best[phase] >= sthresh and nodes.nodes[bndx[phase]].ndx%3 == phase:
                        rreal[rbs[phase]] += 1.0
                        treal[type[phase]] += 1.0
                        if i == 9:
                            count_upstream_composition(seq, tinf, nodes.nodes[bndx[phase]].ndx, strand=-1);
                    best[phase] = 0.0
                    bndx[phase] = -1
                    rbs[phase] = 0
                    type[phase] = 0
                else:
                    if tinf.tinf.rbs_wt[nodes.nodes[j].rbs[0]] > tinf.tinf.rbs_wt[nodes.nodes[j].rbs[1]] + 1.0 or nodes.nodes[j].rbs[1] == 0:
                        max_rb = nodes.nodes[j].rbs[0]
                    elif tinf.tinf.rbs_wt[nodes.nodes[j].rbs[0]] < tinf.tinf.rbs_wt[nodes.nodes[j].rbs[1]] - 1.0 or nodes.nodes[j].rbs[0] == 0:
                        max_rb = nodes.nodes[j].rbs[1]
                    elif nodes.nodes[j].rbs[0] > nodes.nodes[j].rbs[1]:
                        max_rb = nodes.nodes[j].rbs[0]
                    else:
                        max_rb = nodes.nodes[j].rbs[1]
                    if nodes.nodes[j].cscore + wt*tinf.tinf.rbs_wt[max_rb] + wt*tinf.tinf.type_wt[nodes.nodes[j].type] >= best[phase]:
                        best[phase] = nodes.nodes[j].cscore + wt*tinf.tinf.rbs_wt[max_rb] + wt*tinf.tinf.type_wt[nodes.nodes[j].type]
                        bndx[phase] = j
                        type[phase] = nodes.nodes[j].type
                        rbs[phase] = max_rb

        # Update RBS weights
        sum = 0.0;
        for j in range(28):
            sum += rreal[j]
        if sum == 0.0:
            for j in range(28):
                tinf.tinf.rbs_wt[j] = 0.0
        else:
            for j in range(28):
                rreal[j] /= sum
                if rbg[j] != 0:
                    tinf.tinf.rbs_wt[j] = log(rreal[j]/rbg[j])
                else:
                    tinf.tinf.rbs_wt[j] = -4.0
                if tinf.tinf.rbs_wt[j] > 4.0:
                    tinf.tinf.rbs_wt[j] = 4.0
                elif tinf.tinf.rbs_wt[j] < -4.0:
                    tinf.tinf.rbs_wt[j] = -4.0

        # Update type weights
        sum = 0.0
        for j in range(3):
            sum += treal[j]
        if sum == 0.0:
            for j in range(3):
                tinf.tinf.type_wt[j] = 0.0
        else:
            for j in range(3):
                treal[j] /= sum;
                if tbg[j] != 0:
                    tinf.tinf.type_wt[j] = log(treal[j]/tbg[j])
                else:
                    tinf.tinf.type_wt[j] = -4.0
                if tinf.tinf.type_wt[j] > 4.0:
                    tinf.tinf.type_wt[j] = 4.0
                elif tinf.tinf.type_wt[j] < -4.0:
                    tinf.tinf.type_wt[j] = -4.0;
        if sum*2000.0 <= nodes.length:
            sthresh /= 2.0


    # Convert upstream base composition to a log score
    for i in range(32):
        sum = 0.0;
        for j in range(4):
            sum += tinf.tinf.ups_comp[i][j];
        if sum == 0.0:
            for j in range(4):
                tinf.tinf.ups_comp[i][j] = 0.0;
        else:
          for j in range(4):
              tinf.tinf.ups_comp[i][j] /= sum
              if tinf.tinf.gc <= 0.1:
                  if j == 0 or j == 3:
                      tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/0.90)
                  else:
                      tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/0.10)
              elif tinf.tinf.gc >= 0.9:
                  if j == 0 or j == 3:
                      tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/0.10)
                  else:
                      tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/0.90)
              else:
                  if j == 0 or j == 3:
                      tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/(1.0-tinf.tinf.gc))
                  else:
                      tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/tinf.tinf.gc)
              if tinf.tinf.ups_comp[i][j] > 4.0:
                  tinf.tinf.ups_comp[i][j] = 4.0
              if tinf.tinf.ups_comp[i][j] < -4.0:
                  tinf.tinf.ups_comp[i][j] = -4.0

cpdef void train_starts_nonsd(Nodes nodes, Sequence seq, TrainingInfo tinf) nogil:
    cdef int i
    cdef int j
    cdef int k
    cdef int l
    cdef int fr
    cdef int stage
    cdef int bndx[3]
    cdef int mgood[4][4][4096]
    cdef double sum
    cdef double ngenes
    cdef double zbg
    cdef double zreal
    cdef double best[3]
    cdef double tbg[3]
    cdef double treal[3]
    cdef double mbg[4][4][4096]
    cdef double mreal[4][4][4096]
    cdef double wt                = tinf.tinf.st_wt
    cdef double sthresh           = 35.0;
    cdef int    nn                = nodes.length

    for i in range(32):
        for j in range(4):
            tinf.tinf.ups_comp[i][j] = 0.0

    # Build the background of random types
    for i in range(3):
        tinf.tinf.type_wt[i] = 0.0
        tbg[i] = 0.0;
    for i in range(nn):
        if nodes.nodes[i].type == node_type.STOP:
            continue
        tbg[nodes.nodes[i].type] += 1.0
    sum = 0.0
    for i in range(3):
        sum += tbg[i]
    for i in range(3):
        tbg[i] /= sum

    # Iterate 20 times through the list of nodes
    # Converge upon optimal weights for ATG vs GTG vs TTG and RBS motifs
    # (convergence typically takes 4-5 iterations, but we run a few
    # extra to be safe)
    for i in range(20):
        # Determine which stage of motif finding we're in */
        if i < 4:
            stage = 0
        elif i < 12:
            stage = 1
        else:
            stage = 2

        # Recalculate the upstream motif background and set 'real' counts to 0
        for j in range(4):
            for k in range(4):
                for l in range(4096):
                    mbg[j][k][l] = 0.0
        zbg = 0.0
        for j in range(nn):
            if nodes.nodes[j].type == node_type.STOP or nodes.nodes[j].edge:
                continue
            find_best_upstream_motif(nodes, j, seq, tinf, stage);
            update_motif_counts(mbg, &zbg, seq, &nodes.nodes[j], stage)
        sum = 0.0
        for j in range(4):
            for k in range(4):
                for l in range(4096):
                    sum += mbg[j][k][l]
        sum += zbg
        for j in range(4):
            for k in range(4):
                for l in range(4096):
                    mbg[j][k][l] /= sum
        zbg /= sum

        # Reset counts of 'real' motifs/types to 0
        for j in range(4):
            for k in range(4):
                for l in range(4096):
                    mreal[j][k][l] = 0.0
        zreal = 0.0
        for j in range(3):
            treal[j] = 0.0
        ngenes = 0.0

        # Forward strand pass
        for j in range(3):
            best[j] = 0.0
            bndx[j] = -1
        for j in range(nn):
            if nodes.nodes[j].type != node_type.STOP and nodes.nodes[j].edge:
                continue
            if nodes.nodes[j].strand == 1:
                fr = (nodes.nodes[j].ndx)%3;
                if nodes.nodes[j].type == node_type.STOP:
                    if best[fr] >= sthresh:
                        ngenes += 1.0
                        treal[nodes.nodes[bndx[fr]].type] += 1.0
                        update_motif_counts(mreal, &zreal, seq, &nodes.nodes[bndx[fr]], stage)
                        if i == 19:
                            count_upstream_composition(seq, tinf, nodes.nodes[bndx[fr]].ndx, strand=1)
                    best[fr] = 0.0;
                    bndx[fr] = -1;
                else:
                    if nodes.nodes[j].cscore + wt*nodes.nodes[j].mot.score + wt*tinf.tinf.type_wt[nodes.nodes[j].type] >= best[fr]:
                        best[fr] = nodes.nodes[j].cscore + wt*nodes.nodes[j].mot.score + wt*tinf.tinf.type_wt[nodes.nodes[j].type]
                        bndx[fr] = j

        # Reverse strand pass
        for j in range(3):
            best[j] = 0.0
            bndx[j] = -1
        for j in reversed(range(nn)):
            if nodes.nodes[j].type != node_type.STOP and nodes.nodes[j].edge:
                continue
            if nodes.nodes[j].strand == -1:
                fr = (nodes.nodes[j].ndx)%3;
                if nodes.nodes[j].type == node_type.STOP:
                    if best[fr] >= sthresh:
                        ngenes += 1.0;
                        treal[nodes.nodes[bndx[fr]].type] += 1.0;
                        update_motif_counts(mreal, &zreal, seq, &nodes.nodes[bndx[fr]], stage)
                        if i == 19:
                            count_upstream_composition(seq, tinf, nodes.nodes[bndx[fr]].ndx, strand=-1)
                    best[fr] = 0.0
                    bndx[fr] = -1
                else:
                    if nodes.nodes[j].cscore + wt*nodes.nodes[j].mot.score + wt*tinf.tinf.type_wt[nodes.nodes[j].type] >= best[fr]:
                        best[fr] = nodes.nodes[j].cscore + wt*nodes.nodes[j].mot.score + wt*tinf.tinf.type_wt[nodes.nodes[j].type]
                        bndx[fr] = j

        # Update the log likelihood weights for type and RBS motifs
        if stage < 2:
            node.build_coverage_map(mreal, mgood, ngenes, stage)
        sum = 0.0
        for j in range(4):
            for k in range(4):
                for l in range(4096):
                    sum += mreal[j][k][l]
        sum += zreal;
        if sum == 0.0:
            for j in range(4):
                for k in range(4):
                    for l in range(4096):
                        tinf.tinf.mot_wt[j][k][l] = 0.0
            tinf.tinf.no_mot = 0.0;
        else:
            for j in range(4):
                for k in range(4):
                    for l in range(4096):
                        if mgood[j][k][l] == 0:
                            zreal += mreal[j][k][l]
                            zbg += mreal[j][k][l]
                            mreal[j][k][l] = 0.0
                            mbg[j][k][l] = 0.0
                        mreal[j][k][l] /= sum;
                        if mbg[j][k][l] != 0:
                            tinf.tinf.mot_wt[j][k][l] = log(mreal[j][k][l]/mbg[j][k][l])
                        else:
                            tinf.tinf.mot_wt[j][k][l] = -4.0
                        if tinf.tinf.mot_wt[j][k][l] > 4.0:
                            tinf.tinf.mot_wt[j][k][l] = 4.0
                        elif tinf.tinf.mot_wt[j][k][l] < -4.0:
                            tinf.tinf.mot_wt[j][k][l] = -4.0
        zreal /= sum
        if zbg != 0:
            tinf.tinf.no_mot = log(zreal/zbg)
        else:
            tinf.tinf.no_mot = -4.0
        if tinf.tinf.no_mot > 4.0:
            tinf.tinf.no_mot = 4.0
        elif tinf.tinf.no_mot < -4.0:
            tinf.tinf.no_mot = -4.0
        sum = 0.0;
        for j in range(3):
            sum += treal[j]
        if sum == 0.0:
            for j in range(3):
                tinf.tinf.type_wt[j] = 0.0
        else:
            for j in range(3):
                treal[j] /= sum;
                if tbg[j] != 0:
                    tinf.tinf.type_wt[j] = log(treal[j]/tbg[j])
                else:
                    tinf.tinf.type_wt[j] = -4.0;
                if tinf.tinf.type_wt[j] > 4.0:
                    tinf.tinf.type_wt[j] = 4.0
                elif tinf.tinf.type_wt[j] < -4.0:
                    tinf.tinf.type_wt[j] = -4.0
        if sum * 2000.0 <= nn:
            sthresh /= 2.0

    # Convert upstream base composition to a log score
    for i in range(32):
        sum = 0.0
        for j in range(4):
            sum += tinf.tinf.ups_comp[i][j]
        if sum == 0.0:
            for j in range(4):
                tinf.tinf.ups_comp[i][j] = 0.0
        else:
            for j in range(4):
                tinf.tinf.ups_comp[i][j] /= sum
                if tinf.tinf.gc <= 0.1:
                    if j == 0 or j == 3:
                        tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/0.90)
                    else:
                        tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/0.10)
                elif tinf.tinf.gc >= 0.9:
                    if j == 0 or j == 3:
                        tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/0.10)
                    else:
                        tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/0.90)
                else:
                    if j == 0 or j == 3:
                        tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/(1.0-tinf.tinf.gc));
                    else:
                        tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/tinf.tinf.gc);
                if tinf.tinf.ups_comp[i][j] > 4.0:
                    tinf.tinf.ups_comp[i][j] = 4.0
                elif tinf.tinf.ups_comp[i][j] < -4.0:
                    tinf.tinf.ups_comp[i][j] = -4.0

cdef void update_motif_counts(double mcnt[4][4][4096], double *zero, Sequence seq, _node* nod, int stage) nogil:
    cdef int     i
    cdef int     j
    cdef int     k
    cdef int     start
    cdef int     spacendx
    cdef _motif* mot      = &nod.mot

    if nod.type == node_type.STOP or nod.edge == 1:
        return
    if mot.len == 0:
        zero[0] += 1.0
        return

    if nod.strand == 1:
        start = nod.ndx
    else:
        start = seq.slen-1-nod.ndx

    # Stage 0:  Count all motifs.
    # If a motif is detected, it is counted for every distance in stage 0.
    # This is done to make sure off-distance good motifs are recognized.
    if stage == 0:
        for i in reversed(range(4)):
            for j in range(start - 18 - i, start - 5 - i):
                if j < 0:
                    continue
                if j <= start-16-i:
                    spacendx = 3
                elif j <= start-14-i:
                    spacendx = 2
                elif j >= start-7-i:
                    spacendx = 1
                else:
                    spacendx = 0
                for k in range(4):
                    mcnt[i][k][seq._mer_ndx(j, length=i+3, strand=nod.strand)] += 1.0
    # Stage 1:  Count only the best motif, but also count all its sub-motifs.
    elif stage == 1:
        mcnt[mot.len-3][mot.spacendx][mot.ndx] += 1.0;
        for i in range(mot.len - 3):
            for j in range(start-mot.spacer-mot.len, start-mot.spacer-i-2):
                if j < 0:
                    continue
                if j <= start-16-i:
                    spacendx = 3
                elif j <= start-14-i:
                    spacendx = 2
                elif j >= start-7-i:
                    spacendx = 1
                else:
                    spacendx = 0
                mcnt[i][spacendx][seq._mer_ndx(j, length=i+3, strand=nod.strand)] += 1.0
    # Stage 2:  Only count the highest scoring motif.
    elif stage == 2:
        mcnt[mot.len-3][mot.spacendx][mot.ndx] += 1.0

# --- Wrappers ---------------------------------------------------------------

cpdef inline void reset_node_scores(Nodes nodes) nogil:
    node.reset_node_scores(nodes.nodes, nodes.length)

cpdef inline void record_overlapping_starts(Nodes nodes, TrainingInfo tinf, bint is_meta = False) nogil:
    node.record_overlapping_starts(nodes.nodes, nodes.length, tinf.tinf, is_meta)

cpdef inline void eliminate_bad_genes(Nodes nodes, int ipath, TrainingInfo tinf) nogil:
    dprog.eliminate_bad_genes(nodes.nodes, ipath, tinf.tinf)

cpdef inline void tweak_final_starts(Genes genes, Nodes nodes, TrainingInfo tinf) nogil:
    gene.tweak_final_starts(genes.genes, genes.length, nodes.nodes, nodes.length, tinf.tinf)

cpdef inline void record_gene_data(Genes genes, Nodes nodes, TrainingInfo tinf, int sequence_index) nogil:
    gene.record_gene_data(genes.genes, genes.length, nodes.nodes, tinf.tinf, sequence_index)

cpdef inline void determine_sd_usage(TrainingInfo tinf) nogil:
    node.determine_sd_usage(tinf.tinf)

# --- Main functions ---------------------------------------------------------

cpdef TrainingInfo train(Sequence sequence, bint closed=False, bint force_nonsd=False, double start_weight=4.35, int translation_table=11):

    cdef int*         gc_frame
    cdef Nodes        nodes    = Nodes()
    cdef TrainingInfo tinf     = TrainingInfo(sequence.gc, start_weight, translation_table)
    cdef ConnectionScorer scorer = ConnectionScorer()

    with nogil:
        # find all the potential starts and stops
        add_nodes(nodes, sequence, tinf, closed=closed)
        nodes._sort()
        scorer._index(nodes)
        # scan all the ORFs looking for a potential GC bias in a particular
        # codon position, in order to acquire a good initial set of genes
        gc_frame = calc_most_gc_frame(sequence)
        if not gc_frame:
            raise MemoryError()
        node.record_gc_bias(gc_frame, nodes.nodes, nodes.length, tinf.tinf)
        free(gc_frame)
        # do an initial dynamic programming routine with just the GC frame bias
        # used as a scoring function.
        record_overlapping_starts(nodes, tinf, is_meta=False)
        ipath = dynamic_programming(nodes, tinf, scorer, final=False)
        # gather dicodon statistics for the training set
        calc_dicodon_gene(tinf, sequence, nodes, ipath)
        raw_coding_score(nodes, sequence, tinf)
        # determine if this organism uses Shine-Dalgarno and score the node
        rbs_score(nodes, sequence, tinf)
        train_starts_sd(nodes, sequence, tinf)
        if force_nonsd:
            tinf.tinf.uses_sd = False
        else:
            determine_sd_usage(tinf)
        if not tinf.tinf.uses_sd:
            train_starts_nonsd(nodes, sequence, tinf)

    return tinf

cpdef Predictions find_genes_single(Sequence sequence, TrainingInfo tinf, bint closed = False, int sequence_index = 1):

    cdef int   ipath
    cdef Genes genes = Genes()
    cdef Nodes nodes = Nodes()
    cdef ConnectionScorer scorer = ConnectionScorer()

    with nogil:
        # find all the potential starts and stops, and sort them
        add_nodes(nodes, sequence, tinf, closed=closed)
        nodes._sort()
        scorer._index(nodes)
        # second dynamic programming, using the dicodon statistics as the
        # scoring function
        reset_node_scores(nodes)
        score_nodes(nodes, sequence, tinf, closed=closed, is_meta=False)
        record_overlapping_starts(nodes, tinf, is_meta=True)
        ipath = dynamic_programming(nodes, tinf, scorer, final=True)
        # eliminate eventual bad genes in the nodes
        if nodes.length > 0:
            eliminate_bad_genes(nodes, ipath, tinf)
        # record genes
        add_genes(genes, nodes, ipath)
        tweak_final_starts(genes, nodes, tinf)
        record_gene_data(genes, nodes, tinf, sequence_index=sequence_index)

    # return the final predictions
    return Predictions(genes, nodes, sequence, tinf)

cpdef Predictions find_genes_meta(Sequence sequence, bint closed = False, int sequence_index = 1):

    cdef int          i
    cdef double       low
    cdef double       high
    cdef int          ipath
    cdef int          tt        = -1
    cdef Genes        genes     = Genes()
    cdef Nodes        nodes     = Nodes()
    cdef TrainingInfo tinf      = TrainingInfo.__new__(TrainingInfo)
    cdef ConnectionScorer scorer = ConnectionScorer()
    cdef int          max_phase = 0
    cdef double       max_score = -100.0

    with nogil:
        # make sure tinf does not deallocate by mistake
        tinf.owned = False

        # compute the min/max acceptable gc for the sequence to only
        # use appropriate metagenomic bins
        low = 0.88495*sequence.gc - 0.0102337
        high = 0.86596*sequence.gc + 0.1131991
        if low > 0.65:
            low = 0.65
        if high < 0.35:
            high = 0.35

        # check which of the metagenomeic bins gets the best results
        for i in range(NUM_META):
            # check which of the metagenomic bins gets the best results
            if _METAGENOMIC_BINS[i].tinf.gc < low or _METAGENOMIC_BINS[i].tinf.gc > high:
                continue
            # record the training information for the current bin
            tinf.tinf = _METAGENOMIC_BINS[i].tinf
            # recreate the node list if the translation table changed
            if tinf.tinf.trans_table != tt:
                tt = tinf.tinf.trans_table
                nodes._clear()
                add_nodes(nodes, sequence, tinf, closed=closed)
                nodes._sort()
                scorer._index(nodes)
            # compute the score for the current bin
            reset_node_scores(nodes)
            score_nodes(nodes, sequence, tinf, closed=closed, is_meta=True)
            record_overlapping_starts(nodes, tinf, is_meta=True)
            ipath = dynamic_programming(nodes, tinf, scorer, final=True)
            # update genes if the current bin had a better score
            if nodes.length > 0 and nodes.nodes[ipath].score > max_score:
                # record best phase and score
                max_phase = i
                max_score = nodes.nodes[ipath].score
                # eliminate eventual bad genes in the nodes
                eliminate_bad_genes(nodes, ipath, tinf)
                # clear the gene array
                genes._clear()
                # extract the genes from the dynamic programming array
                add_genes(genes, nodes, ipath)
                tweak_final_starts(genes, nodes, tinf)
                record_gene_data(genes, nodes, tinf, sequence_index=sequence_index)

        # recover the nodes corresponding to the best run
        tinf.tinf = _METAGENOMIC_BINS[max_phase].tinf
        nodes._clear()
        add_nodes(nodes, sequence, tinf, closed=closed)
        nodes._sort()
        score_nodes(nodes, sequence, tinf, closed=closed, is_meta=True)

    # return the final predictions
    return Predictions(genes, nodes, sequence, tinf)
