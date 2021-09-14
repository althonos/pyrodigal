# coding: utf-8
# cython: language_level=3, linetrace=True

"""Bindings to Prodigal, an ORF finder for genomes, progenomes and metagenomes.
"""

# ----------------------------------------------------------------------------

from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AsString
from cpython.exc cimport PyErr_CheckSignals
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from cpython.ref cimport Py_INCREF
from cpython.tuple cimport PyTuple_New, PyTuple_SET_ITEM
from libc.math cimport sqrt, log, pow, fmax, fmin
from libc.stdint cimport uint8_t
from libc.stdio cimport printf
from libc.stdlib cimport free, qsort
from libc.string cimport memchr, memset, strstr

from pyrodigal.prodigal cimport bitmap, dprog, gene, node, sequence
from pyrodigal.prodigal.bitmap cimport bitmap_t
from pyrodigal.prodigal.gene cimport _gene
from pyrodigal.prodigal.metagenomic cimport NUM_META, _metagenomic_bin, initialize_metagenomic_bins
from pyrodigal.prodigal.node cimport _node, MIN_EDGE_GENE, MIN_GENE, cross_mask
from pyrodigal.prodigal.sequence cimport calc_most_gc_frame, _mask, node_type, rcom_seq
from pyrodigal.prodigal.training cimport _training
from pyrodigal._utils cimport _mini_training
from pyrodigal._unicode cimport *

# ----------------------------------------------------------------------------

import warnings
import threading

# --- Module-level constants -------------------------------------------------

cdef int    IDEAL_SINGLE_GENOME = 100000
cdef int    MIN_SINGLE_GENOME   = 20000
cdef size_t MIN_GENES_ALLOC     = 8
cdef size_t MIN_NODES_ALLOC     = 8 * MIN_GENES_ALLOC
cdef set   TRANSLATION_TABLES   = set(range(1, 7)) | set(range(9, 17)) | set(range(21, 26))

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


cdef class Sequence:
    """A compressed input sequence.
    """

    def __cinit__(self):
        self.slen = 0
        self.gc = 0.0
        self.seq = self.rseq = self.useq = self.digits = NULL

    def __dealloc__(self):
        PyMem_Free(self.seq)
        PyMem_Free(self.rseq)
        PyMem_Free(self.useq)
        PyMem_Free(self.digits)

    def __len__(self):
        return self.slen

    def __sizeof__(self):
        cdef size_t blen = self.slen//4 + (self.slen%4 != 0)
        cdef size_t ulen = self.slen//8 + (self.slen%8 != 0)
        return (ulen + 2 * blen) * sizeof(unsigned char) + self.slen * sizeof(uint8_t) + sizeof(self)

    cdef int _allocate(self, int slen) except 1:
        cdef size_t blen = slen//4 + (slen%4 != 0)
        cdef size_t ulen = slen//8 + (slen%8 != 0)

        self.slen = slen
        self.seq  = <bitmap_t> PyMem_Malloc(blen * sizeof(unsigned char))
        self.rseq = <bitmap_t> PyMem_Malloc(blen * sizeof(unsigned char))
        self.useq = <bitmap_t> PyMem_Malloc(ulen * sizeof(unsigned char))
        self.digits = <uint8_t*> PyMem_Malloc(slen * sizeof(uint8_t))

        if self.seq == NULL or self.rseq == NULL or self.useq == NULL or self.digits == NULL:
            raise MemoryError()

        with nogil:
            memset(self.seq, 0,  blen * sizeof(unsigned char))
            memset(self.rseq, 0, blen * sizeof(unsigned char))
            memset(self.useq, 0, ulen * sizeof(unsigned char))
            memset(self.digits, 0, slen * sizeof(uint8_t))

        return 0

    @classmethod
    def from_bytes(cls, const unsigned char[:] sequence):
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
                    pass
                elif letter == b'T' or letter == b't':
                    seq.digits[i] = T
                    bitmap.set(seq.seq, j)
                    bitmap.set(seq.seq, j+1)
                elif letter == b'G' or letter == b'g':
                    seq.digits[i] = G
                    bitmap.set(seq.seq, j)
                    gc_count += 1
                elif letter == b'C' or letter == b'c':
                    seq.digits[i] = C
                    bitmap.set(seq.seq, j+1)
                    gc_count += 1
                else:
                    seq.digits[i] = N
                    bitmap.set(seq.seq,  j+1)
                    bitmap.set(seq.useq, i)
            rcom_seq(seq.seq, seq.rseq, seq.useq, seq.slen)
            if seq.slen > 0:
                seq.gc = (<double> gc_count) / (<double> seq.slen)

        return seq

    @classmethod
    def from_string(cls, str sequence):
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
                    pass
                elif letter == u'T' or letter == u't':
                    seq.digits[i] = T
                    bitmap.set(seq.seq, j)
                    bitmap.set(seq.seq, j+1)
                elif letter == u'G' or letter == u'g':
                    seq.digits[i] = G
                    bitmap.set(seq.seq, j)
                    gc_count += 1
                elif letter == u'C' or letter == u'c':
                    seq.digits[i] = C
                    bitmap.set(seq.seq, j+1)
                    gc_count += 1
                else:
                    seq.digits[i] = N
                    bitmap.set(seq.seq,  j+1)
                    bitmap.set(seq.useq, i)
            rcom_seq(seq.seq, seq.rseq, seq.useq, seq.slen)
            if seq.slen > 0:
                seq.gc = (<double> gc_count) / (<double> seq.slen)

        return seq

    cdef inline bint _is_a(self, int i, int strand = 1) nogil:
        cdef uint8_t x

        if strand == 1:
            return self.digits[i] == A
        else:
            return self.digits[self.slen - 1 - i] == _translation[A]

    cdef inline bint _is_g(self, int i, int strand = 1) nogil:
        cdef uint8_t x

        if strand == 1:
            return self.digits[i] == G
        else:
            return self.digits[self.slen - 1 - i] == _translation[G]

    cdef inline bint _is_gc(self, int i, int strand = 1) nogil:
        cdef uint8_t x

        if strand == 1:
            x = self.digits[i]
        else:
            x = self.digits[self.slen - 1 - i]

        return x == C or x == G

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

# --- Nodes ------------------------------------------------------------------

cdef class Motif:
    pass

cdef class Node:

    @property
    def score(self):
        assert self.node != NULL
        return self.node.score

cdef class Nodes:

    def __cinit__(self):
        self.nodes = NULL
        self.capacity = 0
        self.length = 0

    def __init__(self):
        self._clear()

    def __dealloc__(self):
        PyMem_Free(self.nodes)

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
        with nogil:
            self._clear()

    cdef int _sort(self) nogil except 1:
        """Sort all nodes in the vector by their index and strand.
        """
        qsort(self.nodes, self.length, sizeof(_node), node.compare_nodes)

    def sort(self):
        with nogil:
            self._sort()

# --- Genes ------------------------------------------------------------------

cdef class Gene:
    """A single raw gene found by Prodigal within a DNA sequence.
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
        return self.gene.start_ndx

    @property
    def stop_ndx(self):
        return self.gene.stop_ndx

cdef class Genes:
    """A list of raw genes found by Prodigal in a single sequence.
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
        with nogil:
            self._clear()

# --- Training Info ----------------------------------------------------------

cdef class TrainingInfo:

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
        ``AGGA/GGAG/GAGG``, ``GGAG/GAGG``, ``AGGAG/GGAGG``, ``AGGAG``
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

        """

        cdef size_t   i
        cdef size_t   j
        cdef bitmap_t seq
        cdef size_t   begin
        cdef size_t   end
        cdef size_t   unk
        cdef size_t   length
        cdef Py_UCS4  nuc
        cdef _gene*   gene   = self.gene.gene
        cdef int      slen   = self.owner.sequence.slen
        cdef int      strand = self.owner.nodes.nodes[gene.start_ndx].strand
        cdef bitmap_t useq   = self.owner.sequence.useq

        # compute the right length to hold the nucleotides
        length = (<size_t> gene.end) - (<size_t> gene.begin) + 1

        # compute the offsets in the sequence bitmap
        if strand == 1:
            begin = gene.begin - 1
            end = gene.end
            seq = self.owner.sequence.seq
            unk = gene.begin - 1
        else:
            begin = slen - gene.end
            end = slen + 1 - gene.begin
            seq = self.owner.sequence.rseq
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
                if sequence.is_n(useq, unk):
                    nuc = u"N"
                elif sequence.is_a(seq, j):
                    nuc = u"A"
                elif sequence.is_t(seq, j):
                    nuc = u"T"
                elif sequence.is_g(seq, j):
                    nuc = u"G"
                else:
                    nuc = u"C"
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
            `ValueError`: when ``translation_table`` is not a valid number.

        """

        cdef size_t         nucl_length
        cdef size_t         prot_length
        cdef size_t         i
        cdef size_t         j
        cdef _mini_training mini_tinf
        cdef _training*     tinf
        cdef object         protein
        cdef int            kind
        cdef void*          data
        cdef Py_UCS4        aa
        cdef _gene*         gene        = self.gene.gene
        cdef int            slen        = self.owner.sequence.slen
        cdef bitmap_t       useq        = self.owner.sequence.useq
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
            tinf = self.owner.training_info.tinf
        else:
            if translation_table not in TRANSLATION_TABLES:
                raise ValueError(f"{translation_table} is not a valid translation table index")
            mini_tinf.trans_table = translation_table
            tinf = <_training*> &mini_tinf
            if tinf.trans_table != translation_table:
                raise RuntimeError("failed to dynamically change the translation table")

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
                aa = self.owner.sequence._amino(j, tinf.trans_table, strand=strand, is_init=i==0 and not edge)
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

    def __cinit__(self):
        self._num_seq = 1

    def __init__(self, meta=False, closed=False):
        """Instantiate and configure a new ORF finder.

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
                value has been manually selected by the PRODIGAL authors as an
                appropriate value for 99% of genomes.
            translation_table (`int`, optional): The translation table to use.

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
    for i in range(seq.slen-3, -1, -1):
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
    for i in range(seq.slen-3, -1, -1):
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

    cdef int  path = ipath
    cdef int  ng   = 0

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

cpdef int calc_orf_gc(Nodes nodes, Sequence seq, TrainingInfo tinf) nogil except -1:
    cdef int i
    cdef int j
    cdef int last[3]
    cdef int phase
    cdef double gc[3]
    cdef double gsize = 0.0

    # direct strand
    gc[0] = gc[1] = gc[2] = 0.0
    for i in range(nodes.length - 1, -1, -1):
        phase = nodes.nodes[i].ndx %3
        if nodes.nodes[i].strand == 1:
            if nodes.nodes[i].type == node_type.STOP:
                last[phase] = nodes.nodes[i].ndx
                gc[phase] = seq._is_gc(nodes.nodes[i].ndx) + seq._is_gc(nodes.nodes[i].ndx+1) + seq._is_gc(nodes.nodes[i].ndx+2)
            else:
                for j in range(last[phase] - 3, nodes.nodes[i].ndx - 1, -3):
                    gc[phase] = seq._is_gc(j) + seq._is_gc(j+1) + seq._is_gc(j+2)
                gsize = <float> nodes.nodes[i].stop_val - nodes.nodes[i].ndx + 3.0
                nodes.nodes[i].gc_cont = gc[phase] / gsize
                last[phase] = nodes.nodes[i].ndx

    # reverse strand
    gc[0] = gc[1] = gc[2] = 0.0
    for i in range(nodes.length):
        phase = nodes.nodes[i].ndx % 3
        if nodes.nodes[i].strand == -1:
            if nodes.nodes[i].type == node_type.STOP:
                last[phase] = nodes.nodes[i].ndx
                gc[phase] = seq._is_gc(nodes.nodes[i].ndx) + seq._is_gc(nodes.nodes[i].ndx-1) + seq._is_gc(nodes.nodes[i].ndx-2)
            else:
                for j in range(last[phase] + 3, nodes.nodes[i].ndx + 1, 3):
                    gc[phase] += seq._is_gc(j) + seq._is_gc(j+1) + seq._is_gc(j+2)
                gsize = <float> nodes.nodes[i].ndx - nodes.nodes[i].stop_val + 3.0
                nodes.nodes[i].gc_cont = gc[phase] / gsize
                last[phase] = nodes.nodes[i].ndx

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

    for i in range(3, -1, -1):
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
    for i in range(nn-1, -1, -1):
        phase = nodes.nodes[i].ndx%3
        if nodes.nodes[i].strand == 1 and nodes.nodes[i].type == node_type.STOP:
            last[phase] = nodes.nodes[i].ndx
            score[phase] = 0.0
        elif nodes.nodes[i].strand == 1:
            for j in range(last[phase]-3, nodes.nodes[i].ndx - 1, -3):
                score[phase] += tinf.tinf.gene_dc[seq._mer_ndx(j, length=6, strand=1)];
            nodes.nodes[i].cscore = score[phase]
            last[phase] = nodes.nodes[i].ndx
    score[0] = score[1] = score[2] = 0.0
    for i in range(nn):
        phase = nodes.nodes[i].ndx%3
        if nodes.nodes[i].strand == -1 and nodes.nodes[i].type == node_type.STOP:
            last[phase] = nodes.nodes[i].ndx
            score[phase] = 0.0
        elif nodes.nodes[i].strand == -1:
            for j in range(last[phase]+3, nodes.nodes[i].ndx+1, 3):
                score[phase] += tinf.tinf.gene_dc[seq._mer_ndx(seq.slen-1-j, length=6, strand=-1)]
            nodes.nodes[i].cscore = score[phase]
            last[phase] = nodes.nodes[i].ndx

    # Second Pass: Penalize start nodes with ascending coding to their left
    score[0] = score[1] = score[2] = -10000.0;
    for i in range(nn):
        phase = nodes.nodes[i].ndx%3
        if nodes.nodes[i].strand == 1 and nodes.nodes[i].type == node_type.STOP:
            score[phase] = -10000.0
        elif nodes.nodes[i].strand == 1:
          if nodes.nodes[i].cscore > score[phase]:
              score[phase] = nodes.nodes[i].cscore
          else:
              nodes.nodes[i].cscore -= score[phase] - nodes.nodes[i].cscore
    score[0] = score[1] = score[2] = -10000.0
    for i in range(nn-1, -1, -1):
        phase = nodes.nodes[i].ndx%3
        if nodes.nodes[i].strand == -1 and nodes.nodes[i].type == node_type.STOP:
            score[phase] = -10000.0
        elif nodes.nodes[i].strand == -1:
            if nodes.nodes[i].cscore > score[phase]:
                score[phase] = nodes.nodes[i].cscore
            else:
                nodes.nodes[i].cscore -= (score[phase] - nodes.nodes[i].cscore);

    # Third Pass: Add length-based factor to the score
    # Penalize start nodes based on length to their left
    for i in range(nn):
        phase = nodes.nodes[i].ndx%3
        if nodes.nodes[i].strand == 1 and nodes.nodes[i].type == node_type.STOP:
            score[phase] = -10000.0;
        elif nodes.nodes[i].strand == 1:
            gsize = ((<double> nodes.nodes[i].stop_val - nodes.nodes[i].ndx)+3.0)/3.0
            if gsize > 1000.0:
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
                nodes.nodes[i].cscore = 0.5*lfac;
            nodes.nodes[i].cscore += lfac;
    for i in range(nn-1, -1, -1):
        phase = nodes.nodes[i].ndx%3;
        if nodes.nodes[i].strand == -1 and nodes.nodes[i].type == node_type.STOP:
            score[phase] = -10000.0;
        elif nodes.nodes[i].strand == -1:
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

    for i in range(nodes.length):
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
            for j in range(seq.slen - nodes.nodes[i].ndx - 21, seq.slen - nodes.nodes[i].ndx - 6):
                if j > seq.slen-1:
                    continue
                cur_sc[0] = shine_dalgarno_exact(seq, j, seq.slen-1-nodes.nodes[i].ndx, tinf, strand=-1)
                cur_sc[1] = shine_dalgarno_mm(seq, j, seq.slen-1-nodes.nodes[i].ndx, tinf, strand=-1)
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
                for j in range(i-1, -1, -1):
                    if nodes.nodes[j].edge and nodes.nodes[i].stop_val == nodes.nodes[j].stop_val:
                        nodes.nodes[i].uscore += node.EDGE_UPS*tinf.tinf.st_wt
                        break
            elif i >= nodes.length - 500 and nodes.nodes[i].strand == -1:
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
        if start - i < 0:
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
    cdef double match[6]
    cdef double cur_ctr
    cdef double dis_flag

    limit = start - 4 - pos
    if limit > 6:
        limit = 6
    for i in range(limit, 6):
        match[i] = -10.0

    # Compare the 6-base region to AGGAGG
    for i in range(limit):
        if pos + i < 0:
            continue
        if i%3 == 0 and seq._is_a(pos+i, strand=strand):
            match[i] = 2.0
        elif i%3 != 0 and seq._is_g(pos+i, strand=strand):
            match[i] = 3.0
        else:
            match[i] = -10.0

    # Find the maximally scoring motif
    max_val = 0
    for i in range(limit, 2, -1):
        for j in range(limit+1-i):
            cur_ctr = -2.0;
            mism = 0;
            for k in range(j, j+i):
                cur_ctr += match[k]
                if match[k] < 0.0:
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
            if cur_ctr < 6.0:
                cur_val = 0
            elif cur_ctr == 6.0 and dis_flag == 2:
                cur_val = 1
            elif cur_ctr == 6.0 and dis_flag == 3:
                cur_val = 2
            elif cur_ctr == 8.0 and dis_flag == 3:
                cur_val = 3
            elif cur_ctr == 9.0 and dis_flag == 3:
                cur_val = 3
            elif cur_ctr == 6.0 and dis_flag == 1:
                cur_val = 6
            elif cur_ctr == 11.0 and dis_flag == 3:
                cur_val = 10
            elif cur_ctr == 12.0 and dis_flag == 3:
                cur_val = 10
            elif cur_ctr == 14.0 and dis_flag == 3:
                cur_val = 10
            elif cur_ctr == 8.0 and dis_flag == 2:
                cur_val = 11
            elif cur_ctr == 9.0 and dis_flag == 2:
                cur_val = 11
            elif cur_ctr == 8.0 and dis_flag == 1:
                cur_val = 12
            elif cur_ctr == 9.0 and dis_flag == 1:
                cur_val = 12
            elif cur_ctr == 6.0 and dis_flag == 0:
                cur_val = 13
            elif cur_ctr == 8.0 and dis_flag == 0:
                cur_val = 15
            elif cur_ctr == 9.0 and dis_flag == 0:
                cur_val = 16
            elif cur_ctr == 11.0 and dis_flag == 2:
                cur_val = 20
            elif cur_ctr == 11.0 and dis_flag == 1:
                cur_val = 21
            elif cur_ctr == 11.0 and dis_flag == 0:
                cur_val = 22
            elif cur_ctr == 12.0 and dis_flag == 2:
                cur_val = 20
            elif cur_ctr == 12.0 and dis_flag == 1:
                cur_val = 23
            elif cur_ctr == 12.0 and dis_flag == 0:
                cur_val = 24
            elif cur_ctr == 14.0 and dis_flag == 2:
                cur_val = 25
            elif cur_ctr == 14.0 and dis_flag == 1:
                cur_val = 26
            elif cur_ctr == 14.0 and dis_flag == 0:
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
    cdef double match[6]
    cdef double cur_ctr
    cdef double dis_flag

    limit = start - 4 - pos
    if limit > 6:
        limit = 6
    for i in range(limit, 6):
        match[i] = -10.0

    # Compare the 6-base region to AGGAGG
    for i in range(limit):
        if pos + i < 0:
            continue
        if i%3 == 0:
            match[i] = 2.0 if seq._is_a(pos+i, strand=strand) else -3.0
        else:
            match[i] = 3.0 if seq._is_g(pos+i, strand=strand) else -2.0

    # Find the maximally scoring motif
    max_val = 0
    for i in range(limit, 4, -1):
        for j in range(limit+1-i):
            cur_ctr = -2.0;
            mism = 0;
            for k in range(j, j+i):
                cur_ctr += match[k]
                if match[k] < 0.0:
                    mism += 1
                    if k <= j+1 or k >= j+i-2:
                      cur_ctr -= 10.0
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
            if rdis > 15 or cur_ctr < 6.0:
                continue

            # Single-Matching RBS Motifs
            if cur_ctr < 6.0:
                cur_val = 0
            elif cur_ctr == 6.0 and dis_flag == 3:
                cur_val = 2
            elif cur_ctr == 7.0 and dis_flag == 3:
                cur_val = 2
            elif cur_ctr == 9.0 and dis_flag == 3:
                cur_val = 3
            elif cur_ctr == 6.0 and dis_flag == 2:
                cur_val = 4
            elif cur_ctr == 6.0 and dis_flag == 1:
                cur_val = 5
            elif cur_ctr == 6.0 and dis_flag == 0:
                cur_val = 9
            elif cur_ctr == 7.0 and dis_flag == 2:
                cur_val = 7
            elif cur_ctr == 7.0 and dis_flag == 1:
                cur_val = 8
            elif cur_ctr == 7.0 and dis_flag == 0:
                cur_val = 14
            elif cur_ctr == 9.0 and dis_flag == 2:
                cur_val = 17
            elif cur_ctr == 9.0 and dis_flag == 1:
                cur_val = 18
            elif cur_ctr == 9.0 and dis_flag == 0:
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

    node.train_starts_sd(seq.seq, seq.rseq, seq.slen, nodes.nodes, nodes.length, tinf.tinf)

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
            if nodes.nodes[j].type != node_type.STOP and nodes.nodes[j].edge == 1:
                continue
            phase = nodes.nodes[j].ndx % 3
            if nodes.nodes[j].type == node_type.STOP and nodes.nodes[j].strand == 1:
                if best[phase] >= sthresh and nodes.nodes[bndx[phase]].ndx%3 == phase:
                    rreal[rbs[phase]] += 1.0;
                    treal[type[phase]] += 1.0;
                    if i == 9:
                        count_upstream_composition(seq, tinf, nodes.nodes[bndx[phase]].ndx, strand=1);
                best[phase] = 0.0; bndx[phase] = -1; rbs[phase] = 0; type[phase] = 0;
            elif nodes.nodes[j].strand == 1:
                if tinf.tinf.rbs_wt[nodes.nodes[j].rbs[0]] > tinf.tinf.rbs_wt[nodes.nodes[j].rbs[1]]+1.0 or nodes.nodes[j].rbs[1] == 0:
                    max_rb = nodes.nodes[j].rbs[0];
                elif tinf.tinf.rbs_wt[nodes.nodes[j].rbs[0]] < tinf.tinf.rbs_wt[nodes.nodes[j].rbs[1]]-1.0 or nodes.nodes[j].rbs[0] == 0:
                    max_rb = nodes.nodes[j].rbs[1];
                elif nodes.nodes[j].rbs[0] > nodes.nodes[j].rbs[1]:
                    max_rb = nodes.nodes[j].rbs[0]
                else:
                    max_rb = nodes.nodes[j].rbs[1]
                if nodes.nodes[j].cscore + wt*tinf.tinf.rbs_wt[max_rb] + wt*tinf.tinf.type_wt[nodes.nodes[j].type] >= best[phase]:
                    best[phase] = nodes.nodes[j].cscore + wt*tinf.tinf.rbs_wt[max_rb];
                    best[phase] += wt*tinf.tinf.type_wt[nodes.nodes[j].type];
                    bndx[phase] = j;
                    type[phase] = nodes.nodes[j].type;
                    rbs[phase] = max_rb;

        # Reverse strand pass
        for j in range(3):
            best[j] = 0.0
            bndx[j] = -1
            rbs[j] = 0
            type[j] = 0
        for j in range(nn-1, -1, -1):
            if nodes.nodes[j].type != node_type.STOP and nodes.nodes[j].edge:
                continue
            phase =  nodes.nodes[j].ndx % 3
            if nodes.nodes[j].type == node_type.STOP and nodes.nodes[j].strand == -1:
                if best[phase] >= sthresh and nodes.nodes[bndx[phase]].ndx%3 == phase:
                    rreal[rbs[phase]] += 1.0
                    treal[type[phase]] += 1.0
                    if i == 9:
                        count_upstream_composition(seq, tinf, nodes.nodes[bndx[phase]].ndx, strand=-1);
                best[phase] = 0.0
                bndx[phase] = -1
                rbs[phase] = 0
                type[phase] = 0
            elif nodes.nodes[j].strand == -1:
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
              tinf.tinf.ups_comp[i][j] /= sum;
              if tinf.tinf.gc > 0.1 and tinf.tinf.gc < 0.9:
                  if j == 0 or j == 3:
                      tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/(1.0-tinf.tinf.gc))
                  else:
                      tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/tinf.tinf.gc)
              elif tinf.tinf.gc <= 0.1:
                  if j == 0 or j == 3:
                      tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/0.90)
                  else:
                      tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/0.10)
              else:
                  if j == 0 or j == 3:
                      tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/0.10)
                  else:
                      tinf.tinf.ups_comp[i][j] = log(tinf.tinf.ups_comp[i][j]*2.0/0.90)
              if tinf.tinf.ups_comp[i][j] > 4.0:
                  tinf.tinf.ups_comp[i][j] = 4.0
              if tinf.tinf.ups_comp[i][j] < -4.0:
                  tinf.tinf.ups_comp[i][j] = -4.0


# --- Wrappers ---------------------------------------------------------------

cpdef inline void count_upstream_composition(Sequence seq, TrainingInfo tinf, int pos, int strand=1) nogil:
    if strand == 1:
        node.count_upstream_composition(seq.seq, seq.slen, strand, pos, tinf.tinf)
    else:
        node.count_upstream_composition(seq.rseq, seq.slen, strand, pos, tinf.tinf)

cpdef inline void reset_node_scores(Nodes nodes) nogil:
    node.reset_node_scores(nodes.nodes, nodes.length)

cpdef inline void record_overlapping_starts(Nodes nodes, TrainingInfo tinf, bint is_meta = False) nogil:
    node.record_overlapping_starts(nodes.nodes, nodes.length, tinf.tinf, is_meta)

cpdef inline int dynamic_programming(Nodes nodes, TrainingInfo tinf, bint is_meta = False) nogil:
    return dprog.dprog(nodes.nodes, nodes.length, tinf.tinf, is_meta)

cpdef inline void eliminate_bad_genes(Nodes nodes, int ipath, TrainingInfo tinf) nogil:
    dprog.eliminate_bad_genes(nodes.nodes, ipath, tinf.tinf)

cpdef inline void tweak_final_starts(Genes genes, Nodes nodes, TrainingInfo tinf) nogil:
    gene.tweak_final_starts(genes.genes, genes.length, nodes.nodes, nodes.length, tinf.tinf)

cpdef inline void record_gene_data(Genes genes, Nodes nodes, TrainingInfo tinf, int sequence_index) nogil:
    gene.record_gene_data(genes.genes, genes.length, nodes.nodes, tinf.tinf, sequence_index)

cpdef inline void calc_dicodon_gene(TrainingInfo tinf, Sequence sequence, Nodes nodes, int ipath) nogil:
    node.calc_dicodon_gene(tinf.tinf, sequence.seq, sequence.rseq, sequence.slen, nodes.nodes, ipath)

cpdef inline void train_starts_nonsd(Nodes nodes, Sequence sequence, TrainingInfo tinf) nogil:
    node.train_starts_nonsd(sequence.seq, sequence.rseq, sequence.slen, nodes.nodes, nodes.length, tinf.tinf)

cpdef inline void determine_sd_usage(TrainingInfo tinf) nogil:
    node.determine_sd_usage(tinf.tinf)


# --- Main functions ---------------------------------------------------------

cpdef TrainingInfo train(Sequence sequence, bint closed=False, bint force_nonsd=False, double start_weight=4.35, int translation_table=11):

    cdef int*         gc_frame
    cdef Nodes        nodes    = Nodes()
    cdef TrainingInfo tinf     = TrainingInfo(sequence.gc, start_weight, translation_table)

    with nogil:
        # find all the potential starts and stops
        add_nodes(nodes, sequence, tinf, closed=closed)
        nodes._sort()
        # scan all the ORFs looking for a potential GC bias in a particular
        # codon position, in order to acquire a good initial set of genes
        gc_frame = calc_most_gc_frame(sequence.seq, sequence.slen)
        if not gc_frame:
            raise MemoryError()
        node.record_gc_bias(gc_frame, nodes.nodes, nodes.length, tinf.tinf)
        free(gc_frame)
        # do an initial dynamic programming routine with just the GC frame bias
        # used as a scoring function.
        record_overlapping_starts(nodes, tinf, is_meta=False)
        ipath = dynamic_programming(nodes, tinf, is_meta=False)
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

    with nogil:
        # find all the potential starts and stops, and sort them
        add_nodes(nodes, sequence, tinf, closed=closed)
        nodes._sort()
        # second dynamic programming, using the dicodon statistics as the
        # scoring function
        reset_node_scores(nodes)
        score_nodes(nodes, sequence, tinf, closed=closed, is_meta=False)
        record_overlapping_starts(nodes, tinf, is_meta=True)
        ipath = dynamic_programming(nodes, tinf, is_meta=True)
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
            # compute the score for the current bin
            reset_node_scores(nodes)
            score_nodes(nodes, sequence, tinf, closed=closed, is_meta=True)
            record_overlapping_starts(nodes, tinf, is_meta=True)
            ipath = dynamic_programming(nodes, tinf, is_meta=True)
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
