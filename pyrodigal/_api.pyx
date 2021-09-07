from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.stdlib cimport realloc, calloc, malloc, free, qsort
from libc.string cimport memchr, memcmp, memcpy, memset, strcpy, strstr

from pyrodigal.prodigal cimport bitmap, dprog, gene, node, sequence
from pyrodigal.prodigal.bitmap cimport bitmap_t
from pyrodigal.prodigal.gene cimport MAX_GENES, _gene
from pyrodigal.prodigal.metagenomic cimport NUM_META, _metagenomic_bin, initialize_metagenomic_bins
from pyrodigal.prodigal.node cimport _node, MIN_EDGE_GENE, MIN_GENE, cross_mask
from pyrodigal.prodigal.sequence cimport calc_most_gc_frame, gc_content, _mask, node_type, rcom_seq
from pyrodigal.prodigal.training cimport _training
from pyrodigal._utils cimport _mini_training
from pyrodigal._unicode cimport *


# ---

DEF MIN_GENES_ALLOC = 8
DEF MIN_NODES_ALLOC = 8 * MIN_GENES_ALLOC

# ---

cdef class Sequence:
    """A compressed input sequence.
    """
    cdef int      slen
    cdef bitmap_t seq
    cdef bitmap_t rseq
    cdef bitmap_t useq
    cdef double   gc

    def __cinit__(self):
        self.slen = 0
        self.gc = 0.0
        self.seq = self.rseq = self.useq = NULL

    def __dealloc__(self):
        PyMem_Free(self.seq)
        PyMem_Free(self.rseq)
        PyMem_Free(self.useq)

    def __len__(self):
        return self.slen

    cdef int _allocate(self, int slen) except 1:
        cdef size_t blen = slen//4 + (slen%4 != 0)
        cdef size_t ulen = slen//8 + (slen%8 != 0)

        self.slen = slen
        self.seq  = <bitmap_t> PyMem_Malloc(blen * sizeof(unsigned char))
        self.rseq = <bitmap_t> PyMem_Malloc(blen * sizeof(unsigned char))
        self.useq = <bitmap_t> PyMem_Malloc(ulen * sizeof(unsigned char))

        if self.seq == NULL or self.rseq == NULL or self.useq == NULL:
            raise MemoryError()

        with nogil:
            memset(self.seq, 0,  blen * sizeof(unsigned char))
            memset(self.rseq, 0, blen * sizeof(unsigned char))
            memset(self.useq, 0, ulen * sizeof(unsigned char))

        return 0

    @classmethod
    def from_bytes(cls, const unsigned char[:] sequence):
        cdef int           i
        cdef int           j
        cdef unsigned char letter
        cdef Sequence      seq
        cdef int           gc_count

        seq = Sequence.__new__(Sequence)
        seq._allocate(sequence.shape[0])

        with nogil:
            for i, j in enumerate(range(0, seq.slen * 2, 2)):
                letter = sequence[i]
                if letter == b'A' or letter == b'a':
                    pass
                elif letter == b'T' or letter == b't':
                    bitmap.set(seq.seq, j)
                    bitmap.set(seq.seq, j+1)
                elif letter == b'G' or letter == b'g':
                    bitmap.set(seq.seq, j)
                    gc_count += 1
                elif letter == b'C' or letter == b'c':
                    bitmap.set(seq.seq, j+1)
                    gc_count += 1
                else:
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
        cdef int      gc_count

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
                    pass
                elif letter == u'T' or letter == u't':
                    bitmap.set(seq.seq, j)
                    bitmap.set(seq.seq, j+1)
                elif letter == u'G' or letter == u'g':
                    bitmap.set(seq.seq, j)
                    gc_count += 1
                elif letter == u'C' or letter == u'c':
                    bitmap.set(seq.seq, j+1)
                    gc_count += 1
                else:
                    bitmap.set(seq.seq,  j+1)
                    bitmap.set(seq.useq, i)
            rcom_seq(seq.seq, seq.rseq, seq.useq, seq.slen)
            if seq.slen > 0:
                seq.gc = (<double> gc_count) / (<double> seq.slen)

        return seq

# ---

cdef class Motif:
    pass

cdef class Node:
    pass

cdef class NodeVec:
    cdef _node* nodes
    cdef size_t capacity
    cdef size_t length

    def __cinit__(self):
        self.nodes = NULL
        self.capacity = 0
        self.length = 0

    def __dealloc__(self):
        PyMem_Free(self.nodes)

    def __len__(self):
        return self.length

    cdef _node* _add_node(
        self,
        int ndx,
        int type,
        int strand,
        int stop_val,
        bint edge,
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

    cpdef void add_nodes(self, Sequence seq, TrainingInfo tinf, bint closed=False):

        cdef int    i
        cdef int    last[3]
        cdef int    min_dist[3]
        cdef bint   saw_start[3]
        cdef int    slmod        = seq.slen % 3
        # TODO(@althonos): handle region masking
        cdef _mask* mlist = NULL
        cdef int    nm    = 0

        # If sequence is smaller than a codon, there are no nodes to add
        if seq.slen < 3:
            return

        with nogil:
            # Forward strand nodes
            for i in range(3):
                last[(i+slmod)%3] = seq.slen + i
                saw_start[i%3] = False
                min_dist[i%3] = MIN_EDGE_GENE
                if not closed:
                    while last[(i+slmod)%3] + 3 > seq.slen:
                        last[(i+slmod)%3] -= 3
            for i in range(seq.slen-3, -1, -1):
                if sequence.is_stop(seq.seq, i, tinf.tinf):
                    if saw_start[i%3]:
                        self._add_node(
                            ndx = last[i%3],
                            type = node_type.STOP,
                            strand = 1,
                            stop_val = i,
                            edge = not sequence.is_stop(seq.seq, last[i%3], tinf.tinf),
                        )
                    min_dist[i%3] = MIN_GENE
                    last[i%3] = i
                    saw_start[i%3] = False
                    continue
                if last[i%3] >= seq.slen:
                    continue
                if not cross_mask(i, last[i%3], mlist, nm):
                    if last[i%3] - i + 3 >= min_dist[i%3] and sequence.is_start(seq.seq, i, tinf.tinf):
                        if sequence.is_atg(seq.seq, i):
                            saw_start[i%3] = True
                            self._add_node(
                                ndx = i,
                                type = node_type.ATG,
                                stop_val = last[i%3],
                                strand = 1,
                                edge = False
                            )
                        elif sequence.is_ttg(seq.seq, i):
                            saw_start[i%3] = True
                            self._add_node(
                                ndx = i,
                                type = node_type.TTG,
                                stop_val = last[i%3],
                                strand = 1,
                                edge = False
                            )
                        elif sequence.is_gtg(seq.seq, i):
                            saw_start[i%3] = True
                            self._add_node(
                                ndx = i,
                                type = node_type.GTG,
                                stop_val = last[i%3],
                                strand = 1,
                                edge = False
                            )
                    if i <= 2 and not closed and last[i%3] - i > MIN_EDGE_GENE:
                        saw_start[i%3] = True
                        self._add_node(
                            ndx = i,
                            type = node_type.ATG,
                            stop_val = last[i%3],
                            strand = 1,
                            edge = True,
                        )
            for i in range(3):
                if saw_start[i%3]:
                    self._add_node(
                        ndx = last[i%3],
                        type = node_type.STOP,
                        strand = 1,
                        stop_val = i - 6,
                        edge = not sequence.is_stop(seq.seq, last[i%3], tinf.tinf)
                    )
            # Reverse strand nodes
            for i in range(3):
                last[(i + slmod) % 3] = seq.slen + i
                saw_start[i%3] = False
                min_dist[i%3] = MIN_EDGE_GENE
                if not closed:
                    while last[(i+slmod) % 3] + 3 > seq.slen:
                        last[(i+slmod)%3] -= 3
            for i in range(seq.slen-3, -1, -1):
                if sequence.is_stop(seq.rseq, i, tinf.tinf):
                    if saw_start[i%3]:
                        self._add_node(
                            ndx = seq.slen - last[i%3] - 1,
                            type = node_type.STOP,
                            strand = -1,
                            stop_val = seq.slen - i - 1,
                            edge = not sequence.is_stop(seq.rseq, last[i%3], tinf.tinf)
                        )
                    min_dist[i%3] = MIN_GENE
                    last[i%3] = i
                    saw_start[i%3] = False
                    continue
                if last[i%3] >= seq.slen:
                    continue
                if not cross_mask(i, last[i%3], mlist, nm):
                    if last[i%3] - i + 3 >= min_dist[i%3] and sequence.is_start(seq.rseq, i, tinf.tinf):
                        if sequence.is_atg(seq.rseq, i):
                            saw_start[i%3] = True
                            self._add_node(
                                ndx = seq.slen - i - 1,
                                type = node_type.ATG,
                                strand = -1,
                                stop_val = seq.slen - last[i%3] - 1,
                                edge = False
                            )
                        elif sequence.is_gtg(seq.rseq, i):
                            saw_start[i%3] = True
                            self._add_node(
                                ndx = seq.slen - i - 1,
                                type = node_type.GTG,
                                strand = -1,
                                stop_val = seq.slen - last[i%3] - 1,
                                edge = False
                            )
                        elif sequence.is_ttg(seq.rseq, i):
                            saw_start[i%3] = 1
                            self._add_node(
                                ndx = seq.slen - i - 1,
                                type = node_type.TTG,
                                strand = -1,
                                stop_val = seq.slen - last[i%3] - 1,
                                edge = False,
                            )
                if i <= 2 and not closed and last[i%3] - i > MIN_EDGE_GENE:
                    saw_start[i%3] = 1
                    node = self._add_node(
                        ndx = seq.slen - i - 1,
                        type = node_type.ATG,
                        strand = -1,
                        stop_val = seq.slen - last[i%3] - 1,
                        edge = True,
                    )
            for i in range(3):
                if saw_start[i%3]:
                    node = self._add_node(
                        ndx = seq.slen - last[i%3] - 1,
                        type = node_type.STOP,
                        strand = -1,
                        stop_val = seq.slen - i + 5,
                        edge = not sequence.is_stop(seq.rseq, last[i%3], tinf.tinf),
                    )

    cpdef void clear(self):
        """Remove all nodes from the vector.
        """
        cdef size_t old_length
        old_length, self.length = self.length, 0
        with nogil:
            memset(self.nodes, 0, old_length * sizeof(_node))

    cpdef void sort(self):
        """Sort all nodes in the vector by their index and strand.
        """
        assert self.nodes != NULL
        with nogil:
            qsort(self.nodes, self.length, sizeof(_node), node.compare_nodes)

# ---

cdef class Gene:
    pass

cdef class GeneVec:
    cdef _node* nodes
    cdef size_t capacity
    cdef size_t length

# ---


cdef class TrainingInfo:
    cdef _training* tinf


cdef class MetagenomicBin:
    cdef _metagenomic_bin* bin

    @property
    def index(self):
        assert self.bin != NULL
        return self.bin.index

    @property
    def description(self):
        assert self.bin != NULL
        return self.bin.desc.decode('ascii')

    @property
    def training_info(self):
        cdef TrainingInfo tinf
        tinf = TrainingInfo.__new__(TrainingInfo)
        tinf.tinf = self.bin.tinf
        return tinf


cdef _metagenomic_bin _META_BINS[NUM_META]
cdef ssize_t _i
for _i in range(NUM_META):
    memset(&_META_BINS[_i], 0, sizeof(_metagenomic_bin))
    _META_BINS[_i].tinf = <_training*> PyMem_Malloc(sizeof(_training))
    if not _META_BINS[_i].tinf:
        raise MemoryError()
    memset(_META_BINS[_i].tinf, 0, sizeof(_training))
initialize_metagenomic_bins(_META_BINS)


cdef MetagenomicBin bin
METAGENOMIC_BINS = []
for _i in range(NUM_META):
    bin = MetagenomicBin.__new__(MetagenomicBin)
    bin.bin = &_META_BINS[_i]
    METAGENOMIC_BINS.append(bin)
