# coding: utf-8
# cython: language_level=3, linetrace=True

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

cdef size_t MIN_GENES_ALLOC   = 8
cdef size_t MIN_NODES_ALLOC   = 8 * MIN_GENES_ALLOC
cdef set   TRANSLATION_TABLES = set(range(1, 7)) | set(range(9, 17)) | set(range(21, 26))

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
        cdef int           gc_count = 0

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
    cdef Nodes  owner
    cdef _node* node

cdef class Nodes:
    cdef _node* nodes
    cdef size_t capacity
    cdef size_t length

    def __cinit__(self):
        self.nodes = NULL
        self.capacity = 0
        self.length = 0

    def __init__(self):
        self.clear()

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

    cdef _node* _add_node(
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
    cdef Genes  owner
    cdef _gene* gene


cdef class Genes:
    cdef _gene* genes
    cdef size_t capacity
    cdef size_t length

    def __cinit__(self):
        self.genes = NULL
        self.capacity = 0
        self.length = 0

    def __init__(self):
        self.clear()

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

    cdef _gene* _add_gene(
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

    cpdef void clear(self):
        """Remove all genes from the vector.
        """
        cdef size_t old_length
        old_length, self.length = self.length, 0
        with nogil:
            memset(self.genes, 0, old_length * sizeof(_gene))

# ---

cdef class TrainingInfo:
    cdef bint       owned
    cdef _training* tinf

    def __dealloc__(self):
        if self.owned:
            PyMem_Free(self.tinf)

cdef class MetagenomicBin:
    cdef          _metagenomic_bin* bin
    cdef readonly TrainingInfo      training_info

    @property
    def index(self):
        assert self.bin != NULL
        return self.bin.index

    @property
    def description(self):
        assert self.bin != NULL
        return self.bin.desc.decode('ascii')

# ---

cdef class Prediction:
    """A single predicted gene found by Prodigal.
    """
    cdef readonly Predictions owner
    cdef readonly Gene        gene
    # cdef readonly Nodes       nodes
    # cdef readonly Sequence     sequence
    # cdef readonly TrainingInfo training_info

    cpdef unicode translate(
        self,
        object translation_table=None,
        Py_UCS4 unknown_residue=u"X",
    ):
        """Translate the predicted gene into a protein sequence.

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
        cdef size_t         slen        = self.owner.sequence.slen
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
        nucl_length = (<size_t> gene.end) - (<size_t> gene.begin)
        prot_length = nucl_length//3 + (nucl_length%3 != 0)
        # create an empty protein string that we can write to
        # with the appropriate functions
        protein = PyUnicode_New(prot_length, 0x7F)
        kind    = PyUnicode_KIND(protein)
        data    = PyUnicode_DATA(protein)

        # compute the offsets in the sequence bitmaps:
        # - begin is the coordinates of the first nucleotide in the gene
        # - unk is the coordinate of the first nucleotide in the useq bitmap
        if strand == 1:
            begin = gene.begin
            end = gene.end
            seq = self.owner.sequence.seq
            unk = begin
        else:
            begin = slen + 1 - gene.end
            end = slen + 1 - gene.begin
            seq = self.owner.sequence.rseq
            unk = slen + 1 - begin

        # fill the sequence string, replacing residues with any unknown
        # nucleotide in the codon with an "X".
        with nogil:
            for i, j in enumerate(range(begin, end, 3)):
                if bitmap.test(useq, unk-1) or bitmap.test(useq, unk) or bitmap.test(useq, unk+1):
                    aa = "X"
                else:
                    aa = sequence.amino(seq, j-1, tinf, i==0 and edge==0)
                PyUnicode_WRITE(kind, data, i, aa)
                unk += 3 * strand

        # return the string containing the protein sequence
        return protein



cdef class Predictions:
    """A sequence of predictions found by Prodigal in a single sequence.
    """
    cdef readonly Genes        genes
    cdef readonly Nodes        nodes
    cdef readonly Sequence     sequence
    cdef readonly TrainingInfo training_info

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
    bin = MetagenomicBin.__new__(MetagenomicBin, )
    bin.bin = &_META_BINS[_i]
    bin.training_info = TrainingInfo.__new__(TrainingInfo)
    bin.training_info.owned = False
    bin.training_info.tinf = bin.bin.tinf
    METAGENOMIC_BINS.append(bin)
METAGENOMIC_BINS = tuple(METAGENOMIC_BINS)

# ---

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
    # TODO(@althonos): handle region masking
    cdef _mask* mlist = NULL
    cdef int    nm    = 0

    # If sequence is smaller than a codon, there are no nodes to add
    if seq.slen < 3:
        return nn

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
                    nodes._add_node(
                        ndx = last[i%3],
                        type = node_type.STOP,
                        strand = 1,
                        stop_val = i,
                        edge = not sequence.is_stop(seq.seq, last[i%3], tinf.tinf),
                    )
                    nn += 1
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
                        nodes._add_node(
                            ndx = i,
                            type = node_type.ATG,
                            stop_val = last[i%3],
                            strand = 1,
                            edge = False
                        )
                        nn += 1
                    elif sequence.is_ttg(seq.seq, i):
                        saw_start[i%3] = True
                        nodes._add_node(
                            ndx = i,
                            type = node_type.TTG,
                            stop_val = last[i%3],
                            strand = 1,
                            edge = False
                        )
                        nn += 1
                    elif sequence.is_gtg(seq.seq, i):
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
                    edge = not sequence.is_stop(seq.seq, last[i%3], tinf.tinf)
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
            if sequence.is_stop(seq.rseq, i, tinf.tinf):
                if saw_start[i%3]:
                    nodes._add_node(
                        ndx = seq.slen - last[i%3] - 1,
                        type = node_type.STOP,
                        strand = -1,
                        stop_val = seq.slen - i - 1,
                        edge = not sequence.is_stop(seq.rseq, last[i%3], tinf.tinf)
                    )
                    nn += 1
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
                        nodes._add_node(
                            ndx = seq.slen - i - 1,
                            type = node_type.ATG,
                            strand = -1,
                            stop_val = seq.slen - last[i%3] - 1,
                            edge = False
                        )
                        nn += 1
                    elif sequence.is_gtg(seq.rseq, i):
                        saw_start[i%3] = True
                        nodes._add_node(
                            ndx = seq.slen - i - 1,
                            type = node_type.GTG,
                            strand = -1,
                            stop_val = seq.slen - last[i%3] - 1,
                            edge = False
                        )
                        nn += 1
                    elif sequence.is_ttg(seq.rseq, i):
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
                    edge = not sequence.is_stop(seq.rseq, last[i%3], tinf.tinf),
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

    with nogil:
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

cpdef void reset_node_scores(Nodes nodes) nogil:
    with nogil:
        node.reset_node_scores(nodes.nodes, nodes.length)

cpdef void record_overlapping_starts(Nodes nodes, TrainingInfo tinf, bint is_meta = False) nogil:
    with nogil:
        node.record_overlapping_starts(nodes.nodes, nodes.length, tinf.tinf, is_meta)

cpdef void score_nodes(Nodes nodes, Sequence seq, TrainingInfo tinf, bint closed=False, bint is_meta=False) nogil:
    with nogil:
        node.score_nodes(
            seq.seq,
            seq.rseq,
            seq.slen,
            nodes.nodes,
            nodes.length,
            tinf.tinf,
            closed,
            is_meta
        )

cpdef int dynamic_programming(Nodes nodes, TrainingInfo tinf, bint is_meta = False) nogil:
    with nogil:
        return dprog.dprog(nodes.nodes, nodes.length, tinf.tinf, is_meta)

cpdef void eliminate_bad_genes(Nodes nodes, int ipath, TrainingInfo tinf) nogil:
    with nogil:
        dprog.eliminate_bad_genes(nodes.nodes, ipath, tinf.tinf)

cpdef void tweak_final_starts(Genes genes, Nodes nodes, TrainingInfo tinf) nogil:
    with nogil:
        gene.tweak_final_starts(genes.genes, genes.length, nodes.nodes, nodes.length, tinf.tinf)

cpdef void record_gene_data(Genes genes, Nodes nodes, TrainingInfo tinf, int sequence_index) nogil:
    with nogil:
        gene.record_gene_data(genes.genes, genes.length, nodes.nodes, tinf.tinf, sequence_index)

cpdef Predictions find_genes_meta(object sequence, bint closed = False):

    cdef Sequence seq
    cdef double   low
    cdef double   high
    cdef int      i
    cdef Genes    genes     = Genes()
    cdef Nodes    nodes     = Nodes()
    cdef int      max_phase = 0
    cdef double   max_score = -100.0

    if isinstance(sequence, str):
        seq = Sequence.from_string(sequence)
    else:
        seq = Sequence.from_bytes(sequence)

    # compute the min/max acceptable gc for the sequence to only
    # use appropriate metagenomic bins
    low = 0.88495*seq.gc - 0.0102337
    high = 0.86596*seq.gc + 0.1131991
    if low > 0.65:
        low = 0.65
    if high < 0.35:
        high = 0.35

    # check which of the metagenomeic bins gets the best results
    for i in range(NUM_META):
        # check which of the metagenomic bins gets the best results
        if _META_BINS[i].tinf.gc < low or _META_BINS[i].tinf.gc > high:
            continue
        # recreate the node list if the translation table changed
        if i == 0 or _META_BINS[i].tinf.trans_table != _META_BINS[i-1].tinf.trans_table:
            nodes.clear()
            add_nodes(nodes, seq, METAGENOMIC_BINS[i].training_info, closed)
            nodes.sort()
        # compute the score for the current bin
        reset_node_scores(nodes)
        score_nodes(nodes, seq, METAGENOMIC_BINS[i].training_info, closed, is_meta=True)
        record_overlapping_starts(nodes, METAGENOMIC_BINS[i].training_info, is_meta=True)
        ipath = dynamic_programming(nodes, METAGENOMIC_BINS[i].training_info, is_meta=True)
        # update genes if the current bin had a better score
        if nodes.length > 0 and nodes.nodes[ipath].score > max_score:
            # record best phase and score
            max_phase = i
            max_score = nodes.nodes[ipath].score
            # eliminate eventual bad genes in the nodes
            eliminate_bad_genes(nodes, ipath, METAGENOMIC_BINS[i].training_info)
            # clear the gene array
            genes.clear()
            # extract the genes from the dynamic programming array
            add_genes(genes, nodes, ipath)
            tweak_final_starts(genes, nodes, METAGENOMIC_BINS[i].training_info)
            record_gene_data(genes, nodes, METAGENOMIC_BINS[i].training_info, sequence_index=0)

    # recover the nodes corresponding to the best run
    nodes.clear()
    add_nodes(nodes, seq, METAGENOMIC_BINS[max_phase].training_info, closed=closed)
    nodes.sort()
    score_nodes(nodes, seq, METAGENOMIC_BINS[max_phase].training_info, closed=closed, is_meta=True)

    #
    return Predictions(genes, nodes, seq, METAGENOMIC_BINS[max_phase].training_info)
