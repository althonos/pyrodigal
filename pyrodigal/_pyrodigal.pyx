# coding: utf-8
# cython: language_level=3, linetrace=True

"""Bindings to Prodigal, an ORF finder for genomes, progenomes and metagenomes.
"""

# ----------------------------------------------------------------------------

from cpython.exc cimport PyErr_CheckSignals
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from cpython.ref cimport Py_INCREF
from cpython.tuple cimport PyTuple_New, PyTuple_SET_ITEM
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

cdef class Sequence:
    """A compressed input sequence.
    """

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

    def __sizeof__(self):
        cdef size_t blen = self.slen//4 + (self.slen%4 != 0)
        cdef size_t ulen = self.slen//8 + (self.slen%8 != 0)
        return (ulen + 2 * blen) * sizeof(unsigned char) + sizeof(self)

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
            begin = gene.begin - 1
            end = gene.end - 1
            seq = self.owner.sequence.seq
            unk = gene.begin - 1
        else:
            begin = slen - gene.end
            end = slen - gene.begin
            seq = self.owner.sequence.rseq
            unk = gene.end - 3

        with nogil:
            for i, j in enumerate(range(begin, end, 3)):
                if bitmap.test(useq, unk) or bitmap.test(useq, unk+1) or bitmap.test(useq, unk+2):
                    aa = unknown_residue
                else:
                    aa = sequence.amino(seq, j, tinf, i==0 and edge==0)
                PyUnicode_WRITE(kind, data, i, aa)
                unk += 3 * strand

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

# --- C-level API ------------------------------------------------------------

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

cpdef inline void reset_node_scores(Nodes nodes) nogil:
    node.reset_node_scores(nodes.nodes, nodes.length)

cpdef inline void record_overlapping_starts(Nodes nodes, TrainingInfo tinf, bint is_meta = False) nogil:
    node.record_overlapping_starts(nodes.nodes, nodes.length, tinf.tinf, is_meta)

cpdef inline void score_nodes(Nodes nodes, Sequence seq, TrainingInfo tinf, bint closed=False, bint is_meta=False) nogil:
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

cpdef inline void raw_coding_score(Sequence sequence, Nodes nodes, TrainingInfo tinf) nogil:
    node.raw_coding_score(sequence.seq, sequence.rseq, sequence.slen, nodes.nodes, nodes.length, tinf.tinf)

cpdef inline void rbs_score(Sequence sequence, Nodes nodes, TrainingInfo tinf) nogil:
    node.rbs_score(sequence.seq, sequence.rseq, sequence.slen, nodes.nodes, nodes.length, tinf.tinf)

cpdef inline void train_starts_sd(Sequence sequence, Nodes nodes, TrainingInfo tinf) nogil:
    node.train_starts_sd(sequence.seq, sequence.rseq, sequence.slen, nodes.nodes, nodes.length, tinf.tinf)

cpdef inline void determine_sd_usage(TrainingInfo tinf) nogil:
    node.determine_sd_usage(tinf.tinf)

cpdef inline void train_starts_nonsd(Sequence sequence, Nodes nodes, TrainingInfo tinf) nogil:
    node.train_starts_nonsd(sequence.seq, sequence.rseq, sequence.slen, nodes.nodes, nodes.length, tinf.tinf)

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
        raw_coding_score(sequence, nodes, tinf)
        # determine if this organism uses Shine-Dalgarno and score the node
        rbs_score(sequence, nodes, tinf)
        train_starts_sd(sequence, nodes, tinf)
        if force_nonsd:
            tinf.tinf.uses_sd = False
        else:
            determine_sd_usage(tinf)
        if not tinf.tinf.uses_sd:
            train_starts_nonsd(sequence, nodes, tinf)

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
