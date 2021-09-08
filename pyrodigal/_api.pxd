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
    cdef          int      slen
    cdef          bitmap_t seq
    cdef          bitmap_t rseq
    cdef          bitmap_t useq
    cdef readonly double   gc

    cdef int _allocate(self, int slen) except 1

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

    cdef inline _node* _add_node(
        self,
        const int  ndx,
        const int  type,
        const int  strand,
        const int  stop_val,
        const bint edge,
    ) nogil except NULL

    cdef int _clear(self) nogil except 1
    cdef int _sort(self) nogil except 1

# ---

cdef class Gene:
    cdef Genes  owner
    cdef _gene* gene

cdef class Genes:
    cdef _gene* genes
    cdef size_t capacity
    cdef size_t length

    cdef inline _gene* _add_gene(
        self,
        const int begin,
        const int end,
        const int start_ndx,
        const int stop_ndx,
    ) nogil except NULL

    cdef int _clear(self) nogil except 1

# ---

cdef class TrainingInfo:
    cdef bint       owned
    cdef _training* tinf

cdef class MetagenomicBin:
    cdef          _metagenomic_bin* bin
    cdef readonly TrainingInfo      training_info

# ---

cdef class Prediction:
    cdef readonly Predictions owner
    cdef readonly Gene        gene

    cpdef unicode translate(
        self,
        object translation_table=?,
        Py_UCS4 unknown_residue=?,
    )

cdef class Predictions:
    cdef readonly Genes        genes
    cdef readonly Nodes        nodes
    cdef readonly Sequence     sequence
    cdef readonly TrainingInfo training_info

# ---

cdef _metagenomic_bin _METAGENOMIC_BINS[NUM_META]

# ---

cpdef int add_nodes(Nodes nodes, Sequence seq, TrainingInfo tinf, bint closed=*) nogil except -1
cpdef int add_genes(Genes genes, Nodes nodes, int ipath) nogil except -1

cpdef void reset_node_scores(Nodes nodes) nogil
cpdef void record_overlapping_starts(Nodes nodes, TrainingInfo tinf, bint is_meta=*) nogil
cpdef void score_nodes(Nodes nodes, Sequence seq, TrainingInfo tinf, bint closed=*, bint is_meta=*) nogil
cpdef int dynamic_programming(Nodes nodes, TrainingInfo tinf, bint is_meta=*) nogil
cpdef void eliminate_bad_genes(Nodes nodes, int ipath, TrainingInfo tinf) nogil
cpdef void tweak_final_starts(Genes genes, Nodes nodes, TrainingInfo tinf) nogil
cpdef void record_gene_data(Genes genes, Nodes nodes, TrainingInfo tinf, int sequence_index) nogil

cpdef Predictions find_genes_meta(Sequence seq, bint closed=*, int sequence_index = *)
