# coding: utf-8
# cython: language_level=3, linetrace=True

# ----------------------------------------------------------------------------

from pyrodigal.prodigal.bitmap cimport bitmap_t
from pyrodigal.prodigal.gene cimport _gene
from pyrodigal.prodigal.metagenomic cimport NUM_META, _metagenomic_bin
from pyrodigal.prodigal.node cimport _node
from pyrodigal.prodigal.sequence cimport _mask
from pyrodigal.prodigal.training cimport _training

# --- Module-level constants -------------------------------------------------

cdef        int    MIN_SINGLE_GENOME
cdef        int    IDEAL_SINGLE_GENOME
cdef        size_t MIN_GENES_ALLOC
cdef        size_t MIN_NODES_ALLOC
cdef public set    TRANSLATION_TABLES

# --- Input sequence ---------------------------------------------------------

cdef class Sequence:
    cdef          int      slen
    cdef          bitmap_t seq
    cdef          bitmap_t rseq
    cdef          bitmap_t useq
    cdef readonly double   gc

    cdef int _allocate(self, int slen) except 1

# --- Nodes ------------------------------------------------------------------

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

# --- Genes ------------------------------------------------------------------

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

# --- Training Info ----------------------------------------------------------

cdef class TrainingInfo:
    cdef bint       owned
    cdef _training* tinf

# --- Metagenomic Bins -------------------------------------------------------

cdef class MetagenomicBin:
    cdef          _metagenomic_bin* bin
    cdef readonly TrainingInfo      training_info

cdef _metagenomic_bin _METAGENOMIC_BINS[NUM_META]

# --- Predictions ------------------------------------------------------------

cdef class Prediction:
    cdef readonly Predictions owner
    cdef readonly Gene        gene

    cpdef double confidence(self)
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

# --- Pyrodigal --------------------------------------------------------------

cdef class Pyrodigal:
    cdef readonly size_t       _num_seq
    cdef readonly bint         closed
    cdef readonly object       lock
    cdef readonly bint         meta
    cdef readonly TrainingInfo training_info

    cpdef Predictions  find_genes(self, object sequence)
    cpdef TrainingInfo train(self, object sequence, bint force_nonsd=*, double st_wt=*, int translation_table=*)

# --- C-level API ------------------------------------------------------------

cpdef int add_nodes(Nodes nodes, Sequence seq, TrainingInfo tinf, bint closed=*) nogil except -1
cpdef int add_genes(Genes genes, Nodes nodes, int ipath) nogil except -1

cpdef void reset_node_scores(Nodes nodes) nogil
cpdef void record_overlapping_starts(Nodes nodes, TrainingInfo tinf, bint is_meta=*) nogil
cpdef void score_nodes(Nodes nodes, Sequence seq, TrainingInfo tinf, bint closed=*, bint is_meta=*) nogil
cpdef int  dynamic_programming(Nodes nodes, TrainingInfo tinf, bint is_meta=*) nogil
cpdef void eliminate_bad_genes(Nodes nodes, int ipath, TrainingInfo tinf) nogil
cpdef void tweak_final_starts(Genes genes, Nodes nodes, TrainingInfo tinf) nogil
cpdef void record_gene_data(Genes genes, Nodes nodes, TrainingInfo tinf, int sequence_index) nogil
cpdef void calc_dicodon_gene(TrainingInfo tinf, Sequence sequence, Nodes nodes, int ipath) nogil
cpdef void raw_coding_score(Sequence sequence, Nodes nodes, TrainingInfo tinf) nogil
cpdef void rbs_score(Sequence sequence, Nodes nodes, TrainingInfo tinf) nogil
cpdef void train_starts_sd(Sequence sequence, Nodes nodes, TrainingInfo tinf) nogil
cpdef void determine_sd_usage(TrainingInfo tinf) nogil
cpdef void train_starts_nonsd(Sequence sequence, Nodes nodes, TrainingInfo tinf) nogil

# --- Main functions ---------------------------------------------------------

cpdef TrainingInfo train(Sequence sequence, bint closed=*, bint force_nonsd=*, double start_weight=*, int translation_table=*)
cpdef Predictions find_genes_single(Sequence sequence, TrainingInfo tinf, bint closed=*, int sequence_index=*)
cpdef Predictions find_genes_meta(Sequence seq, bint closed=*, int sequence_index=*)
