# coding: utf-8
# cython: language_level=3, linetrace=True

# ----------------------------------------------------------------------------

from libc.stdint cimport int8_t, uint8_t

from pyrodigal.prodigal.bitmap cimport bitmap_t
from pyrodigal.prodigal.metagenomic cimport NUM_META, _metagenomic_bin
from pyrodigal.prodigal.node cimport _node, _motif
from pyrodigal.prodigal.sequence cimport _mask
from pyrodigal.prodigal.training cimport _training

# --- Module-level constants -------------------------------------------------

cdef        int    MIN_SINGLE_GENOME
cdef        int    IDEAL_SINGLE_GENOME
cdef        size_t MIN_GENES_ALLOC
cdef        size_t MIN_NODES_ALLOC
cdef public set    _TRANSLATION_TABLES


# --- Sequence mask ----------------------------------------------------------

cdef class Mask:
    cdef _mask  _data
    cdef Masks  owner
    cdef _mask* mask

    @staticmethod
    cdef bint _intersects(_mask* mask, int begin, int end) nogil

    cpdef bint intersects(self, int begin, int end)

cdef class Masks:
    cdef _mask* masks
    cdef size_t capacity
    cdef size_t length

    cdef int _allocate(self, size_t capacity) except 1
    cdef inline _mask* _add_mask(
        self,
        const int  begin,
        const int  end,
    ) nogil except NULL
    cdef int _clear(self) nogil

    cpdef size_t __sizeof__(self)
    cpdef list __getstate__(self)
    cpdef object __setstate__(self, list state)

    cpdef Masks copy(self)
    cpdef void clear(self)


# --- Input sequence ---------------------------------------------------------

cdef class Sequence:
    cdef          Py_ssize_t slen
    cdef          uint8_t*   digits
    cdef readonly double     gc
    cdef readonly Masks      masks

    @staticmethod
    cdef int _build(
        const int      kind,
        const void*    data,
        const size_t   length,
              double*  gc,
              uint8_t* digits,
              Masks    masks,
    ) nogil except 1

    cdef int _allocate(self, int slen) except 1
    cdef char _amino(
        self,
        int i,
        int tt,
        int strand=*,
        bint is_init=*,
        char unknown_residue=*
    ) nogil
    cdef int _shine_dalgarno_exact(
        self,
        const int pos,
        const int start,
        const _training* tinf,
        const int strand
    ) nogil
    cdef int _shine_dalgarno_mm(
        self,
        const int pos,
        const int start,
        const _training* tinf,
        const int strand
    ) nogil

    cpdef size_t __sizeof__(self)
    cpdef dict __getstate__(self)
    cpdef object __setstate__(self, dict state)

    cpdef int shine_dalgarno(
        self,
        int pos,
        int start,
        TrainingInfo training_info,
        int strand=*,
        bint exact=*
    ) except -1


# --- Connection Scorer ------------------------------------------------------

cdef class ConnectionScorer:
    # which SIMD backend to use
    cdef uint8_t  backend
    # capacity of bypassing buffers
    cdef size_t   capacity
    # connection skip lookup table
    cdef uint8_t* skip_connection
    cdef uint8_t* skip_connection_raw
    # aligned storage of node types
    cdef uint8_t* node_types
    cdef uint8_t* node_types_raw
    # aligned storage of node strands
    cdef int8_t*  node_strands
    cdef int8_t*  node_strands_raw
    # aligned storage of node frame
    cdef uint8_t* node_frames
    cdef uint8_t* node_frames_raw

    cpdef size_t __sizeof__(self)

    cdef int _index(self, Nodes nodes) nogil except -1
    cdef int _compute_skippable(
        self,
        const int min,
        const int i
    ) nogil
    cdef void _score_connections(
        self,
        Nodes nodes,
        const int min,
        const int i,
        const _training* tinf,
        const bint final
    ) nogil


# --- Nodes ------------------------------------------------------------------

cdef class Motif:
    cdef Node    owner
    cdef _motif* motif

cdef class Node:
    cdef Nodes  owner
    cdef _node* node

    @staticmethod
    cdef void _find_best_upstream_motif(
        _node* node,
        Sequence seq,
        const _training* tinf,
        const int stage
    ) nogil
    @staticmethod
    cdef void _score_upstream_composition(
        _node* node,
        Sequence seq,
        const _training* tinf,
    ) nogil

cdef class Nodes:
    # contiguous array of nodes, with capacity and length
    cdef          _node* nodes
    cdef readonly size_t capacity
    cdef readonly size_t length

    cpdef size_t __sizeof__(self)
    cpdef list __getstate__(self)
    cpdef object __setstate__(self, list state)

    cdef int _allocate(self, size_t capacity) except 1
    cdef inline _node* _add_node(
        self,
        const int  ndx,
        const int  type,
        const int  strand,
        const int  stop_val,
        const bint edge,
    ) nogil except NULL
    cdef int _calc_orf_gc(self, Sequence seq) nogil except -1
    cdef int _clear(self) nogil except 1
    cdef int _dynamic_programming(
        self,
        const _training* tinf,
        ConnectionScorer scorer,
        const bint final
    ) nogil
    cdef int _extract(
        self,
        Sequence sequence,
        const int translation_table,
        const bint closed,
        const int min_gene,
        const int min_edge_gene
    ) nogil except -1
    cdef int _raw_coding_score(
        self,
        Sequence seq,
        const _training* tinf
    ) nogil except -1
    cdef int _rbs_score(
        self,
        Sequence seq,
        const _training* tinf
    ) nogil except -1
    cdef void _record_overlapping_starts(
        self,
        const _training* tinf,
        const int flag,
        const int max_sam_overlap,
    ) nogil
    cdef int _reset_scores(self) nogil except 1
    cdef int _score(
        self,
        Sequence seq,
        const _training* tinf,
        const bint closed,
        const bint is_meta
    ) nogil except -1
    cdef int _sort(self) nogil except 1

    cpdef Nodes copy(self)


# --- Genes ------------------------------------------------------------------

cdef struct _gene:
    int begin
    int end
    int start_ndx
    int stop_ndx

cdef class Gene:
    cdef Genes  owner
    cdef _gene* gene

    cpdef double confidence(self)
    cpdef str sequence(self)
    cpdef str translate(
        self,
        object translation_table=?,
        char unknown_residue=*,
    )

cdef class Genes:
    # Raw gene array
    cdef          _gene*        genes
    cdef          size_t       capacity
    cdef          size_t       length
    # References to source data
    cdef          size_t       _num_seq
    cdef readonly Nodes        nodes
    cdef readonly Sequence     sequence
    cdef readonly TrainingInfo training_info

    cpdef size_t __sizeof__(self)
    cpdef dict __getstate__(self)
    cpdef object __setstate__(self, dict state)

    cdef int _allocate(self, size_t capacity) except 1
    cdef inline _gene* _add_gene(
        self,
        const int begin,
        const int end,
        const int start_ndx,
        const int stop_ndx,
    ) nogil except NULL
    cdef int _extract(self, Nodes nodes, int ipath) nogil except -1
    cdef void _tweak_final_starts(
        self,
        Nodes nodes,
        const _training* tinf,
        const int max_sam_overlap,
    ) nogil
    cdef int _clear(self) nogil except 1

    cpdef ssize_t write_gff(self, object file, str prefix=*) except -1
    cpdef ssize_t write_genes(self, object file, str prefix=*, object width=*) except -1
    cpdef ssize_t write_translations(self, object file, str prefix=*, object width=*, object translation_table=?) except -1
    cpdef ssize_t write_scores(self, object file, bint header=*) except -1


# --- Training Info ----------------------------------------------------------

cdef class TrainingInfo:
    cdef bint       owned
    cdef _training* tinf

    cpdef size_t __sizeof__(self)
    cpdef dict __getstate__(self)
    cpdef object __setstate__(self, dict state)

    @staticmethod
    cdef void _update_motif_counts(double mcnt[4][4][4096], double *zero, Sequence seq, _node* nod, int stage) nogil

    cdef void _calc_dicodon_gene(self, Sequence seq, _node* nodes, int ipath) nogil
    cdef void _count_upstream_composition(self, Sequence seq, int pos, int strand=*) nogil
    cdef void _train_starts_nonsd(self, Nodes nodes, Sequence seq) nogil
    cdef void _train_starts_sd(self, Nodes nodes, Sequence seq) nogil

    cpdef object dump(self, object fp)


# --- Metagenomic Bins -------------------------------------------------------

cdef class MetagenomicBin:
    cdef          _metagenomic_bin* bin
    cdef readonly TrainingInfo      training_info

cdef _metagenomic_bin _METAGENOMIC_BINS[NUM_META]


# --- OrfFinder --------------------------------------------------------------

cdef class OrfFinder:
    cdef readonly size_t       _num_seq
    cdef readonly bint         closed
    cdef readonly object       lock
    cdef readonly bint         meta
    cdef readonly bint         mask
    cdef readonly int          min_gene
    cdef readonly int          min_edge_gene
    cdef readonly int          max_overlap
    cdef readonly TrainingInfo training_info

    cpdef dict __getstate__(self)
    cpdef object __setstate__(self, dict state)

    cdef int _train(
        self,
        Sequence sequence,
        Nodes nodes,
        ConnectionScorer scorer,
        TrainingInfo tinf,
        bint force_nonsd,
    ) nogil except -1
    cdef int _find_genes_single(
        self,
        Sequence sequence,
        TrainingInfo tinf,
        ConnectionScorer scorer,
        Nodes nodes,
        Genes genes,
    ) nogil except -1
    cdef int _find_genes_meta(
        self,
        Sequence sequence,
        ConnectionScorer scorer,
        Nodes nodes,
        Genes genes,
    ) nogil except -1

    cpdef Genes find_genes(self, object sequence)


# --- C-level API reimplementation -------------------------------------------

cdef int* calc_most_gc_frame(Sequence seq) nogil except NULL
