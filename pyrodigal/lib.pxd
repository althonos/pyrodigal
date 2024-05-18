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

cdef        size_t MIN_MASKS_ALLOC
cdef        size_t MIN_GENES_ALLOC
cdef        size_t MIN_NODES_ALLOC

cdef public int    _MIN_SINGLE_GENOME
cdef public int    _IDEAL_SINGLE_GENOME
cdef public set    _TRANSLATION_TABLES
cdef public str    _PRODIGAL_VERSION

# --- Sequence mask ----------------------------------------------------------

cdef class Mask:
    cdef _mask  _data
    cdef Masks  owner
    cdef _mask* mask

    cpdef size_t __sizeof__(self)

    @staticmethod
    cdef bint _intersects(_mask* mask, int begin, int end) noexcept nogil

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
    ) except NULL nogil
    cdef int _clear(self) noexcept nogil

    cpdef size_t __sizeof__(self)

    cpdef Masks copy(self)
    cpdef void clear(self)


# --- Input sequence ---------------------------------------------------------

cdef class Sequence:
    cdef          Py_ssize_t slen
    cdef          uint8_t*   digits
    cdef readonly double     gc
    cdef readonly double     gc_known
    cdef readonly size_t     unknown
    cdef readonly Masks      masks

    cdef int _build(
        self,
        const int      kind,
        const void*    data,
        const size_t   length,
    ) except 1 nogil
    cdef int _mask(self, const size_t mask_size) except 1 nogil
    cdef int _allocate(self, int slen) except 1
    cdef int* _max_gc_frame_plot(self, int window_size) except NULL nogil
    cdef char _amino(
        self,
        int i,
        int tt,
        int strand=*,
        bint is_init=*,
        char unknown_residue=*,
        bint strict=*,
    ) noexcept nogil
    cdef int _shine_dalgarno_exact(
        self,
        const int pos,
        const int start,
        const double rbs_wt[28],
        const int strand
    ) noexcept nogil
    cdef int _shine_dalgarno_mm(
        self,
        const int pos,
        const int start,
        const double rbs_wt[28],
        const int strand
    ) noexcept nogil

    cpdef size_t __sizeof__(self)

    cpdef object max_gc_frame_plot(self, int window_size=*)
    cpdef int shine_dalgarno(
        self,
        int pos,
        int start,
        TrainingInfo training_info,
        int strand=*,
        bint exact=*
    ) except -1

    cpdef double start_probability(self) noexcept
    cpdef double stop_probability(self) noexcept
        
        

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

    cdef int _index(self, Nodes nodes) except -1 nogil
    cdef int _compute_skippable(
        self,
        const int min,
        const int i
    ) noexcept nogil
    
    cdef void _score_node_connections(
        self,
        Nodes nodes,
        const int min,
        const int i,
        const _training* tinf,
        const bint final
    ) noexcept nogil
    cdef void _score_connections(
        self,
        Nodes nodes,
        const _training* tinf,
        const bint final
    ) noexcept nogil
    cdef int _find_max_index(
        self, 
        Nodes nodes
    ) noexcept nogil
    cdef void _disentangle_overlaps(
        self, 
        Nodes nodes, 
        int max_index
    ) noexcept nogil
    cdef void _max_forward_pointers(
        self, 
        Nodes nodes, 
        int max_index
    ) noexcept nogil
    cdef int _dynamic_programming(
        self,
        Nodes nodes,
        const _training* tinf,
        const bint final
    ) noexcept nogil

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
    ) noexcept nogil
    @staticmethod
    cdef void _score_upstream_composition(
        _node* node,
        Sequence seq,
        const _training* tinf,
    ) noexcept nogil

cdef class Nodes:
    # contiguous array of nodes, with capacity and length
    cdef          void*  nodes_raw
    cdef          _node* nodes
    cdef readonly size_t capacity
    cdef readonly size_t length

    cpdef size_t __sizeof__(self)

    cdef int _allocate(self, size_t capacity) except 1
    cdef inline _node* _add_node(
        self,
        const int  ndx,
        const int  type,
        const int  strand,
        const int  stop_val,
        const bint edge,
    ) except NULL nogil
    cdef int _calc_orf_gc(self, Sequence seq) except -1 nogil
    cdef int _clear(self) except 1 nogil
    cdef int _extract(
        self,
        Sequence sequence,
        const int translation_table,
        const bint closed,
        const int min_gene,
        const int min_edge_gene
    ) except -1 nogil
    cdef int _raw_coding_score(
        self,
        Sequence seq,
        const _training* tinf
    ) except -1 nogil
    cdef int _rbs_score(
        self,
        Sequence seq,
        const _training* tinf
    ) except -1 nogil
    cdef void _record_overlapping_starts(
        self,
        const _training* tinf,
        const int flag,
        const int max_sam_overlap,
    ) noexcept nogil
    cdef int _reset_scores(self) except 1 nogil
    cdef int _score(
        self,
        Sequence seq,
        const _training* tinf,
        const bint closed,
        const bint is_meta
    ) except -1 nogil
    cdef int _sort(self) except 1 nogil

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

    cpdef str _gene_data(self, object sequence_id)
    cpdef str _score_data(self)

    cpdef double confidence(self)
    cpdef str sequence(self)
    cpdef str translate(
        self,
        object translation_table=?,
        char unknown_residue=*,
        bint include_stop=*,
        bint strict=*,
    )

cdef class Genes:
    # Raw gene array
    cdef          _gene*         genes
    cdef          size_t         capacity
    cdef          size_t         length
    # Dynamic programming path
    cdef          int            ipath
    # References to source data
    cdef          size_t         _num_seq
    cdef readonly bint           meta
    cdef readonly Nodes          nodes
    cdef readonly Sequence       sequence
    cdef readonly TrainingInfo   training_info
    cdef readonly MetagenomicBin metagenomic_bin

    cpdef size_t __sizeof__(self)

    cdef int _allocate(self, size_t capacity) except 1
    cdef inline _gene* _add_gene(
        self,
        const int begin,
        const int end,
        const int start_ndx,
        const int stop_ndx,
    ) except NULL nogil
    cdef int _extract(self, Nodes nodes, int ipath) except -1 nogil
    cdef void _tweak_final_starts(
        self,
        Nodes nodes,
        const _training* tinf,
        const int max_sam_overlap,
    ) noexcept nogil
    cdef int _clear(self) except 1 nogil

    cpdef ssize_t write_genbank(self, object file, str sequence_id, str division=*, object date=*, object translation_table=?, bint strict_translation=*) except -1
    cpdef ssize_t write_gff(self, object file, str sequence_id, bint header=*, bint include_translation_table=*, bint full_id=*) except -1
    cpdef ssize_t write_genes(self, object file, str sequence_id, object width=*, bint full_id=*) except -1
    cpdef ssize_t write_translations(self, object file, str sequence_id, object width=*, object translation_table=?, bint include_stop=*, bint strict_translation=*, bint full_id=*) except -1
    cpdef ssize_t write_scores(self, object file, str sequence_id, bint header=*) except -1


# --- Training Info ----------------------------------------------------------

cdef class TrainingInfo:
    cdef bint       owned
    cdef _training* tinf

    cpdef size_t __sizeof__(self)

    @staticmethod
    cdef void _update_motif_counts(
        double mcnt[4][4][4096],
        const double *zero,
        Sequence seq,
        _node* nod,
        int stage
    ) noexcept nogil

    cdef void _calc_dicodon_gene(self, Sequence seq, _node* nodes, int ipath) noexcept nogil
    cdef void _count_upstream_composition(self, Sequence seq, int pos, int strand=*) noexcept nogil
    cdef void _train_starts_nonsd(self, Nodes nodes, Sequence seq) noexcept nogil
    cdef void _train_starts_sd(self, Nodes nodes, Sequence seq) noexcept nogil

    cpdef dict to_dict(self)
    cpdef object dump(self, object fp)


# --- Metagenomic Bins -------------------------------------------------------

cdef class MetagenomicBins:
    cdef readonly tuple              _objects
    cdef          _metagenomic_bin** bins
    cdef          size_t length

    @staticmethod
    cdef MetagenomicBins from_array(_metagenomic_bin* bins, size_t length) except *
    @staticmethod
    cdef MetagenomicBins from_initializer(
        void (*initializer)(_metagenomic_bin*),
        size_t length,
    ) except *

cdef class MetagenomicBin:
    cdef          _metagenomic_bin* bin
    cdef readonly TrainingInfo      training_info

cdef _metagenomic_bin _METAGENOMIC_BINS[NUM_META]


# --- GeneFinder --------------------------------------------------------------

cdef class GeneFinder:
    cdef readonly size_t          _num_seq
    cdef readonly str             backend
    cdef readonly bint            closed
    cdef readonly object          lock
    cdef readonly bint            mask
    cdef readonly int             min_mask
    cdef readonly int             max_overlap
    cdef readonly bint            meta
    cdef readonly MetagenomicBins metagenomic_bins
    cdef readonly int             min_gene
    cdef readonly int             min_edge_gene
    cdef readonly TrainingInfo    training_info

    cdef int _train(
        self,
        Sequence sequence,
        Nodes nodes,
        ConnectionScorer scorer,
        TrainingInfo tinf,
        bint force_nonsd,
    ) except -1 nogil
    cdef int _find_genes_single(
        self,
        Sequence sequence,
        TrainingInfo tinf,
        ConnectionScorer scorer,
        Nodes nodes,
        Genes genes,
    ) except -1 nogil
    cdef ssize_t _find_genes_meta(
        self,
        Sequence sequence,
        ConnectionScorer scorer,
        Nodes nodes,
        Genes genes,
    ) except? -1 nogil

    cpdef Genes find_genes(self, object sequence)
