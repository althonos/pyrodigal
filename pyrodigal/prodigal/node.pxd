from pyrodigal.prodigal.bitmap cimport bitmap_t
from pyrodigal.prodigal.sequence cimport _mask
from pyrodigal.prodigal.training cimport _training


cdef extern from "node.h":

    cdef const size_t STT_NOD = 100000
    cdef const size_t MIN_GENE = 90
    cdef const size_t MIN_EDGE_GENE = 60
    cdef const size_t MAX_SAM_OVLP = 60
    cdef const size_t ST_WINDOW = 60
    cdef const size_t OPER_DIST = 60
    cdef const double EDGE_BONUS = 0.74
    cdef const double EDGE_UPS = -1.0
    cdef const double META_PEN = 7.5

    cdef struct _motif:
        int ndx
        int len
        int spacer
        int spacendx
        double score

    cdef struct _node:
        int type
        int edge
        int ndx
        int strand
        int stop_val
        int star_ptr[3]
        int gc_bias
        double gc_score
        double cscore
        double gc_cont
        int rbs[2]
        _motif mot
        double uscore
        double tscore
        double rscore
        int traceb
        int tracef
        int ov_mark
        double score
        int elim

    int add_nodes(bitmap_t seq, bitmap_t rseq, int slen, _node* nodes, bint closed, _mask* mlist, int nm, _training* tinf) nogil
    cdef void reset_node_scores(_node*, int) nogil
    cdef int compare_nodes(const void*, const void*) nogil
    int stopcmp_nodes(const void*, const void*) nogil

    void record_overlapping_starts(_node*, int, _training*, int) nogil
    void record_gc_bias(int*, _node*, int, _training*) nogil

    void calc_dicodon_gene(_training*, unsigned char*, unsigned char*, int, _node*, int) nogil
    void calc_amino_bg(_training*, unsigned char*, unsigned char*, int, _node*, int) nogil

    void score_nodes(unsigned char*, unsigned char*, int, _node*, int, _training*, int, int) nogil
