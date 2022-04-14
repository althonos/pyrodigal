from libc.stdint cimport uint8_t
from pyrodigal.prodigal.node cimport _node
from pyrodigal.prodigal.training cimport _training


cdef extern from "_connection.h" nogil:

    cdef double _intergenic_mod_diff(const _node* n1, const _node* n2, const double start_weight)
    cdef double _intergenic_mod_same(const _node* n1, const _node* n2, const double start_weight)
    cdef double _intergenic_mod(const _node* n1, const _node* n2, const double start_weight)

    cdef void _score_connection_forward_start(
        const _node*     nodes,
        const _node*     n1,
              _node*     n2,
        const _training* tinf,
        const int        final,
    )
    cdef void _score_connection_forward_stop(
        const _node*     nodes,
        const _node*     n1,
              _node*     n2,
        const _training* tinf,
        const int        final,
    )
    cdef void _score_connection_backward_start(
        const _node*     nodes,
        const _node*     n1,
              _node*     n2,
        const _training* tinf,
        const int        final,
    )
    cdef void _score_connection_backward_stop(
        const _node*     nodes,
        const _node*     n1,
              _node*     n2,
        const _training* tinf,
        const int        final,
    )

    cdef void _score_connection(
        const _node*     nodes,
        const _node*     n1,
              _node*     n2,
        const _training* tinf,
        const int        final,
    )
