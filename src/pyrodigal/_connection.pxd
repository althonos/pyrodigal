from libc.stdint cimport uint8_t, int8_t
from pyrodigal.prodigal.node cimport _node
from pyrodigal.prodigal.training cimport _training


cdef extern from "_connection.h" nogil:


    cdef double _intergenic_mod_diff(const _node* n1, const _node* n2, const double start_weight)
    cdef double _intergenic_mod_same(const _node* n1, const _node* n2, const double start_weight)
    cdef double _intergenic_mod(const _node* n1, const _node* n2, const double start_weight)

    cdef void _score_connections(
        const uint8_t*   skip_connection,
        const uint8_t*   node_types,
        const int8_t*    node_strands,
              _node*     nodes,
        const int        min,
        const int        i,
        const _training* tinf,
        const int        final
    )

    # ctypedef void (*connection_function)(
    #     const _node*,
    #     const _node*,
    #     _node*,
    #     const _training*,
    #     const bint,
    # )
    # cdef connection_function CONNECTION_FUNCTIONS[4]

    cdef void _score_connection_forward_start(
        const _node*     nodes,
        const _node*     n1,
              _node*     n2,
        const _training* tinf,
        const bint       final,
    )
    cdef void _score_connection_forward_stop(
        const _node*     nodes,
        const _node*     n1,
              _node*     n2,
        const _training* tinf,
        const bint       final,
    )
    cdef void _score_connection_backward_start(
        const _node*     nodes,
        const _node*     n1,
              _node*     n2,
        const _training* tinf,
        const bint       final,
    )
    cdef void _score_connection_backward_stop(
        const _node*     nodes,
        const _node*     n1,
              _node*     n2,
        const _training* tinf,
        const bint       final,
    )
