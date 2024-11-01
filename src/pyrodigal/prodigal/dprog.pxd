from pyrodigal.prodigal.node cimport _node
from pyrodigal.prodigal.training cimport _training


cdef extern from "dprog.h" nogil:

    const int MAX_SAM_OVLP
    const int MAX_OPP_OVLP
    const int MAX_NODE_DIST

    int dprog(_node*, int, _training*, int) noexcept
    void score_connection(_node*, int, int, _training*, int) noexcept
    void eliminate_bad_genes(_node*, int, _training*) noexcept
