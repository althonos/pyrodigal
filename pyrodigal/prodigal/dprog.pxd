from pyrodigal.prodigal.node cimport _node
from pyrodigal.prodigal.training cimport _training


cdef extern from "dprog.h" nogil:

    const int MAX_SAM_OVLP  = 60
    const int MAX_OPP_OVLP  = 200
    const int MAX_NODE_DIST = 500

    int dprog(_node*, int, _training*, int)
    void score_connection(_node*, int, int, _training*, int)
    void eliminate_bad_genes(_node*, int, _training*)
