from pyrodigal.prodigal.node cimport _node
from pyrodigal.prodigal.training cimport _training


cdef extern from "dprog.h" nogil:
    int dprog(_node*, int, _training*, int)
    void score_connection(_node*, int, int, _training*, int)
    void eliminate_bad_genes(_node*, int, _training*)
