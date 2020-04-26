from pyrodigal.prodigal.node cimport _node
from pyrodigal.prodigal.training cimport _training


cdef extern from "dprog.h":
    cdef int dprog(_node*, int, _training*, int)
    cdef void score_connection(_node*, int, int, _training*, int);
    cdef void eliminate_bad_genes(_node*, int, _training*);
