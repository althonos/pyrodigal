from libc.stdio cimport FILE

from pyrodigal.prodigal.bitmap cimport bitmap_t
from pyrodigal.prodigal.node cimport _node
from pyrodigal.prodigal.training cimport _training


cdef extern from "gene.h" nogil:

    const size_t MAX_GENES = 30000

    struct _gene:
        int begin
        int end
        int start_ndx
        int stop_ndx
        char gene_data[500]
        char score_data[500]

    int add_genes(_gene*, _node*, int)
    void record_gene_data(_gene*, int, _node*, _training*, int)
    void tweak_final_starts(_gene*, int, _node*, int, _training*)

    void write_translations(FILE *fh, _gene* genes, int ng, _node* nod, bitmap_t seq, bitmap_t rseq, bitmap_t useq, int slen, _training* tinf, int sctr, char* short_hdr)
    void print_genes(FILE*, _gene*, int, _node*, int, int, int, int, char*, _training*, char*, char*, char*)
