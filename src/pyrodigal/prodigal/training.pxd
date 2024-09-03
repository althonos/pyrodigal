cdef extern from "training.h":

    struct _training:
        double gc
        int trans_table
        double st_wt
        double bias[3]
        double type_wt[3]
        int uses_sd
        double rbs_wt[28]
        double ups_comp[32][4]
        double mot_wt[4][4][4096]
        double no_mot
        double gene_dc[4096]
