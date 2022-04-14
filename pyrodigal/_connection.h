#ifndef _PYRODIGAL_CONNECTION_H
#define _PYRODIGAL_CONNECTION_H

#include <stdlib.h>
#include <stdint.h>
#include "node.h"
#include "training.h"


static inline double _intergenic_mod_diff(
    const struct _node* n1,
    const struct _node* n2,
    const double        start_weight
) {
    return -0.15 * start_weight;
}


static inline double _intergenic_mod_same(
    const struct _node* n1,
    const struct _node* n2,
    const double        start_weight
) {
    int    dist    = abs(n1->ndx - n2->ndx);
    int    overlap = n1->ndx + 2*n1->strand >= n2->ndx;
    double rval    = 0.0;

    if ((n1->ndx + 2 == n2->ndx) || (n1->ndx == n2->ndx + 1)) {
        if (n1->strand == 1) {
            if (n2->rscore < 0) rval -= n2->rscore;
            if (n2->uscore < 0) rval -= n2->uscore;
        } else {
            if (n1->rscore < 0) rval -= n1->rscore;
            if (n1->uscore < 0) rval -= n1->uscore;
        }
    }

    if (dist > 3 * OPER_DIST) {
        rval -= 0.15 * start_weight;
    } else if (((dist <= OPER_DIST) && !overlap) || (dist * 4 < OPER_DIST)) {
        rval += (2.0 - ((double) dist / OPER_DIST)) * 0.15 * start_weight;
    }

    return rval;
}


static inline double _intergenic_mod(
    const struct _node* n1,
    const struct _node* n2,
    const double        start_weight
) {
    if (n1->strand == n2->strand) {
        return _intergenic_mod_same(n1, n2, start_weight);
    } else {
        return _intergenic_mod_diff(n1, n2, start_weight);
    }
}


static void _score_connection_forward_start(
    const struct _node*              nodes,
    const struct _node*     restrict n1,
          struct _node*     restrict n2,
    const struct _training*          tinf,
    const int                        final
) {
    int ovlp  = 0;
    int maxfr = -1;
    int left  = n1->ndx;
    int right = n2->ndx;

    double score   = 0.0;
    double scr_mod = 0.0;

    // --- Edge Artifacts ---
    if ((n1->traceb == -1) && (n1->strand == 1) && (n1->type == STOP)) {
        return;
    } else if ((n1->traceb == -1) && (n1->strand == -1) && (n1->type != STOP)) {
        return;
    }

    // --- Intergenic space (Noncoding) ---
    // 3'fwd->5'fwd
    if ((n1->strand == 1) && (n1->type == STOP)) {
        left += 2;
        if (left >= right)
            return;
        if (final)
            score = _intergenic_mod_same(n1, n2, tinf->st_wt);
    // 5'rev->5'fwd
    } else if ((n1->strand == -1) && (n1->type != STOP)) {
        if (left >= right)
            return;
        if (final)
            score = _intergenic_mod_diff(n1, n2, tinf->st_wt);
    }

    if (!final) {
        score = ((double) (right - left + 1 - ovlp*2)) * scr_mod;
    }
    if (n1->score + score >= n2->score) {
        n2->score = n1->score + score;
        n2->traceb = n1 - nodes;
        n2->ov_mark = maxfr;
    }
}


static void _score_connection_forward_stop(
    const struct _node*              nodes,
    const struct _node*     restrict n1,
          struct _node*     restrict n2,
    const struct _training*          tinf,
    const int                        final
) {
    const struct _node* n3;

    int ovlp  = 0;
    int maxfr = -1;
    int left  = n1->ndx;
    int right = n2->ndx;

    double score   = 0.0;
    double scr_mod = 0.0;

    // --- Edge Artifacts ---
    if ((n1->traceb == -1) && (n1->strand == 1) && (n1->type == STOP)) {
        return;
    } else if ((n1->traceb == -1) && (n1->strand == -1) && (n1->type != STOP)) {
        return;
    }

    // --- Genes ---
    // 5'fwd->3'fwd
    if ((n1->strand == 1) && (n1->type != STOP)) {
        if (n2->stop_val >= n1->ndx)
            return;
        right += 2;
        if (final)
            score = n1->cscore + n1->sscore;
        else
            scr_mod = tinf->bias[0]*n1->gc_score[0] + tinf->bias[1]*n1->gc_score[1] + tinf->bias[2]*n1->gc_score[2];

    // --- Possible Operons */ ---
    // 3'fwd->3'fwd, check for a start just to left of first 3'
    } else if ((n1->strand == 1) && (n1->type == STOP)) {
        if (n2->stop_val >= n1->ndx)
            return;
        if (n1->star_ptr[n2->ndx%3] == -1)
            return;
        n3 = &nodes[n1->star_ptr[n2->ndx%3]];
        left = n3->ndx;
        right += 2;
        if (final)
            score = n3->cscore + n3->sscore + _intergenic_mod(n1, n3, tinf->st_wt);
        else
            scr_mod = tinf->bias[0]*n3->gc_score[0] + tinf->bias[1]*n3->gc_score[1] + tinf->bias[2]*n3->gc_score[2];
    }

    if (!final) {
        score = ((double) (right - left + 1 - ovlp*2)) * scr_mod;
    }
    if (n1->score + score >= n2->score) {
        n2->score = n1->score + score;
        n2->traceb = n1 - nodes;
        n2->ov_mark = maxfr;
    }
}


static void _score_connection_backward_start(
    const struct _node*              nodes,
    const struct _node*     restrict n1,
          struct _node*     restrict n2,
    const struct _training*          tinf,
    const int                        final
) {
    int bnd;
    int ovlp  = 0;
    int maxfr = -1;
    int left  = n1->ndx;
    int right = n2->ndx;

    double score   = 0.0;
    double scr_mod = 0.0;

    // --- Edge Artifacts ---
    if ((n1->traceb == -1) && (n1->strand == 1) && (n1->type == STOP)) {
        return;
    } else if ((n1->traceb == -1) && (n1->strand == -1) && (n1->type != STOP)) {
        return;
    }

    // --- Genes ---
    // 3'rev->5'rev
    if ((n1->strand == -1) && (n1->type == STOP)) {
        if (n1->stop_val <= n2->ndx)
            return;
        left -= 2;
        if (final)
            score = n2->cscore + n2->sscore;
        else
            scr_mod = tinf->bias[0]*n2->gc_score[0] + tinf->bias[1]*n2->gc_score[1] + tinf->bias[2]*n2->gc_score[2];

    // --- Overlapping Opposite Strand 3' Ends ---
    // 3'for->5'rev
    } else if ((n1->strand == 1) && (n1->type == STOP)) {
        if (n2->stop_val - 2 >= n1->ndx + 2)
            return;
        ovlp = (n1->ndx+2) - (n2->stop_val-2) + 1;
        if (ovlp >= MAX_OPP_OVLP)
            return;
        if ((n1->ndx - n2->stop_val) >= (n2->ndx - n1->ndx + 3))
            return;
        bnd = (n1->traceb == -1) ? 0 : nodes[n1->traceb].ndx;
        if ((n1->ndx - n2->stop_val) >= (n2->stop_val - 3 - bnd))
            return;
        left = n2->stop_val-2;
        if (final)
            score = n2->cscore + n2->sscore + _intergenic_mod_diff(n1, n2, tinf->st_wt);
        else
            scr_mod = tinf->bias[0]*n2->gc_score[0] + tinf->bias[1]*n2->gc_score[1] + tinf->bias[2]*n2->gc_score[2];
    }

    if (!final) {
        score = ((double) (right - left + 1 - ovlp*2)) * scr_mod;
    }
    if (n1->score + score >= n2->score) {
        n2->score = n1->score + score;
        n2->traceb = n1 - nodes;
        n2->ov_mark = maxfr;
    }
}


static void _score_connection_backward_stop(
    const struct _node*              nodes,
    const struct _node*     restrict n1,
          struct _node*     restrict n2,
    const struct _training*          tinf,
    const int                        final
) {
    const struct _node* n3;

    int i;
    int ovlp  = 0;
    int maxfr = -1;
    int left  = n1->ndx;
    int right = n2->ndx;

    double maxval;
    double score   = 0.0;
    double scr_mod = 0.0;

    // --- Edge Artifacts ---
    if ((n1->traceb == -1) && (n1->strand == 1) && (n1->type == STOP)) {
        return;
    } else if ((n1->traceb == -1) && (n1->strand == -1) && (n1->type != STOP)) {
        return;
    }

    // --- Intergenic space (Noncoding) ---
    // 3'fwd->3'rev
    if ((n1->strand == 1) && (n1->type == STOP)) {
        left += 2;
        right -= 2;
        if (left >= right)
            return;
        // Overlapping Gene Case 2: Three consecutive overlapping genes f r r
        maxfr = -1;
        maxval = 0.0;
        if (n1->traceb != -1) {
            for (i = 0; i < 3; i++) {
                if (n2->star_ptr[i] == -1)
                    continue;
                n3 = &nodes[n2->star_ptr[i]];
                ovlp = left - n3->stop_val + 3;
                if ((ovlp <= 0) || (ovlp >= MAX_OPP_OVLP))
                    continue;
                if (ovlp >= n3->ndx - left)
                    continue;
                if (ovlp >= n3->stop_val - nodes[n1->traceb].ndx - 2)
                    continue;
                // record max frame
                if (final)
                    score = n3->cscore + n3->sscore + _intergenic_mod(n3, n2, tinf->st_wt);
                else
                    score = tinf->bias[0]*n3->gc_score[0] + tinf->bias[1]*n3->gc_score[1] + tinf->bias[2]*n3->gc_score[2];
                if (score > maxval) {
                    maxfr = i;
                    maxval = score;
                }
            }
        }
        if (maxfr != -1) {
            n3 = &nodes[n2->star_ptr[maxfr]];
            if (final)
                score = n3->cscore + n3->sscore + _intergenic_mod(n3, n2, tinf->st_wt);
            else
                scr_mod = tinf->bias[0]*n3->gc_score[0] + tinf->bias[1]*n3->gc_score[1] + tinf->bias[2]*n3->gc_score[2];
        } else if (final) {
            score = _intergenic_mod_diff(n1, n2, tinf->st_wt);
        }
    // 5'rev->3'rev
    } else if ((n1->strand == -1) && (n1->type != STOP)) {
        right -= 2;
        if (left >= right)
            return;
        if (final)
            score = _intergenic_mod_same(n1, n2, tinf->st_wt);

    // --- Possible Operons */ ---
    // 3'rev->3'rev, check for a start just to right of second 3'
    } else if ((n1->strand == -1) && (n1->type == STOP)) {
        if (n1->stop_val <= n2->ndx)
            return;
        if (n2->star_ptr[n1->ndx%3] == -1)
            return;
        n3 = &nodes[n2->star_ptr[n1->ndx%3]];
        left -= 2;
        right = n3->ndx;
        if (final)
            score = n3->cscore + n3->sscore + _intergenic_mod(n3, n2, tinf->st_wt);
        else
            scr_mod = tinf->bias[0]*n3->gc_score[0] + tinf->bias[1]*n3->gc_score[1] + tinf->bias[2]*n3->gc_score[2];
    }

    if (!final) {
        score = ((double) (right - left + 1 - ovlp*2)) * scr_mod;
    }
    if (n1->score + score >= n2->score) {
        n2->score = n1->score + score;
        n2->traceb = n1 - nodes;
        n2->ov_mark = maxfr;
    }
}


typedef void(*connection_function)(
    const struct _node*,
    const struct _node*,
          struct _node*,
    const struct _training*,
    const int
);

static connection_function CONNECTION_FUNCTIONS[4] = {
    _score_connection_forward_start,
    _score_connection_forward_stop,
    _score_connection_backward_start,
    _score_connection_backward_stop,
};


static inline void _score_connections(
    const uint8_t*          skip_connection,
    const uint8_t*          node_types,
    const int8_t*           node_strands,
          struct _node*     nodes,
    const int               min,
    const int               i,
    const struct _training* tinf,
    const int               final
) {
    int j;
    int kind;
    kind = 2*(node_strands[i] == -1) + 1*(node_types[i] == STOP);
    for (j = min; j < i; j++)
        if (!skip_connection[j])
            CONNECTION_FUNCTIONS[kind](nodes, &nodes[j], &nodes[i], tinf, final);
}


#endif
