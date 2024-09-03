/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#ifndef _NODE_H
#define _NODE_H

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "sequence.h"
#include "training.h"

#define STT_NOD 100000
#define MIN_GENE 90
#define MIN_EDGE_GENE 60
#define MAX_SAM_OVLP 60
#define ST_WINDOW 60
#define OPER_DIST 60
#define EDGE_BONUS 0.74
#define EDGE_UPS -1.00
#define META_PEN 7.5

struct _motif {
  double score;        /* Score for the motif */
  uint16_t ndx: 12;    /* 6mer of the best motif for this node */
  uint8_t spacer: 4;   /* Spacer between coding start and the motif */
  uint8_t len: 3;      /* Length of the motif */
  uint8_t spacendx: 2; /* Index for this spacer length */
};

struct _node {
  struct _motif mot;   /* Upstream motif information for this node */
  
  double gc_score[3];  /* % GC content in different codon positions */
  double cscore;       /* Coding score for this node (based on 6-mer usage) */
  double uscore;       /* Score for the upstream -1/-2, -15to-45 region */
  double tscore;       /* Score for the ATG/GTG/TTG value */
  double rscore;       /* Score for the RBS motif */
  double sscore;       /* Score for the strength of the start codon */
  double score;        /* Score of total solution to this point */
  float gc_cont;      /* GC Content for the node */

  int star_ptr[3];     /* Array of starts w/in MAX_SAM_OVLP bases of a stop in 3 frames */ 
  int traceb;          /* Traceback to connecting node */
  int tracef;          /* Forward trace */ 
  int ndx;             /* position in the sequence of the node */
  int stop_val;        /* For a stop, record previous stop; for start, record its stop */

  int8_t  ov_mark;     /* Marker to help untangle overlapping genes */
  int8_t  strand;      /* 1 = forward, -1 = reverse */
  uint8_t rbs[2];     /* SD RBS score for this node (based on binding energy)
                          rbs[0] = best motif with exact match, rbs[1] = with
                          mismatches */
  
  uint8_t edge;     /* Runs off the edge; 0 = normal, 1 = edge node */
  uint8_t elim;     /* If set to 1, eliminate this gene from the model */
  uint8_t gc_bias;  /* Frame of highest GC content within this node */
  uint8_t type;     /* 0=ATG, 1=GTG, 2=TTG/Other, 3=Stop */
};

int add_nodes(unsigned char *, unsigned char *, int, struct _node *, int,
              mask *, int, struct _training *);
void reset_node_scores(struct _node *, int);
int compare_nodes(const void *, const void *);
int stopcmp_nodes(const void *, const void *);

void record_overlapping_starts(struct _node *, int, struct _training *, int);
void record_gc_bias(int *, struct _node *, int, struct _training *);

void calc_dicodon_gene(struct _training *, unsigned char *, unsigned char *,
                       int, struct _node *, int);
void calc_amino_bg(struct _training *, unsigned char *, unsigned char *, int,
                   struct _node *, int);

void score_nodes(unsigned char *, unsigned char *, int, struct _node *, int,
                 struct _training *, int, int);
void raw_coding_score(unsigned char *, unsigned char *, int, struct _node *,
                      int, struct _training *);
void calc_orf_gc(unsigned char *, unsigned char *, int, struct _node *, int, 
                 struct _training *);
void rbs_score(unsigned char *, unsigned char *, int, struct _node *, int,
               struct _training *);
void score_upstream_composition(unsigned char *, int, struct _node *, 
                                struct _training *);

void determine_sd_usage(struct _training *);

double intergenic_mod(struct _node *, struct _node *, struct _training *);

void train_starts_sd(unsigned char *, unsigned char *, int, struct _node *,
                        int, struct _training *);
void train_starts_nonsd(unsigned char *, unsigned char *, int, struct _node *,
                        int, struct _training *);

void count_upstream_composition(unsigned char *, int, int, int, 
                                struct _training *);

void build_coverage_map(double [4][4][4096], int [4][4][4096], double, int);
void find_best_upstream_motif(struct _training *, unsigned char *, unsigned
                              char *, int, struct _node *, int);
void update_motif_counts(double [4][4][4096], double *, unsigned char *,
                         unsigned char *, int, struct _node *, int);

void write_start_file(FILE *, struct _node *, int, struct _training *, int,
                      int, int, char *, char *, char *);

int cross_mask(int, int, mask *, int);

double dmax(double, double);
double dmin(double, double);

#endif
