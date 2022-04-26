---
title: 'Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes'
tags:
  - Python
  - Cython
  - bioinformatics
  - gene prediction
authors:
  - name: Martin Larralde
    orcid: 0000-0002-3947-4444
    affiliation: 1
affiliations:
  - name: Structural and Computational Biology Unit, EMBL, Heidelberg, Germany
    index: 1
date: 24 February 2022
bibliography: paper.bib
---

# Summary

Improvements in sequencing technologies have seen the amount of available
genomic data expand considerably over the last twenty years. One of the key
steps for analysing is the prediction of protein-coding regions in genomic
sequences, known as Open Reading Frames (ORFs), which span between
a start and a stop codon. A recent comparison of several ORF prediction methods
[@AssessORF:2019] has shown that Prodigal [@Prodigal:2010], a prokaryotic gene
finder that uses dynamic programming, is one of the highest performing
*ab initio* ORF finders. Pyrodigal is a Python package that provides bindings
and an interface to Prodigal to make it easier to use in Python applications.


# Statement of need

Prodigal is used in thousands of applications as the gene calling method
for processing genomic sequences. It is implemented in ANSI C, making it
extremely efficient and relatively easy to compile on different platforms.
However, the only way to use Prodigal is through an executable, making it
inconvenient for bioinformaticians who rely on the increasingly popular
Python language. In addition to the hassle caused by the invocation of
the executable from Python code, the distribution of Python programs relying
on Prodigal is also problematic, since they now require an external binary
that cannot be installed from the [Python Package Index](https://pypi.org) (PyPI).

To address these issues, we developed Pyrodigal, a Python module implemented
in Cython [@Cython] that binds to the Prodigal internals, resulting in identical
predictions and similar performance through a friendly object-oriented interface.
The predicted genes are returned as Python objects, with properties for retrieving
confidence scores or coordinates, and methods for translating the gene sequence
with a default or user-provided translation table. Prodigal is compiled from
source and statically linked into a compiled Python extension, which allows
it to be installed with a single `pip install` command, even on a target
machine that requires compilation.

Pyrodigal has already been used as the implementation for the initial ORF
finding stage in several domains, including biosynthetic gene cluster
prediction [@GECCO:2021], prophage identification [@PhageBoost:2021; @hafeZ:2021],
and pangenome analysis [@AlphaMine:2021].


# Method

Internally, Prodigal identifies start and stop codons throughout a genomic
sequence, which are represented as nodes. It then uses a dynamic programming
approach to compute scores for every pair of nodes (i.e. putative genes) as
shown in \autoref{fig:method}. Scores assigned to predicted genes are based
primarily on the frequency of nucleotide hexamers inside the gene sequence.

![Graphical depiction of the Prodigal method for identifying genes in a sequence.
First, the sequence is analysed to find start and stop codons in the 6 reading frames (1).
Dynamic programming nodes are then created for each codon, each storing the strand
(shown with the triangle direction, right for the direct strand, left for the reverse strand)
and the type (green for start, red for stop) of the codon they were obtained from.
Nodes are then scored on biological criteria (2). Putative genes are identified between all
start and stop codons in a given window (a gene cannot span between any pair of nodes;
some invalid connections are shown with dashed lines as examples). Then putative
genes are scored using a dynamic programming approach (3). Once all connections
have been processed, the dynamic programming matrix is traversed to find the highest scoring
path, giving the final predictions (4). \label{fig:method}](figure1.svg){width=100%}

Pyrodigal adapts the first two steps so that the dynamic programming nodes
can be extracted directly from a Python string containing the sequence data,
rather than requiring formatting to an external file. In addition, the node
storage has been reworked to use reallocating buffers, saving memory on
smaller sequences. The node and connection scoring steps use the original
Prodigal code.


# Optimization

Prodigal was profiled with Valgrind [@Valgrind:2007] to identify critical
parts of the original code. Using bacterial genomes from the proGenomes v2.1
database [@proGenomes2:2020], we found that for long enough sequences,
a large fraction of the CPU cycles was spent in the `score_connection`
function.

There are pairs of codons between which a gene can never span, such as two
stop codons, or a forward start codon and a reverse stop codon, as shown in
\autoref{fig:method}. Upon inspection, we realized the `score_connection`
was still called in invalid cases that could be labelled as such beforehand.
Identifying these invalid connections is feasible by checking the strand, type
and reading frame of a node pair. Considering two nodes $i$ and $j$, the
connection between them is invalid if any of these boolean equations is true:

- $(T_i \ne STOP) \land (T_j \ne STOP) \land (S_i = S_j)$
- $(S_i = 1) \land (T_i \ne STOP) \land (S_j = -1)$
- $(S_i = -1) \land (T_i = STOP) \land (S_j = 1)$
- $(S_i = -1) \land (T_i \ne STOP) \land (S_j = 1) \land (T_j = STOP)$
- $(S_i = S_j) \land (S_i = 1) \land (T_i \ne STOP) \land (T_j = STOP) \land (F_i \ne F_j)$
- $(S_i = S_j) \land (S_i = -1) \land (T_i = STOP) \land (T_j \ne STOP) \land (F_i \ne F_j)$

where $T_i$, $S_i$ and $F_i$ are respectively the type, strand, and reading
frame of the node $i$, with forward and reverse strands encoded as $+1$ and
$-1$ respectively.

We developed a heuristic filter to quickly identify node pairs forming an
invalid connection prior to the scoring step using the above formulas.
Since all these attributes have a small number of possible values
($+1$ or $-1$ for the forward or reverse strand; $ATG$, $GTG$, $TTG$ or $STOP$
for the codon type; $-1$, $-2$, $-3$, $+1$, $+2$, $+3$ for the reading frame),
they can all be stored in a single byte. Use of the SIMD features of modern CPUs
allows several nodes to be processed at once (8 nodes with NEON and SSE2 features,
16 nodes with AVX2). This first pass produces a look-up table used to bypass the
scoring of invalid connections.

The performance of the connection scoring was evaluated on 50 bacterial
sequences of various length, as shown in \autoref{fig:benchmark}. It suggests
that even with the added cost of the additional pass for each node, enabling
the heuristic filter in Pyrodigal saves about half of the time needed to score
connections between all the nodes of a sequence.

![Evaluation of the connection scoring performance with different heuristic
filter SIMD backends (SSE2 or AVX2), with a generic backend (Generic) or without enabling 
the filter (None).
*Each sequence was processed 10 times on a quiet i7-10710U CPU @ 1.10GHz*. \label{fig:benchmark}](figure2.svg){width=100%}


# Availability

Pyrodigal is distributed on PyPI under the
[GNU Lesser General Public License](https://www.gnu.org/licenses/lgpl-3.0).
Pre-compiled distributions are provided for MacOS, Linux and Windows x86-64,
as well as Linux Aarch64 machines. A [Conda](https://conda.io/) package is
also available in the Bioconda channel [@Bioconda:2018].

The source code is available in a git repository on [GitHub](https://github.com/althonos/pyrodigal),
and features a Continuous Integration workflow to run integration tests on changes.
Documentation is hosted on [ReadTheDocs](https://pyrodigal.readthedocs.io) and
built for each new release.


# Acknowledgments

We thank Laura M. Carroll for her input on the redaction of this manuscript,
and Georg Zeller for his supervision. This work was funded by the European
Molecular Biology Laboratory and the German Research Foundation
(Deutsche Forschungsgemeinschaft, DFG, grant no. 395357507 – [SFB 1371](https://www.sfb1371.tum.de/)).


# References
