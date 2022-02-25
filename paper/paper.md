---
title: 'Pyrodigal: Python bindings and interface to Prodigal, an efficient ORF finder for prokaryotes.'
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
date: 24 February 2022
bibliography: paper.bib
---

# Summary

Improvements in sequencing technologies have seen the amount of available
genomic data expand considerably over the last twenty years. One of the key
step for analysing this data is the prediction of protein-coding regions in
the genomic sequences, known as Open Reading Frames (ORFs), which span between
a start and a stop codon.

A recent comparison of several ORF prediction methods [@AssessORF:2019]
has shown that Prodigal [@Prodigal:2010], a prokaryotic gene finder using
dynamic programming, is one of the highest performing *ab initio* ORF finder.


# Statement of need

Prodigal is used in thousands of applications as the gene calling method
for processing genomic sequences. It is implemented in ANSI C, making it
extremely efficient and relatively easy to compile on different platforms.
However, the only way to use Prodigal is through an executable making it
inconvenient for less knowledgeable bioinformaticians only familiar with
the Python language.

In this work we present Pyrodigal, a Python module implemented in Cython
that binds to the Prodigal internals, resulting in identical predictions and
similar performance through a friendly object-oriented interface. The
predicted genes are returned as Python objects, with properties for retrieving
confidence scores or coordinates, and methods for translating the gene sequence.

`pyrodigal` has already been used in several works in preparation [@hafeZ:2021, @GECCO:2021]
and scientific publications [@PhageBoost:2021].


# Improvements

Although using the same data structures and scoring method as Prodigal,
Pyrodigal improves on several aspects of the original software by
re-implementing peripheral parts of the original software:

- Sequence data is not stored as a bitmap but as a byte array, which comes at
  the cost of slightly increased memory consumption in exchange of faster memory
  access. Memory profiling has revealed that the bulk of memory consumption
  in Prodigal is not caused by sequence data but node data, so this trade-off
  is acceptable.
- Memory management has been reworked to use buffers growing on insertion,
  making the memory consumption slightly more conservative on smaller sequences
  and addressing edges cases on sequences with a large number of start and
  STOP codons.
- Use of local buffers allows for thread-locality in the `OrfFinder.find_genes`
  method. In addition, `pyrodigal` makes use of the Cython feature for
  releasing the Python Global Interpreter Lock (GIL), which makes the ORF
  finder class usable in multi-threaded code.


# Optimization

In addition to the changes to the data storage, Pyrodigal implements a
SIMD-accelerated heuristic filter for skipping the scoring of connections
between dynamic programming nodes.

Prodigal works internally with a dynamic programming approach to score all the
connections nodes, which represent start and stop codon throughout the sequence.
The dynamic programming algorithm assigns a score for a gene spanning between
two codons based mostly on the frequency of nucleotide hexamers inside the
gene sequence. There are however pairs of codons between which a gene can
never span, such as two stop codons, or a forward start codon and a reverse
stop codon.

Identifying these invalid connections is done by checking the strand, type and
frame of the pair of nodes. Considering two nodes $i$ and $j$, the connection
between them is invalid if any of these boolean equations evaluates to true:

- $(type_i = STOP) \and (type_j \ne STOP) \and (strand_i = strand_j)$
- $(strand_i = 1) \and (type_i \ne STOP) \and (strand_j = -1)$
- $(strand_i = -1) \and (type_i = STOP) \and (strand_j = 1)$
- $(strand_i = -1) \and (type_i \ne STOP) \and (strand_j = 1) \and (type_j = STOP)$
- $(strand_i = strand_j) \and (strand_i = 1) \and (type_i \ne STOP) \and (type_j = STOP) \and (frame_i \ne frame_j)$
- $(strand_i = strand_j) \and (strand_i = -1) \and (type_i = STOP) \and (type_j \ne STOP) \and (frame_i \ne frame_j)$

<!-- Since all these attributes have a small number of possible values ($+1$ or $-1$ for the strand;
$ATG$, $GTG$, $TTG$ or $STOP$ for the type; $-1$, $-2$, $-3$, $+1$, $+2$, $+3$ for the frame),
they can each be stored in an single byte. -->

Before scoring the connections between node $i$ and all subsequent nodes,
Pyrodigal makes a first pass to check which connections are invalid based on the
equations above. Using SIMD features of modern CPUs allows processing several
nodes at once (8 nodes with NEON and SSE features, 16 nodes with AVX). This
first pass produces a look-up table that allows bypassing the scoring
of invalid connection.


# Distribution

Pyrodigal is distributed on the Python Package Index (PyPI), which make it
extremely simple to install with `pip`. Pre-compiled distributions are made
available for MacOS, Linux and Windows x86-64, as well as Linux Aarch64 machines.
Overall, distribution through PyPI makes it easier to develop a Python
workflow or application relying on Pyrodigal for the gene calling, since
end-users are not required to handle the setup of the Prodigal binaries
themselves.


# Acknowledgments

This work was funded by the European Molecular Biology Laboratory and the
German Research Foundation (Deutsche Forschungsgemeinschaft, DFG, grant no. 395357507 – SFB 1371).


# References
