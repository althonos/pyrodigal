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
genomic data expand considerably in the last twenty years. One of the key
step for analysing this data is the prediction of protein-coding regions in
the genomic sequences, known as Open Reading Frames (ORFs), which span between
a start and a STOP codon.


# Statement of need

Prodigal [@Prodigal:2010] is an ORF finder for prokaryotes used in thousands
of applications as the gene calling method for processing genomic sequences.
It is implemented in ANSI C, making it extremely efficient and relatively
easy to compile on different platforms. However, the only way to use Prodigal
is through an executable making it less convenient for less knowledgeable
bioinformaticians only familiar with the Python language.

`pyrodigal` is a Python module implemented in Cython that binds to the
Prodigal interinals, resulting in identical predictions and similar performance
through a friendly object-oriented interface. It removes the need for the
Prodigal binary to be installed on the host machine. The predicted genes are
returned as Python objects, with properties for retrieving confidence scores
or coordinates, and methods for translating the gene sequence.

`pyrodigal` has already been used in several works in preparation [@hafeZ:2021, @GECCO:2021]
and scientific publications [@PhageBoost:2021].


# Performance

Although using the same data structures and scoring method as Prodigal,
`pyrodigal` improves on several aspects of the original software by
re-implementing key parts of the original software:

- Sequence data is not stored as a bitmap anymore, which comes at the cost
  of slightly increased memory consumption in exchange of faster memory
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

In addition to these changes, `pyrodigal` implements a SIMD-accelerated
heuristic filter for scoring connections between dynamic programming nodes.


# Distribution

`pyrodigal` is distributed on the Python Package Index (PyPI), which make it
extremely simple to install. Pre-compiled distributions are made available
for MacOS, Linux and Windows x86-64, as well as Linux Aarch64 machines.
Overall, distribution through PyPI makes it easier to develop a Python
workflow or application relying on `pyrodigal` for the gene calling, since
end-users are not required to handle the setup of the Prodigal binaries
themselves.


# Acknowledgments

This work was funded by the European Molecular Biology Laboratory and the
German Research Foundation (Deutsche Forschungsgemeinschaft, DFG, grant no. 395357507 – SFB 1371).
