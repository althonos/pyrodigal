# 🔥 Pyrodigal [![Stars](https://img.shields.io/github/stars/althonos/pyrodigal.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyrodigal/stargazers)

*Cython bindings and Python interface to [Prodigal](https://github.com/hyattpd/Prodigal/), an ORF
finder for genomes and metagenomes. **Now with SIMD!***

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pyrodigal/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pyrodigal/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyrodigal?style=flat-square&maxAge=3600&logo=codecov)](https://codecov.io/gh/althonos/pyrodigal/)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![PyPI](https://img.shields.io/pypi/v/pyrodigal.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pyrodigal)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pyrodigal?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/pyrodigal)
[![AUR](https://img.shields.io/aur/version/python-pyrodigal?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pyrodigal)
[![Wheel](https://img.shields.io/pypi/wheel/pyrodigal.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyrodigal/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyrodigal.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pyrodigal/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyrodigal.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pyrodigal/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyrodigal/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/pyrodigal.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyrodigal/issues)
[![Docs](https://img.shields.io/readthedocs/pyrodigal/latest?style=flat-square&maxAge=600)](https://pyrodigal.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyrodigal/blob/main/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/pyrodigal?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/pyrodigal)
[![Paper](https://img.shields.io/badge/paper-JOSS-9400ff?style=flat-square&maxAge=86400)](https://doi.org/10.21105/joss.04296)
[![Citations](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fbadge.dimensions.ai%2Fdetails%2Fid%2Fpub.1147419140%2Fmetadata.json&query=%24.times_cited&style=flat-square&label=citations&maxAge=86400)](https://badge.dimensions.ai/details/id/pub.1147419140)

## 🗺️ Overview

Pyrodigal is a Python module that provides bindings to Prodigal using
[Cython](https://cython.org/). It directly interacts with the Prodigal
internals, which has the following advantages:

- **single dependency**: Pyrodigal is distributed as a Python package, so you
  can add it as a dependency to your project, and stop worrying about the
  Prodigal binary being present on the end-user machine.
- **no intermediate files**: Everything happens in memory, in a Python object
  you fully control, so you don't have to invoke the Prodigal CLI using a
  sub-process and temporary files. Sequences can be passed directly as
  strings or bytes, which avoids the overhead of formatting your input to
  FASTA for Prodigal.
- **better memory usage**: Pyrodigal uses more compact data structures compared
  to the original Prodigal implementation, allowing to save memory to store 
  the same information. A heuristic is used to estimate the number of nodes
  to allocate based on the sequence GC% in order to minimize reallocations.
- **better performance**: Pyrodigal uses *SIMD* instructions to compute which
  dynamic programming nodes can be ignored when scoring connections. This can
  save from a third to half the runtime depending on the sequence. The [Benchmarks](https://pyrodigal.readthedocs.io/en/stable/benchmarks.html) page of the documentation contains comprehensive comparisons. See the [JOSS paper](https://doi.org/10.21105/joss.04296)
  for details about how this is achieved.
- **same results**: Pyrodigal is tested to make sure it produces
  exactly the same results as Prodigal `v2.6.3+31b300a`. *This was verified
  extensively by [Julian Hahnfeld](https://github.com/jhahnfeld) and can be
  checked with his [comparison repository](https://github.com/jhahnfeld/prodigal-pyrodigal-comparison).*

### 📋 Features

The library now features everything from the original Prodigal CLI:

- **run mode selection**: Choose between *single* mode, using a training
  sequence to count nucleotide hexamers, or *metagenomic* mode, using
  pre-trained data from different organisms (`prodigal -p`).
- **region masking**: Prevent genes from being predicted across regions
  containing unknown nucleotides  (`prodigal -m`).
- **closed ends**: Genes will be identified as running over edges if they
  are larger than a certain size, but this can be disabled (`prodigal -c`).
- **training configuration**: During the training process, a custom
  translation table can be given (`prodigal -g`), and the Shine-Dalgarno motif
  search can be forcefully bypassed (`prodigal -n`)
- **output files**: Output files can be written in a format mostly
  compatible with the Prodigal binary, including the protein translations
  in FASTA format (`prodigal -a`), the gene sequences in FASTA format
  (`prodigal -d`), or the potential gene scores in tabular format
  (`prodigal -s`).
- **training data persistence**: Getting training data from a sequence and
  using it for other sequences is supported; in addition, a training data
  file can be saved and loaded transparently (`prodigal -t`).

In addition, the **new** features are available:

- **custom gene size threshold**: While Prodigal uses a minimum gene size
  of 90 nucleotides (60 if on edge), Pyrodigal allows to customize this
  threshold, allowing for smaller ORFs to be identified if needed.
- **custom metagenomic models**: Since `v3.0.0`, you can use your own 
  metagenomic models to run Pyrodigal in *meta*-mode. *Check for instance
  [`pyrodigal-gv`](https://github.com/althonos/pyrodigal-gv), which 
  provides additional models for giant viruses and gut phages.*

### 🐏 Memory

Pyrodigal makes several changes compared to the original Prodigal binary
regarding memory management:

* Sequences are stored as raw bytes instead of compressed bitmaps. This means
  that the sequence itself takes 3/8th more space, but since the memory used
  for storing the sequence is often negligible compared to the memory used to
  store dynamic programming nodes, this is an acceptable trade-off for better
  performance when extracting said nodes.
* Node fields use smaller data types to fit into 128 bytes, compared to the 
  176 bytes of the original Prodigal data structure.
* Node arrays are pre-allocated based on the sequence GC% to extrapolate the
  probability to find a start or stop codon.
* Genes are stored in a more compact data structure than in Prodigal (which
  reserves a buffer to store string data), saving around 1KiB per gene.


### 🧶 Thread-safety

[`pyrodigal.GeneFinder`](https://pyrodigal.readthedocs.io/en/stable/api/orf_finder.html#pyrodigal.GeneFinder)
instances are thread-safe. In addition, the
[`find_genes`](https://pyrodigal.readthedocs.io/en/stable/api/gene_finder.html#pyrodigal.GeneFinder.find_genes)
method is re-entrant. This means you can train an
[`GeneFinder`](https://pyrodigal.readthedocs.io/en/stable/api/gene_finder.html#pyrodigal.GeneFinder)
instance once, and then use a pool to process sequences in parallel:
```python
import multiprocessing.pool
import pyrodigal

gene_finder = pyrodigal.GeneFinder()
gene_finder.train(training_sequence)

with multiprocessing.pool.ThreadPool() as pool:
    predictions = pool.map(orf_finder.find_genes, sequences)
```

## 🔧 Installing

Pyrodigal can be installed directly from [PyPI](https://pypi.org/project/pyrodigal/),
which hosts some pre-built wheels for the x86-64 architecture (Linux/MacOS/Windows)
and the Aarch64 architecture (Linux/MacOS), as well as the code required to compile
from source with Cython:
```console
$ pip install pyrodigal
```

Otherwise, Pyrodigal is also available as a [Bioconda](https://bioconda.github.io/)
package:
```console
$ conda install -c bioconda pyrodigal
```

Check the [*install* page](https://pyrodigal.readthedocs.io/en/stable/install.html)
of the documentation for other ways to install Pyrodigal on your machine.

## 💡 Example

Let's load a sequence from a
[GenBank](http://www.insdc.org/files/feature_table.html) file, use an `GeneFinder`
to find all the genes it contains, and print the proteins in two-line FASTA
format.

### 🔬 [Biopython](https://github.com/biopython/biopython)

To use the [`GeneFinder`](https://pyrodigal.readthedocs.io/en/stable/api/orf_finder.html#pyrodigal.GeneFinder)
in single mode (corresponding to `prodigal -p single`, the default operation mode of Prodigal),
you must explicitly call the
[`train`](https://pyrodigal.readthedocs.io/en/stable/api/orf_finder.html#pyrodigal.GeneFinder.train) method
with the sequence you want to use for training before trying to find genes,
or you will get a [`RuntimeError`](https://docs.python.org/3/library/exceptions.html#RuntimeError):
```python
import Bio.SeqIO
import pyrodigal

record = Bio.SeqIO.read("sequence.gbk", "genbank")

orf_finder = pyrodigal.GeneFinder()
orf_finder.train(bytes(record.seq))
genes = orf_finder.find_genes(bytes(record.seq))
```

However, in `meta` mode (corresponding to `prodigal -p meta`), you can find genes directly:
```python
import Bio.SeqIO
import pyrodigal

record = Bio.SeqIO.read("sequence.gbk", "genbank")

orf_finder = pyrodigal.GeneFinder(meta=True)
for i, pred in enumerate(orf_finder.find_genes(bytes(record.seq))):
    print(f">{record.id}_{i+1}")
    print(pred.translate())
```

*On older versions of Biopython (before 1.79) you will need to use
`record.seq.encode()` instead of `bytes(record.seq)`*.


### 🧪 [Scikit-bio](https://github.com/biocore/scikit-bio)

```python
import skbio.io
import pyrodigal

seq = next(skbio.io.read("sequence.gbk", "genbank"))

orf_finder = pyrodigal.GeneFinder(meta=True)
for i, pred in enumerate(orf_finder.find_genes(seq.values.view('B'))):
    print(f">{record.id}_{i+1}")
    print(pred.translate())
```

*We need to use the [`view`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.view.html)
method to get the sequence viewable by Cython as an array of `unsigned char`.*


## 🔖 Citation

Pyrodigal is scientific software, with a
[published paper](https://doi.org/10.21105/joss.04296)
in the [Journal of Open-Source Software](https://joss.theoj.org/). Please
cite both [Pyrodigal](https://doi.org/10.21105/joss.04296)
and [Prodigal](https://doi.org/10.1186/1471-2105-11-119) if you are using it in
an academic work, for instance as:

> Pyrodigal (Larralde, 2022), a Python library binding to Prodigal (Hyatt *et al.*, 2010).

Detailed references are available on the [Publications page](https://pyrodigal.readthedocs.io/en/stable/publications.html) of the
[online documentation](https://pyrodigal.readthedocs.io/).

## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/althonos/pyrodigal/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/pyrodigal/blob/main/CONTRIBUTING.md)
for more details.

## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/althonos/pyrodigal/blob/main/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.


## ⚖️ License

This library is provided under the [GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/).
The Prodigal code was written by [Doug Hyatt](https://github.com/hyattpd) and is distributed under the
terms of the GPLv3 as well. See `vendor/Prodigal/LICENSE` for more information.

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the [original Prodigal authors](https://github.com/hyattpd). It was developed
by [Martin Larralde](https://github.com/althonos/) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*
