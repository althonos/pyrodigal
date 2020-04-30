# 🔥 Pyrodigal

*Python interface to [Prodigal](https://github.com/hyattpd/Prodigal/), an ORF
finder for genomes, progenomes and metagenomes.*

[![TravisCI](https://img.shields.io/travis/althonos/pyrodigal/master.svg?logo=travis&maxAge=600&style=flat-square)](https://travis-ci.com/althonos/pyrodigal/branches)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyrodigal?style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/pyrodigal/)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![PyPI](https://img.shields.io/pypi/v/pyrodigal.svg?style=flat-square&maxAge=600)](https://pypi.org/project/pyrodigal)
[![Wheel](https://img.shields.io/pypi/wheel/pyrodigal.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyrodigal/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyrodigal.svg?style=flat-square&maxAge=600)](https://pypi.org/project/pyrodigal/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyrodigal.svg?style=flat-square&maxAge=600)](https://pypi.org/project/pyrodigal/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyrodigal/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/pyrodigal.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyrodigal/issues)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyrodigal.py/blob/master/CHANGELOG.md)

<!-- [![AppVeyor](https://img.shields.io/appveyor/ci/althonos/pyrodigal/master?logo=appveyor&style=flat-square&maxAge=600)](https://ci.appveyor.com/project/althonos/pyrodigal) -->


## 🗺️ Overview

### 📋 Features

The library now features everything needed to run Prodigal in metagenomic mode.
It does not yet support single mode, which requires more configuration from
the user but offers more flexibility.

**Roadmap**:

- [x] Metagenomic mode
- [ ] Single mode
- [ ] External training file support (`-t` flag)
- [ ] Region masking (`-m` flag)

### 📋 Memory

Contrary to the Prodigal command line, Pyrodigal attempts to be more conservative
about memory usage. This means that most of the allocations will be lazy, and
that some functions will reallocate their results to exact-sized arrays when
it's possible. This leads to Pyrodigal using about 30% less memory, but with
some more overhead

### 🧶 Thread-safety

`pyrodigal.Pyrodigal` instances are not thread-safe: concurrent `find_genes`
calls will overwrite the internal memory used for dynamic programming and
could lead to unexpected crashes. A solution to process sequences in parallel
is to use a consumer/worker pattern, and have on `Pyrodigal` instance in each
worker. Using a pool spawning `Pyrodigal` instances on the fly is also fine,
but prevents recycling internal buffers:
```python
with multiprocessing.pool.ThreadPool() as pool:
    pool.map(lambda s: Pyrodigal(meta=True).find_genes(s), sequences)
```

## 💡 Example

Using [Biopython](https://biopython.org/), load a sequence from a GenBank file,
use Prodigal to find all genes it contains, and print the proteins in FASTA
format:
```python
record = Bio.SeqIO.read("sequence.fa", "genbank")
p = pyrodigal.Pyrodigal(meta=True)

for i, gene in enumerate(p.find_genes(str(record.seq))):
    print(f"> {record.id}_{i+1}")
    print(textwrap.fill(record.translate()))
```

## 📜 License

This library, like the original Prodigal software, is provided under the
[GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/).
