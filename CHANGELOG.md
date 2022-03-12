# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyrodigal/compare/v0.7.0...HEAD


## [v0.7.0] - 2022-03-12
[v0.7.0]: https://github.com/althonos/pyrodigal/compare/v0.6.4...v0.7.0

### Added
- Support for setting a custom minimum gene length in `pyrodigal.OrfFinder`.
- `Genes.write_scores` method to write the node scores to a file.
- `Gene.__repr__` and `Node.__repr__` methods to display some useful attributes.
- `Sequence.__str__` method to get back a nucleotide string from a `Sequence` object.

### Changed
- Use a more compact data structure to store `Gene` data.

### Fixed
- `Nodes._calc_orf_gc` reading nucleotides after the sequence end when computing GC content for edge nodes.

### Removed
- `pyrodigal.Pyrodigal` class (use `pyrodigal.OrfFinder` instead).
- `pyrodigal.Predictions` class (functionality merged into `pyrodigal.Genes`).


## [v0.6.4] - 2021-12-23
[v0.6.4]: https://github.com/althonos/pyrodigal/compare/v0.6.3...v0.6.4

### Added
- `load` and `dump` methods to `TrainingInfo` for storing and loading a raw training info structure.
- Support for creating an `OrfFinder` pre-configured with a training info.
- `-t` and `-n` flags to the CLI.


## [v0.6.3] - 2021-12-23
[v0.6.3]: https://github.com/althonos/pyrodigal/compare/v0.6.2...v0.6.3

### Added
- `pyrodigal` command line script exposing a CLI mimicking the original `prodigal` binary.
- `write_gff`, `write_genes` and `write_translations` methods to `pyrodigal.Predictions` to write the predictions results to a file in different formats.
- Implementation for masking regions of unknown nucleotides in input sequences.

### Changed
- Renamed `pyrodigal.Pyrodigal` class to `pyrodigal.OrfFinder`.

### Fixed
- `setup.py` build different SIMD implementations with the same set of feature flags, causing compilers to re-optimize the SIMD implementations.


## [v0.6.2] - 2021-09-25
[v0.6.2]: https://github.com/althonos/pyrodigal/compare/v0.6.1...v0.6.2

### Added
- Sphinx documentation with small install guide and API reference.

### Fixed
- `setup.py` not detecting SSE2 and AVX2 build support because of a linker error.

### Changed
- Build OSX extension without AVX2 support since runtime detection of AVX2 to avoid the `Illegal Instruction: 4` bug on older CPUs.


## [v0.6.1] - 2021-09-24
[v0.6.1]: https://github.com/althonos/pyrodigal/compare/v0.6.0...v0.6.1

### Fixed
- Source distribution lacking C files necessary for building `cpu_features`.


## [v0.6.0] - 2021-09-23
[v0.6.0]: https://github.com/althonos/pyrodigal/compare/v0.5.4...v0.6.0

### Added
- SIMD code to build an index of which connections can be skipped when scoring node connections in the dynamic programming routine ([#6](https://github.com/althonos/pyrodigal/pull/6)).


## [v0.5.4] - 2021-09-18
[v0.5.4]: https://github.com/althonos/pyrodigal/compare/v0.5.3...v0.5.4

### Added
- `Prediction.confidence` method to compute the confidence for a prediction like reported in Prodigal's GFF output.
- `Prediction.sequence` method get the nucleotide sequence of a predicted gene ([#4](https://github.com/althonos/pyrodigal/issues/4)).

### Changed
- Replaced internal storage of input sequences to use a byte array instead of a bitmap.

### Fixed
- Extract `Prediction.gc_cont` number directly from the start node instead of the text representation to get full accuracy.
- Prodigal bug causing nodes on the reverse strand to always receive a penalty instead of penalizing only small ORFs ([hyattpd/Prodigal#88](https://github.com/hyattpd/Prodigal/pull/88)).


## [v0.5.3] - 2021-09-12
[v0.5.3]: https://github.com/althonos/pyrodigal/compare/v0.5.2...v0.5.3

### Fixed
- `Prediction.translate` not translating the last unknown codon properly for genes on the direct strand.


## [v0.5.2] - 2021-09-11
[v0.5.2]: https://github.com/althonos/pyrodigal/compare/v0.5.1...v0.5.2

### Changed
- Make `Pyrodigal.train` return a reference to the newly created `TrainingInfo` for inspection if needed.
- Reimplement `add_nodes` and `add_genes` to use a growable array instead of counting and pre-allocating the C arrays.

### Fixed
- Inconsistent handling of unknown nucleotides in input sequences and gene translations.


## [v0.5.1] - 2021-09-04
[v0.5.1]: https://github.com/althonos/pyrodigal/compare/v0.5.0...v0.5.1

### Added
- Additional `Gene` properties to access the score

### Changed
- Use more efficient `PyUnicode` macros when reading or creating a string containing a nucleotide or a protein sequence.
- Release the GIL when creating a bitmap for an `str` given as input to `Pyrodigal.find_genes`.
- Release the GIL when creating the protein sequence returned by `Gene.translate`.

### Fixed
- `Pyrodigal.find_genes` and `Gene.translate` not behaving like Prodigal when handling sequences with unknown nucleotides.


## [v0.5.0] - 2021-06-15
[v0.5.0]: https://github.com/althonos/pyrodigal/compare/v0.4.7...v0.5.0

### Added
- `pyrodigal.TrainingInfo` class exposing variables obtained during training as an attribute to `Pyrodigal`, `Gene` and `Genes` instance.
- Support for passing objects implementing the buffer protocol to `Pyrodigal.find_genes` and `Pyrodigal.train` instead of requiring `str` sequences.

### Fixed
- Potential data race on training info in case a `Gene.translate` with a non-default translation table was being translated at the same time as a `Pyrodigal.find_genes` call.
- Spurious handling of Unicode strings causing potential issues on platform using a different base encoding.


## [v0.4.7] - 2021-04-09
[v0.4.7]: https://github.com/althonos/pyrodigal/compare/v0.4.6...v0.4.7

### Fixed
- `Pyrodigal.find_genes` segfaulting on some sequences when called in `single` mode ([#2](https://github.com/althonos/pyrodigal/issues/2)).
- `MemoryError` potentially not being properly raised on allocation issues for sequence bitmaps.


## [v0.4.6] - 2021-03-05
[v0.4.6]: https://github.com/althonos/pyrodigal/compare/v0.4.5...v0.4.6

### Changed
- Tests are now in the `pyrodigal.tests` module and can be run after a site install.

### Fixed
- `Pyrodigal.find_genes` stalling on sequences shorter than 3 nucleotides.


## [v0.4.5] - 2021-03-03
[v0.4.5]: https://github.com/althonos/pyrodigal/compare/v0.4.4...v0.4.5

### Fixed
- Compilation of OSX and Windows wheels.


## [v0.4.4] - 2021-03-03
[v0.4.4]: https://github.com/althonos/pyrodigal/compare/v0.4.3...v0.4.4

### Fixed
- Mark package as OS-independent.

### Added
- Support for Python 3.5.
- Compilation of PyPy wheels on OSX.


## [v0.4.3] - 2021-03-01
[v0.4.3]: https://github.com/althonos/pyrodigal/compare/v0.4.2...v0.4.3

### Fixed
- Buffer overflow when running in `meta` mode on a sequence too small to have any dynamic programming nodes.


## [v0.4.2] - 2021-02-07
[v0.4.2]: https://github.com/althonos/pyrodigal/compare/v0.4.1...v0.4.2

### Fixed
- Buffer overflow coming from the node array, caused by an incorrect
  estimation of the node count from the sequence length.


## [v0.4.1] - 2021-01-07
[v0.4.1]: https://github.com/althonos/pyrodigal/compare/v0.4.0...v0.4.1

### Removed
- Python 3.5 from the project metadata (the code was only compatible with
  Python 3.6+ already because of *f-strings*).

### Fixed
- Broken linking of static `libprodigal` against the `_pyrodigal` extension
  on some OSX environments ([bioconda/bioconda-recipes#25568](https://github.com/bioconda/bioconda-recipes/pull/25568)).


## [v0.4.0] - 2021-01-06
[v0.4.0]: https://github.com/althonos/pyrodigal/compare/v0.3.2...v0.4.0

### Changed
- `trans_table` keyword argument to `Pyrodigal.train` has been renamed
  to `translation_table`.

### Added
- Option to change the translation table to any allowed number in `Gene.translate`
  ([#1](https://github.com/althonos/pyrodigal/issues/1)).


## [v0.3.2] - 2020-11-27
[v0.3.2]: https://github.com/althonos/pyrodigal/compare/v0.3.1...v0.3.2

### Fixed
- Broken compilation of PyPy wheels in Travis-CI.


## [v0.3.1] - 2020-11-27

[v0.3.1]: https://github.com/althonos/pyrodigal/compare/v0.3.0...v0.3.1

### Added
- Link to Zenodo record in `README.md`.
- `Typing :: Typed` classifier to the PyPI metadata.
- Explicit support for Python 3.9.

### Changed
- Streamlined compilation process when building from source distribution.



## [v0.3.0] - 2020-09-07

[v0.3.0]: https://github.com/althonos/pyrodigal/compare/v0.2.4...v0.3.0

### Added
- Thread-safety for all `Pyrodigal` methods

### Fixed
- Reduced total amount of memory used to allocated dynamic programming
  nodes for a given sequence.


## [v0.2.4] - 2020-09-04

[v0.2.4]: https://github.com/althonos/pyrodigal/compare/v0.2.3...v0.2.4

### Added
- Precompiled wheels for Windows x86-64 platform.

### Changed
- Compilation of large `Prodigal/training.c` file is now done in chunks
  and uses `static const` to reduce build time.


## [v0.2.3] - 2020-08-09

[v0.2.3]: https://github.com/althonos/pyrodigal/compare/v0.2.2...v0.2.3

### Fixed
- Buffer overflow issue with Pyrodigal in `closed=False` mode.


## [v0.2.2] - 2020-07-14

[v0.2.2]: https://github.com/althonos/pyrodigal/compare/v0.2.0...v0.2.2

### Added
- Access to the translation table of a `Gene` object.


## [v0.2.1] - 2020-05-29

[v0.2.1]: https://github.com/althonos/pyrodigal/compare/v0.2.0...v0.2.1

### Fixed
- Memory issues causing PyPy to crash when using `Pyrodigal` in single mode.


## [v0.2.0] - 2020-05-28

[v0.2.0]: https://github.com/althonos/pyrodigal/compare/v0.1.1...v0.2.0

### Added
- Support for Prodigal's *single* mode.


## [v0.1.1] - 2020-04-30

[v0.1.1]: https://github.com/althonos/pyrodigal/compare/v0.1.0...v0.1.1

### Added
- Distribution of CPython wheels for ManyLinux2010 and OSX platforms.


## [v0.1.0] - 2020-04-27
[v0.1.0]: https://github.com/althonos/pyrodigal/compare/0a90bf9...v0.1.0

Initial release.
