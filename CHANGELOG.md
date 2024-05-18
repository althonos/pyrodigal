# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyrodigal/compare/v3.4.0...HEAD


## [v3.4.0] - 2024-05-19
[v3.4.0]: https://github.com/althonos/pyrodigal/compare/v3.3.0...v3.4.0

### Added
- `strict` argument to `Gene.translate` to control translation of ambiguous codons with unambiguous translation ([#54](https://github.com/althonos/pyrodigal/issues/54)).
- `strict_translation` argument to `Genes.write_genbank` and `Genes.write_translation`.
- Support for translation tables 26 to 33 in `Gene.translate`.
- Support for translation tables 26, 29, 30, 32 and 33 in `GeneFinder.train`.
- `Genes.score` property to count the total score of all extracted genes.
- `full_id` parameter to `Genes.write_gff`, `Genes.write_translation` and `Genes.write_genes` to control the `ID` field written for each gene ([#53](https://github.com/althonos/pyrodigal/issues/53)).

### Changed
- `Gene.translate` now raises a warning when called with a translation table incompatible with the training info.

### Fixed
- Bug in code for masking trailing nucleotides ([#55](https://github.com/althonos/pyrodigal/issues/55)).


## [v3.3.0] - 2024-01-24
[v3.3.0]: https://github.com/althonos/pyrodigal/compare/v3.2.2...v3.3.0

### Added
- CLI option to disable translation of stop codons ([#51](https://github.com/althonos/pyrodigal/pull/51), by [@zclaas](https://github.com/zclaas)).

### Changed
- `Scorer` internal API to separate connection scoring and overlap disentangling.

### Fixed
- Bug with computation of minimum node in connection scoring loop ([hyattpd/Prodigal#108](https://github.com/hyattpd/Prodigal/pull/108)).
- Out-of-bounds sequence access in `_shine_dalgarno_exact` and `_shine_dalgarno_mm` methods of `Sequence`.
- Memory leak in `Nodes.__setstate__` caused by incorrect reallocation.


## [v3.2.2] - 2024-01-21
[v3.2.2]: https://github.com/althonos/pyrodigal/compare/v3.2.1...v3.2.2

### Fixed
- Always mark SSE2 support on x86-64 CPUs independently of `archspec`-detected features ([#49](https://github.com/althonos/pyrodigal/issues/49)).


## [v3.2.1] - 2023-11-27
[v3.2.1]: https://github.com/althonos/pyrodigal/compare/v3.2.0...v3.2.1

### Added
- Option to change argument parser in `pyrodigal.cli.main`.


## [v3.2.0] - 2023-11-27
[v3.2.0]: https://github.com/althonos/pyrodigal/compare/v3.1.1...v3.2.0

### Added
- AVX-512 implementation of the SIMD pre-filter.
- Additional support for reading `lz4` and `xz` and `zstd`-compressed input in the CLI.
- Option to change gene finder type in `pyrodigal.cli.main`.


## [v3.1.1] - 2023-11-06
[v3.1.1]: https://github.com/althonos/pyrodigal/compare/v3.1.0...v3.1.1

### Fixed
- Incorrect unpickling of `GeneFinder` causing crashes with multiprocessing ([#46](https://github.com/althonos/pyrodigal/issues/46)).


## [v3.1.0] - 2023-10-22
[v3.1.0]: https://github.com/althonos/pyrodigal/compare/v3.0.1...v3.1.0

### Added
- Support for Python 3.12.
- `min_mask` argument to `GeneFinder` to control the minimum lenght of masked regions on `mask=True`.


## [v3.0.1] - 2023-09-27
[v3.0.1]: https://github.com/althonos/pyrodigal/compare/v3.0.0...v3.0.1

### Fixed
- `Genes.write_scores` and `Genes.write_gff` crashing on empty `Genes` ([#44](https://github.com/althonos/pyrodigal/issues/44)).


## [v3.0.0] - 2023-09-17
[v3.0.0]: https://github.com/althonos/pyrodigal/compare/v2.3.0...v3.0.0

### Added
- `MetagenomicBins` collection to store a dense array of `MetagenomicBin` objects.
- `metagenomic_bins` keyword argument to `GeneFinder` allowing to control which models are used when running gene finding in *meta* mode ([#24](https://github.com/althonos/pyrodigal/issues/24)).
- `metagenomic_bin` attribute to `Genes` referencing the metagenomic model with which the genes were predicted, if in *meta* mode.
- Additional `TrainingInfo` properties (`missing_motif_weight`, `coding_statistics`).
- Setters for all remaining `TrainingInfo` properties.
- Proper `TrainingInfo` constructor with configuration option for all attributes.
- `TrainingInfo.to_dict` method to extract all parameters from a `TrainingInfo`.
- `Genes.write_genbank` method to write a GenBank record with all predicted genes from a sequence.
- `include_stop` flag to `Gene.translate` and `Genes.write_translations` to allow excluding the stop codon from the translated sequence.
- `include_translation_table` flag to `Genes.write_gff` to include the translation table to the GFF attributes of each gene.
- `gbk` output format to the Pyrodigal CLI.
- `Sequence.unknown` property exposing the number of unknown nucleotides in the sequence.
- `Sequence.start_probability` and `Sequence.stop_probability` to estimate the probability of encountering a start and a stop codon based on the GC%.

### Fixed
- `Genes.write_gff` not properly reporting the number of bytes written.
- Merge several `nogil` sections in `Sequence` constructor.
- Several Cython functions missing a `noexcept` qualifier.

### Changed
- **BREAKING**: Rename `OrfFinder` to `GeneFinder` for consistency.
- **BREAKING**: Use `memoryview` to expose all `TrainingInfo` attributes instead manually building lists or tuples.
- Reorganize memory management of the built-in metagenomic models.
- Make the internal Cython model public (`pyrodigal.lib`) to allow importing the underlying classes in other Cython projects.
- Use `typing.Literal` for allowed translation table values in `pyrodigal.lib` annotations
- Cache intermediate log-odds in `Nodes._raw_coding_score` to reduce calls to `pow` and `log` functions.
- Inline connection scoring functions to reduce function call overhead.
- Reorganize `struct _node` fields to reduce size in memory.
- Make `GeneFinder.find_genes` and `GeneFinder.train` reserve memory for the `Nodes` based on the GC% of the input sequence.
- Avoid storing temporary results in the generic implementation of `ConnectionScorer.compute_skippable`.
- Use Cython `freelist` for allocating `Node`, `Gene`, `MetagenomicBin` and `Mask`.
- Increase minimum allocation for `Genes` and `Nodes` to reduce early reallocations.

### Removed
- **BREAKING**: `metagenomic_bin` attribute of `TrainingInfo`.


## [v2.3.0] - 2023-07-20
[v2.3.0]: https://github.com/althonos/pyrodigal/compare/v2.2.0...v2.3.0

### Changed
- Bump Cython to `v3.0.0`.


## [v2.2.0] - 2023-06-19
[v2.2.0]: https://github.com/althonos/pyrodigal/compare/v2.1.0...v2.2.0

### Changed
- Release GIL while masking sequence regions in `Sequence.__init__`.
- Use [`archspec`](https://pypi.org/project/archspec) instead of `cpu_features` for runtime feature detection.

### Added
- Support for reading `gzip` and `bz2`-compressed input in the CLI.
- CLI flag to run ORF detection in parallel when input contains several contigs.

### Removed
- Support for Python 3.5.


## [v2.1.0] - 2023-02-20
[v2.1.0]: https://github.com/althonos/pyrodigal/compare/v2.0.4...v2.1.0

### Changed
- Update Prodigal to `v2.6.3+c1e2d36` to fix a bug with Shine-Dalgarno detection on reverse contig edge ([hyattpd/Prodigal#100](https://github.com/hyattpd/Prodigal/pull/100)).

### Added
- CLI flags to set the minimum gene size ([#32](https://github.com/althonos/pyrodigal/pull/32), by [@cjprybol](https://github.com/cjprybol)).

### Fixed
- ArchLinux User Repository package generation in CI.


## [v2.0.4] - 2023-01-09
[v2.0.4]: https://github.com/althonos/pyrodigal/compare/v2.0.3...v2.0.4

### Fixed
- GC% computation and RBS scoring for reverse strand nodes close to the contig edge ([#27](https://github.com/althonos/pyrodigal/issues/27)).


## [v2.0.3] - 2022-12-20
[v2.0.3]: https://github.com/althonos/pyrodigal/compare/v2.0.2...v2.0.3

### Fixed
- `OrfFinder(mask=True)` ignoring the minimum mask size when masking regions ([#26](https://github.com/althonos/pyrodigal/issues/26)).

### Changed
- Use `cibuildhweel` for building wheel distributions.

### Added
- Wheels for MacOS Aarch64 platforms.


## [v2.0.2] - 2022-11-01
[v2.0.2]: https://github.com/althonos/pyrodigal/compare/v2.0.1...v2.0.2

### Fixed
- Syntax issue in Cython files failing build on Bioconda runner.


## [v2.0.1] - 2022-11-01
[v2.0.1]: https://github.com/althonos/pyrodigal/compare/v2.0.0...v2.0.1

### Fixed
- Syntax issue in Cython files failing build on some environments.


## [v2.0.0] - 2022-11-01
[v2.0.0]: https://github.com/althonos/pyrodigal/compare/v1.1.2...v2.0.0

### Added
- MMX implementation of the SIMD prefilter.
- Proper GFF headers and metadata section to GFF output.
- `Sequence.gc_frame_plot` method to compute the max GC frame profile from Python.
- `metagenomic_bin` property to `TrainingInfo` to support recovering the object corresponding to a pre-trained model.
- `meta` attribute to `Genes` to store whether genes were predicted in single or in meta mode.
- `pyrodigal.PRODIGAL_VERSION` constant storing the wrapped Prodigal version.
- `pyrodigal.MIN_SINGLE_GENOME` and `pyrodigal.IDEAL_SINGLE_GENOME` constants storing the minimum and recommended sequence sizes for training.

### Changed
- Make all write methods of `Genes` objects require a ``sequence_id`` argument instead of using the internal sequence number.
- Rewrite SIMD prefilter using a generic template with C macros.
- Make `Mask` record coordinates in start-inclusive end-exclusive mode to follow Python conventions.
- Make connection scoring tests only score some randomly selected node pairs for faster runs.
- Rewrite tests to use `importlib.resources` for managing test data.

### Removed
- `from_bytes` and `from_string` constructors of `Sequence` objects.

### Fixed
- Duplicate extraction of start codons located on contig edges inside `Nodes._extract` ([#21](https://github.com/althonos/pyrodigal/issues/21)).
- Pickling and unpickling of `TrainingInfo` objects corresponding to pre-trained models.
- Implementation of `calc_most_gc_frame` being inconsistent with the Prodigal implementation.
- Implementation of the maximum search in `score_connection_forward_start` not following the (weird?) behaviour from Prodigal ([#21](https://github.com/althonos/pyrodigal/issues/21)).
- Gene identifier being used instead of the sequence identifier in the GFF output ([#18](https://github.com/althonos/pyrodigal/issues/18)).
- Out of bound access to sequence data in `Sequence._shine_dalgarno_mm` and `Sequence._shine_dalgarno_exact`.


## [v1.1.2] - 2022-08-31
[v1.1.2]: https://github.com/althonos/pyrodigal/compare/v1.1.1...v1.1.2

### Changed
- Use the `vbicq` Arm intrinsic in the NEON implementation to combine `vandq` and `vmvnq`.

### Fixed
- Prevent direct instantiation of `Node` and `Gene` objects from Python code.
- Configuration of platform-specific NEON flags in `setup.py` not being applied to the linker.


## [v1.1.1] - 2022-07-08
[v1.1.1]: https://github.com/althonos/pyrodigal/compare/v1.1.0...v1.1.1

### Fixed
- Some `cpu_features` source files not being included in source distribution.


## [v1.1.0] - 2022-06-09
[v1.1.0]: https://github.com/althonos/pyrodigal/compare/v1.0.2...v1.1.0

### Changed
- `OrfFinder.train` can now be given more than one sequence argument to train on contigs from an unclosed genome.
- Updated `cpu_features` to `v0.7.0` and added hardware detection of NEON features on Linux Aarch64 platforms.


## [v1.0.2] - 2022-05-13
[v1.0.2]: https://github.com/althonos/pyrodigal/compare/v1.0.1...v1.0.2

### Fixed
- Detection of Arm64 platform in `setup.py` ([#16](https://github.com/althonos/pyrodigal/issues/16)).


## [v1.0.1] - 2022-04-28
[v1.0.1]: https://github.com/althonos/pyrodigal/compare/v1.0.0...v1.0.1

### Changed
- `pyrodigal.cli` now concatenates training sequences the same way as Prodigal does.


## [v1.0.0] - 2022-04-20
[v1.0.0]: https://github.com/althonos/pyrodigal/compare/v0.7.3...v1.0.0

Stable version, to be published in the [Journal of Open-Source Software](https://joss.theoj.org/).

### Added
- `pickle` protocol implementation for `Nodes`, `TrainingInfo`, `OrfFinder`, `Sequence`, `Masks` and `Genes` objects.
- Buffer protocol implementation for `Sequence`, allowing access to raw digits.
- `__eq__` and `__repr__` magic methods to `Mask` objects.

### Changed
- Optimized code used for region masking to avoid searching for the same mask repeatedly.
- `TRANSLATION_TABLES` and `METAGENOMIC_BINS` are now exposed as constants in the top `pyrodigal` module.
- Refactored connection scoring into different functions based on the type (start/stop) and strand (direct/reverse) of the node being scored.
- Changed the growth factor for dynamic arrays to be the same as the one used in CPython `list` buffers.


## [v0.7.3] - 2022-04-06
[v0.7.3]: https://github.com/althonos/pyrodigal/compare/v0.7.2...v0.7.3

### Added
- `Gene.score` property to get the gene score as reported in the score data string.

### Fixed
- `OrfFinder.find_genes` not producing consistent results across runs in *meta* mode ([#13](https://github.com/althonos/pyrodigal/issues/13)).
- `OrfFinder.find_genes` returning `Nodes` with incomplete score information.


## [v0.7.2] - 2022-03-15
[v0.7.2]: https://github.com/althonos/pyrodigal/compare/v0.7.1...v0.7.2

### Changed
- Improve performance of `mer_ndx` and `score_connection` using dedicated implementations with better branch prediction.
- Mark arguments as `const` in C code where possible.

### Fixed
- Signatures of Cython classes not displaying properly because of the `embedsignature` directive.
- `_sequence.h` functions not being inlined as expected.


## [v0.7.1] - 2022-03-14
[v0.7.1]: https://github.com/althonos/pyrodigal/compare/v0.7.0...v0.7.1

### Changed
- Rewrite internal `Sequence` code using inlined functions to increase performance when the strand is known.

### Fixed
- `Nodes.copy` potentially failing on empty collections after trying to allocate 0 bytes.
- `TestGenes.test_write_scores` failing on some machines because of float rounding issues.
- `Gene.translate` ignoring the `unknown_residue` argument value and always using `"X"`.
- Memory leak in `Pyrodigal.train` cause by memory not being freed after building the GC frame plot.


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
