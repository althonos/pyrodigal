# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyrodigal/compare/v0.4.6...HEAD


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
