# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyrodigal/compare/v0.3.0...HEAD


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
