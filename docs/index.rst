Pyrodigal |Stars|
=================

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pyrodigal.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pyrodigal/stargazers

*Cython bindings and Python interface to* `Prodigal <https://github.com/hyattpd/Prodigal/>`_,
*an ORF finder for genomes and metagenomes*. **Now with SIMD!**

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads| |Paper|

.. |Actions| image:: https://img.shields.io/github/workflow/status/althonos/pyrodigal/Test/main?logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pyrodigal/actions

.. |GitLabCI| image:: https://img.shields.io/gitlab/pipeline/larralde/pyrodigal/main?gitlab_url=https%3A%2F%2Fgit.embl.de&logo=gitlab&style=flat-square&maxAge=600
   :target: https://git.embl.de/larralde/pyrodigal/-/pipelines

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pyrodigal?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pyrodigal/

.. |PyPI| image:: https://img.shields.io/pypi/v/pyrodigal.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pyrodigal

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pyrodigal?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pyrodigal

.. |AUR| image:: https://img.shields.io/aur/version/python-pyrodigal?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pyrodigal

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pyrodigal?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyrodigal/#files

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pyrodigal.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyrodigal/#files

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pyrodigal.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pyrodigal/#files

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/mit/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyrodigal/

.. |Mirror| image:: https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400
   :target: https://git.embl.de/larralde/pyrodigal/

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pyrodigal.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pyrodigal/issues

.. |Docs| image:: https://img.shields.io/readthedocs/pyrodigal?style=flat-square&maxAge=3600
   :target: http://pyrodigal.readthedocs.io/en/stable/?badge=stable

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyrodigal/blob/main/CHANGELOG.md

.. |Downloads| image:: https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fpyrodigal
   :target: https://pepy.tech/project/pyrodigal

.. |Paper| image:: https://img.shields.io/badge/paper-JOSS-9400ff?style=flat-square&maxAge=86400
   :target: https://doi.org/10.21105/joss.04296


Overview
--------

Pyrodigal is a Python module that provides bindings to Prodigal using
`Cython <https://cython.org/>`_. It directly interacts with the Prodigal
internals, which has the following advantages:

- **single dependency**: Pyrodigal is distributed as a Python package, so you
  can add it as a dependency to your project, and stop worrying about the
  Prodigal binary being present on the end-user machine.
- **no intermediate files**: Everything happens in memory, in a Python object
  you fully control, so you don't have to invoke the Prodigal CLI using a
  sub-process and temporary files. Sequences can be passed directly as
  strings or bytes, which avoids the overhead of formatting your input to
  FASTA for Prodigal.
- **lower memory usage**: Pyrodigal is slightly more conservative when it comes
  to using memory, which can help process very large sequences. It also lets
  you save some more memory when running several *meta*-mode analyses
- **better performance**: Pyrodigal uses *SIMD* instructions to compute which
  dynamic programming nodes can be ignored when scoring connections. This can
  save from a third to half the runtime depending on the sequence.

Setup
-----

Run ``pip install pyrodigal`` in a shell to download the latest release and all
its dependencies from PyPi, or have a look at the
:doc:`Installation page <install>` to find other ways to install ``pyrodigal``.


Library
-------

.. toctree::
   :maxdepth: 2

   Installation <install>
   Contributing <contributing>
   Publications <publications>
   API Reference <api/index>
   Changelog <changes>


License
-------

This library is provided under the `GNU General Public License v3.0 <https://choosealicense.com/licenses/gpl-3.0/>`_.
The Prodigal code was written by `Doug Hyatt <https://github.com/hyattpd>`_ and is distributed under the
terms of the GPLv3 as well. See ``vendor/Prodigal/LICENSE`` for more information. The ``cpu_features`` library was
written by `Guillaume Chatelet <https://github.com/gchatelet>`_ and is licensed under the terms of the
`Apache License 2.0 <https://choosealicense.com/licenses/apache-2.0/>`_. See ``vendor/cpu_features/LICENSE``
for more information.

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the original* `Prodigal`_ *authors. It was developed by* `Martin Larralde <https://github.com/althonos>`_ *during his
PhD project at the* `European Molecular Biology Laboratory <https://www.embl.de/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.
