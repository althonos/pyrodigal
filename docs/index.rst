|Logo| Pyrodigal |Stars|
========================

.. |Logo| image:: /_images/logo.png
   :scale: 40%
   :class: dark-light

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pyrodigal.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pyrodigal/stargazers
   :class: dark-light

*Cython bindings and Python interface to* `Prodigal <https://github.com/hyattpd/Prodigal/>`_,
*an ORF finder for genomes and metagenomes*. **Now with SIMD!**

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads| |Paper| |Citations|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/pyrodigal/test.yml?branch=main&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pyrodigal/actions
   :class: dark-light

.. |GitLabCI| image:: https://img.shields.io/gitlab/pipeline/larralde/pyrodigal/main?gitlab_url=https%3A%2F%2Fgit.embl.de&logo=gitlab&style=flat-square&maxAge=600
   :target: https://git.embl.de/larralde/pyrodigal/-/pipelines
   :class: dark-light

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pyrodigal?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pyrodigal/
   :class: dark-light

.. |PyPI| image:: https://img.shields.io/pypi/v/pyrodigal.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pyrodigal
   :class: dark-light

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pyrodigal?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pyrodigal
   :class: dark-light

.. |AUR| image:: https://img.shields.io/aur/version/python-pyrodigal?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pyrodigal
   :class: dark-light

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pyrodigal?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyrodigal/#files
   :class: dark-light

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pyrodigal.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyrodigal/#files
   :class: dark-light

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pyrodigal.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pyrodigal/#files
   :class: dark-light

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/mit/
   :class: dark-light

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=3600&style=flat-square
   :target: https://github.com/althonos/pyrodigal/
   :class: dark-light

.. |Mirror| image:: https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=3600
   :target: https://git.embl.de/larralde/pyrodigal/
   :class: dark-light

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pyrodigal.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pyrodigal/issues
   :class: dark-light

.. |Docs| image:: https://img.shields.io/readthedocs/pyrodigal?style=flat-square&maxAge=3600
   :target: http://pyrodigal.readthedocs.io/en/stable/?badge=stable
   :class: dark-light

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=3600&style=flat-square
   :target: https://github.com/althonos/pyrodigal/blob/main/CHANGELOG.md
   :class: dark-light

.. |Downloads| image:: https://img.shields.io/pypi/dm/pyrodigal?style=flat-square&color=303f9f&maxAge=86400&label=downloads
   :target: https://pepy.tech/project/pyrodigal
   :class: dark-light

.. |Paper| image:: https://img.shields.io/badge/paper-JOSS-9400ff?style=flat-square&maxAge=86400
   :target: https://doi.org/10.21105/joss.04296
   :class: dark-light

.. |Citations| image:: https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fbadge.dimensions.ai%2Fdetails%2Fid%2Fpub.1147419140%2Fmetadata.json&query=%24.times_cited&style=flat-square&label=citations&maxAge=86400
   :target: https://badge.dimensions.ai/details/id/pub.1147419140
   :class: dark-light


Overview
--------

Pyrodigal is a Python module that provides bindings to Prodigal using
`Cython <https://cython.org/>`_. It directly interacts with the Prodigal
internals, which has the following advantages:


.. grid:: 1 2 3 3
   :gutter: 1

   .. grid-item-card:: :fas:`battery-full` Batteries-included

      Just add ``pyrodigal`` as a ``pip`` or ``conda`` dependency, no need
      for the Prodigal binary or any external dependency.

   .. grid-item-card:: :fas:`screwdriver-wrench` Flexible I/O

      Directly pass sequences to process as Python `str` objects, no 
      need for intermediate files.

   .. grid-item-card:: :fas:`memory` Memory-efficient

      Benefit from conservative memory allocation and a reworked data layout 
      for candidate nodes.

   .. grid-item-card:: :fas:`microchip` Faster computation

      Use the full power of your CPU with :wiki:`SIMD` instructions to 
      filter out candidate genes prior to the scoring stage.

   .. grid-item-card:: :fas:`check` Consistent results 

      Get the same results as Prodigal ``v2.6.3+31b300a``, with additional
      bug fixes compared to the latest stable Prodigal version.

   .. grid-item-card:: :fas:`toolbox` Feature-complete

      Access all the features of the original CLI through the :doc:`Python API <api/index>` 
      or a :doc:`drop-in CLI replacement <guide/cli>`.


Features
--------

The library now features everything from the original Prodigal CLI:

- **run mode selection**: Choose between *single* mode, using a training
  sequence to count nucleotide hexamers, or *metagenomic* mode, using
  pre-trained data from different organisms (``prodigal -p``).
- **region masking**: Prevent genes from being predicted across regions
  containing unknown nucleotides  (``prodigal -m``).
- **closed ends**: Genes will be identified as running over edges if they
  are larger than a certain size, but this can be disabled (``prodigal -c``).
- **training configuration**: During the training process, a custom
  translation table can be given (``prodigal -g``), and the Shine-Dalgarno motif
  search can be forcefully bypassed (``prodigal -n``)
- **output files**: Output files can be written in a format mostly
  compatible with the Prodigal binary, including the protein translations
  in FASTA format (``prodigal -a``), the gene sequences in FASTA format
  (``prodigal -d``), or the potential gene scores in tabular format
  (``prodigal -s``). See the :doc:`Output Formats <guide/outputs>` section 
  for supported formats.
- **training data persistence**: Getting training data from a sequence and
  using it for other sequences is supported; in addition, a training data
  file can be saved and loaded transparently (``prodigal -t``).

In addition, the **new** features are available:

- **custom gene size threshold**: While Prodigal uses a minimum gene size
  of 90 nucleotides (60 if on edge), Pyrodigal allows to customize this
  threshold, allowing for smaller ORFs to be identified if needed.

Several changes were done regarding **memory management**:

- **digitized sequences**: Sequences are stored as raw bytes instead of compressed 
  bitmaps. This means that the sequence itself takes 3/8th more space, but since 
  the memory used for storing the sequence is often negligible compared to the 
  memory used to store dynamic programming nodes, this is an acceptable 
  trade-off for better performance when extracting said nodes.
- **node buffer growth**: Node arrays are dynamically allocated and grow 
  exponentially instead of being pre-allocated with a large size. On small 
  sequences, this leads to Pyrodigal using about 30% less memory.
- **lightweight genes**: Genes are stored in a more compact data structure than in 
  Prodigal (which reserves a buffer to store string data), saving around 1KiB 
  per gene.


Setup
-----

Run ``pip install pyrodigal`` in a shell to download the latest release and all
its dependencies from PyPi, or have a look at the
:doc:`Installation page <guide/install>` to find other ways to install ``pyrodigal``.


Citation
--------

Pyrodigal is scientific software, with a
`published paper <https://doi.org/10.21105/joss.04296>`_
in the `Journal of Open-Source Software <https://joss.theoj.org/>`_. Check the 
:doc:`Publications page <guide/publications>` to see how to cite Pyrodigal properly.


Library
-------

Check the following pages of the user guide or the API reference for more
in-depth reference about library setup, usage, and rationale:

.. toctree::
   :maxdepth: 2

   User Guide <guide/index>
   API Reference <api/index>


Related Projects
----------------

The following Python libraries may be of interest for bioinformaticians.

.. include:: related.rst


License
-------

This library is provided under the `GNU General Public License v3.0 <https://choosealicense.com/licenses/gpl-3.0/>`_.
The Prodigal code was written by `Doug Hyatt <https://github.com/hyattpd>`_ and is distributed under the
terms of the GPLv3 as well. See the :doc:`Copyright Notice <guide/copyright>` section
for the full GPLv3 license.

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the original* `Prodigal`_ *authors. It was developed by* `Martin Larralde <https://github.com/althonos>`_ *during his
PhD project at the* `European Molecular Biology Laboratory <https://www.embl.de/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.

The project icon was derived from `UXWing <https://uxwing.com/>`_ and is re-used
under `their permissive license <https://uxwing.com/license/>`_.
