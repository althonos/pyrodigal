from collections.abc import (
    Sequence as _Sequence,
    Sized as _Sized,
)

from ._version import __version__

from . import lib
from .lib import (
    Gene,
    Genes,
    GeneFinder,
    Mask,
    Masks,
    Node,
    Nodes,
    Sequence,
    TrainingInfo,
    MetagenomicBin,
    MetagenomicBins,
    MIN_SINGLE_GENOME,
    IDEAL_SINGLE_GENOME,
    METAGENOMIC_BINS,
    PRODIGAL_VERSION,
    TRANSLATION_TABLES,
)

__doc__ = lib.__doc__
__all__ = [
    "Gene",
    "Genes",
    "GeneFinder",
    "Mask",
    "Masks",
    "Node",
    "Nodes",
    "Sequence",
    "TrainingInfo",
    "MetagenomicBin",
    "MetagenomicBins",
    "IDEAL_SINGLE_GENOME",
    "METAGENOMIC_BINS",
    "MIN_SINGLE_GENOME",
    "PRODIGAL_VERSION",
    "TRANSLATION_TABLES",
]

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"

_Sized.register(Sequence)
_Sequence.register(Genes)
_Sequence.register(Masks)
_Sequence.register(Nodes)
_Sequence.register(MetagenomicBins)

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version
    of the library on
    `Read The Docs <https://pyrodigal.readthedocs.io/en/v{}/>`_.

    """.format(
        __version__
    )
