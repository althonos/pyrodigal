from collections.abc import (
    Sequence as _Sequence,
    Sized as _Sized,
)

from ._version import __version__

from . import _pyrodigal
from ._pyrodigal import (
    Gene,
    Genes,
    Mask,
    Masks,
    Node,
    Nodes,
    OrfFinder,
    Sequence,
    TrainingInfo,
    MetagenomicBin,
)

__doc__ = _pyrodigal.__doc__
__all__ = [
    "Gene",
    "Genes",
    "Mask",
    "Masks",
    "Node",
    "Nodes",
    "OrfFinder",
    "Sequence",
    "TrainingInfo",
    "MetagenomicBin"
]

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"

_Sized.register(Sequence)
_Sequence.register(Genes)
_Sequence.register(Masks)
_Sequence.register(Nodes)

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version
    of the library on
    `Read The Docs <https://pyrodigal.readthedocs.io/en/v{}/>`_.

    """.format(__version__)
