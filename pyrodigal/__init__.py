from collections.abc import (
    Sequence as _Sequence,
    Sized as _Sized,
)

from . import _pyrodigal
from ._pyrodigal import (
    Gene,
    Genes,
    Node,
    Nodes,
    Prediction,
    Predictions,
    Pyrodigal,
    Sequence,
    TrainingInfo,
    MetagenomicBin,
)

__doc__ = _pyrodigal.__doc__
__all__ = [
    "Gene",
    "Genes",
    "Node",
    "Nodes",
    "Prediction",
    "Predictions",
    "Pyrodigal",
    "Sequence",
    "TrainingInfo",
    "MetagenomicBin"
]

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"
__version__ = "0.6.2"

_Sized.register(Sequence)
_Sequence.register(Genes)
_Sequence.register(Nodes)
_Sequence.register(Predictions)

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version of the
    library on `Read The Docs <https://pyrodigal.readthedocs.io/en/v{}/>`_.

    """.format(__version__)
