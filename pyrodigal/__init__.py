from collections.abc import Sequence as _Sequence

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
__version__ = "0.5.2"

_Sequence.register(Genes)
_Sequence.register(Nodes)
_Sequence.register(Predictions)
