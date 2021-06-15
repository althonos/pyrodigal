from collections.abc import Sequence as _Sequence

from . import _pyrodigal
from ._pyrodigal import Gene, Genes, Pyrodigal, TrainingInfo

__all__ = ["Gene", "Genes", "Pyrodigal", "TrainingInfo"]
__doc__ = _pyrodigal.__doc__

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"
__version__ = "0.4.7"

_Sequence.register(Genes)
