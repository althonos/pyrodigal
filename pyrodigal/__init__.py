import pkg_resources as _pkg_resources
from collections.abc import Sequence as _Sequence

from . import _pyrodigal
from ._pyrodigal import Gene, Genes, Pyrodigal

__all__ = ["Gene", "Genes", "Pyrodigal"]
__doc__ = _pyrodigal.__doc__

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"
__version__ = (
    _pkg_resources
    .resource_string(__name__, "_version.txt")
    .decode("utf-8")
    .strip()
)

_Sequence.register(Genes)
