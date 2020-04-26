# 🔥 Pyrodigal

*Python interface to Prodigal, an ORF finder for genomes, progenomes and metagenomes.*

[![TravisCI](https://img.shields.io/travis/althonos/pyrodigal/master.svg?logo=travis&maxAge=600&style=flat-square)](https://travis-ci.com/althonos/pyrodigal/branches)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyrodigal?style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/pyrodigal/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/pyrodigal.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyrodigal/issues)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyrodigal/)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyrodigal.py/blob/master/CHANGELOG.md)

<!-- [![AppVeyor](https://img.shields.io/appveyor/ci/althonos/pyrodigal/master?logo=appveyor&style=flat-square&maxAge=600)](https://ci.appveyor.com/project/althonos/pyrodigal) -->
<!-- [![PyPI](https://img.shields.io/pypi/v/pyrodigal.svg?style=flat-square&maxAge=600)](https://pypi.org/project/pyrodigal) -->
<!-- [![Wheel](https://img.shields.io/pypi/wheel/pyrodigal.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyrodigal/#files) -->
<!-- [![Python Versions](https://img.shields.io/pypi/pyversions/pyrodigal.svg?style=flat-square&maxAge=600)](https://pypi.org/project/pyrodigal/#files) -->
<!-- [![Python Implementations](https://img.shields.io/pypi/implementation/pyrodigal.svg?style=flat-square&maxAge=600)](https://pypi.org/project/pyrodigal/#files) -->

## 💡 Example

Using [Biopython](https://biopython.org/), load a sequence from a FASTA file,
use Prodigal to find all genes it contains, and create a list of `SeqRecord`
containing each of the proteins:
```python
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.SeqIO import read
from Bio.SeqRecord import SeqRecord
from pyrodigal import Pyrodigal

record = read("sequence.fa", "fasta")
proteins = []

p = Pyrodigal(meta=True)
for i, gene in enumerate(p.find_genes(str(record.seq))):
    protein = SeqRecord(
        Seq(gene.translate(), alphabet=generic_protein),
        id=f"{record.id}_{i+1}",
        name=f"{record.name}_{i+1}",
        description=f"{record.description} # {gene.start} # {gene.stop} # {gene.strand}"
    )
    proteins.append(protein)
```

## 📜 License

This library, like the original Prodigal software, is provided under the
[GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/).
