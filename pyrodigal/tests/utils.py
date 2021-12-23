import gzip
import os

from .fasta import parse


def load_record(name):
    data = os.path.realpath(os.path.join(__file__, "..", "data"))
    fna = os.path.join(data, "{name}.fna.gz".format(name=name))
    with gzip.open(fna, "rt") as f:
        return next(parse(f))

def load_proteins(name, mode):
    data = os.path.realpath(os.path.join(__file__, "..", "data"))
    faa = os.path.join(data, "{name}.{mode}.faa.gz".format(name=name, mode=mode))
    with gzip.open(faa, "rt") as f:
        return list(parse(f))

def load_genes(name, mode):
    data = os.path.realpath(os.path.join(__file__, "..", "data"))
    fna = os.path.join(data, "{name}.{mode}.fna.gz".format(name=name, mode=mode))
    with gzip.open(fna, "rt") as f:
        return list(parse(f))
