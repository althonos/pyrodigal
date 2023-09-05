import bz2
import collections
import io
import os

try:
    from isal import igzip as gzip
except ImportError:
    import gzip  # type: ignore

_GZIP_MAGIC = b'\x1f\x8b'
_BZ2_MAGIC = b'BZh'


class Record(collections.namedtuple("Record", ["id", "seq", "description"])):
    pass


def parse(path):
    try:
        path = os.fsencode(path)
        file = open(path, "rb")
        b = file.peek()
        if b.startswith(_GZIP_MAGIC):
            file = gzip.open(file, mode="rt")
        elif b.startswith(_BZ2_MAGIC):
            file = bz2.open(file, mode="rt")
        else:
            file = io.TextIOWrapper(file)
    except TypeError:
        file = path

    with file:
        # parse file
        id_ = None
        seq = []
        for line in file:
            l = line.strip()
            if line.startswith(">"):
                if id_ is not None:
                    yield Record(id_, "".join(seq), desc)
                id_ = line[1:].split()[0].strip()
                desc = " ".join(line[1:].split(maxsplit=1))
                seq = []
            elif l:
                seq.append(l)
        if id_ is not None:
            yield Record(id_, "".join(seq), desc)
        elif seq:
            raise ValueError("not in FASTA format")
