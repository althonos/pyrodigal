import bz2
import collections
import contextlib
import io
import os



_BZ2_MAGIC = b"BZh"
_GZIP_MAGIC = b"\x1f\x8b"
_XZ_MAGIC = b"\xfd7zXZ"
_LZ4_MAGIC = b"\x04\x22\x4d\x18"
_ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"


@contextlib.contextmanager
def zopen(file, mode='r', encoding=None, errors=None, newline=None):
    with contextlib.ExitStack() as ctx:

        try:
            path = os.fsencode(file)
            file = ctx.enter_context(open(path, mode="rb"))
        except TypeError:
            raise

        peek = file.peek()
        if peek.startswith(_GZIP_MAGIC):
            try:
                from isal import igzip as gzip
            except ImportError:
                import gzip  # type: ignore
            file = ctx.enter_context(gzip.open(file, mode="rb"))
        elif peek.startswith(_BZ2_MAGIC):
            import bz2
            file = ctx.enter_context(bz2.open(file, mode="rb"))
        elif peek.startswith(_XZ_MAGIC):
            import lzma
            file = ctx.enter_context(lzma.open(file, mode="rb"))
        elif peek.startswith(_LZ4_MAGIC):
            try:
                import lz4.frame
            except ImportError as err:
                raise RuntimeError("File compression is LZ4 but lz4 is not installed") from err
            file = ctx.enter_context(lz4.frame.open(file))
        elif peek.startswith(_ZSTD_MAGIC):
            try:
                import zstandard
            except ImportError as err:
                raise RuntimeError("File compression is ZSTD but zstandard is not installed") from err
            decompressor = zstandard.ZstdDecompressor()
            file = decompressor.stream_reader(file)
        if mode == "r":
            file = io.TextIOWrapper(file, encoding=encoding, errors=errors, newline=newline)
        yield file


class Record(collections.namedtuple("Record", ["id", "seq", "description"])):
    pass


def parse(path):
    with contextlib.ExitStack() as ctx:
        try:
            path = os.fsencode(path)
            file = ctx.enter_context(zopen(path, "r"))
        except TypeError:
            file = path

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
