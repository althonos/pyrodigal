import contextlib
import gzip
import io

try:
    try:
        from importlib.resources import files
    except ImportError:
        from importlib_resources import files  # type: ignore
except ImportError:
    files = None  # type: ignore

from ..fasta import parse

@contextlib.contextmanager
def load(name, mode="rt"):
     with files(__name__).joinpath(name).open("rb") as src:
        if name.endswith(".gz"):
            src = gzip.open(src, mode=mode)
        elif mode != "rb":
            src = io.TextIOWrapper(src)
        yield src 

def load_record(name):
    with load(name) as f:
        return next(parse(f))

def load_records(name):
    with load(name) as f:
        return list(parse(f))

def load_text(name):
    with load(name) as f:
        return f.read()