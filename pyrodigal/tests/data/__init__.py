import contextlib
import gzip
import io

try:
    try:
        from importlib import resources
    except ImportError:
        import importlib_resources as resources
except ImportError:
    resources = None

from ..fasta import parse

@contextlib.contextmanager
def load(name, mode="rt"):
     with resources.open_binary(__name__, name) as src:
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