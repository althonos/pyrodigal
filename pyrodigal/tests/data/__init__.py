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

from ..fasta import parse, zopen

@contextlib.contextmanager
def load(name, mode="rb"):
     with zopen(files(__name__).joinpath(name), mode=mode) as src:
        yield src

def load_record(name):
    with load(name, mode="r") as f:
        return next(parse(f))

def load_records(name):
    with load(name, mode="r") as f:
        return list(parse(f))

def load_text(name):
    with load(name, mode="r") as f:
        return f.read()