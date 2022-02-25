# Contributing to Pyrodigal

For bug fixes or new features, please file an issue before submitting a
pull request. If the change isn't trivial, it may be best to wait for
feedback.

## Setting up a local repository

Make sure you clone the repository in recursive mode, so you also get the
wrapped code of Prodigal which is exposed as a ``git`` submodule:

```console
$ git clone --recursive https://github.com/althonos/pyrodigal
```

## Running tests

Tests are written as usual Python unit tests with the `unittest` module of
the standard library. Running them requires the extension to be built
locally:

```console
$ python setup.py build_ext --debug --inplace
$ python -m unittest discover -vv
```

## Running benchmarks

### Connection scoring

The `benches` folder contains benchmarks for evaluating the performance of
the node connection scoring step, essentially to make sure that the SIMD
code makes it faster. To run the benchmark you need a CPU that supports
SSE2 and AVX2.

Start by building `pyrodigal` locally:
```console
$ python setup.py build_ext --debug --inplace
```

Then make sure you have the required packages and data:
```console
$ pip install --user -r benches/connection_scoring/requirements.txt
$ python benches/connection_scoring/data/download.py
```

Finally, run the benchmarks and plot the results:
```console
$ python benches/connection_scoring/bench.py -d benches/connection_scoring/data/ -o times.json
$ python benches/connection_scoring/plot.py -i times.json --show
```

## Coding guidelines

This project targets Python 3.5 or later.

Python objects should be typed; since it is not supported by Cython,
you must manually declare types in type stubs (`.pyi` files). In Python
files, you can add type annotations to function signatures (supported in
Python 3.5) but not in variable assignments (supported only from Python
3.6 onward). However, Cython allows you to use [f-strings](https://www.python.org/dev/peps/pep-0498/)
even when compiling the code for Python 3.5.

### Interfacing with C

When interfacing with C, and in particular with pointers, use assertions
everywhere you assume the pointer to be non-NULL. Also consider using
assertions when accessing raw C arrays, if applicable.
