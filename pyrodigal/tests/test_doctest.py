# coding: utf-8
"""Test doctest contained tests in every file of the module.
"""

import configparser
import doctest
import importlib
import json
import gzip
import os
import pkgutil
import re
import shutil
import sys
import types
import warnings
from unittest import mock

import pyrodigal

from .fasta import parse


def _load_tests_from_module(tests, module, globs, setUp=None, tearDown=None):
    """Load tests from module, iterating through submodules."""
    for attr in (getattr(module, x) for x in dir(module) if not x.startswith("_")):
        if isinstance(attr, types.ModuleType):
            suite = doctest.DocTestSuite(
                attr,
                globs,
                setUp=setUp,
                tearDown=tearDown,
                optionflags=+doctest.ELLIPSIS,
            )
            tests.addTests(suite)
    return tests


def load_tests(loader, tests, ignore):
    """`load_test` function used by unittest to find the doctests."""
    _current_cwd = os.getcwd()
    # demonstrate how to load sequences with Biopython without requiring Biopython
    Bio = mock.Mock()
    Bio.SeqIO = mock.Mock()
    Bio.SeqIO.read = lambda file, format: next(parse(file))

    def setUp(self):
        warnings.simplefilter("ignore")
        os.chdir(os.path.realpath(os.path.join(__file__, os.path.pardir, "data")))
        sys.modules["Bio"] = Bio
        sys.modules["Bio.SeqIO"] = Bio.SeqIO

    def tearDown(self):
        os.chdir(_current_cwd)
        warnings.simplefilter(warnings.defaultaction)
        sys.modules.pop("Bio")
        sys.modules.pop("Bio.SeqIO")

    # doctests are not compatible with `green`, so we may want to bail out
    # early if `green` is running the tests
    if sys.argv[0].endswith("green"):
        return tests

    # doctests require Bio.SeqIO for loading data, which may not be available
    if Bio is None:
        return tests

    # recursively traverse all library submodules and load tests from them
    packages = [None, pyrodigal]
    for pkg in iter(packages.pop, None):
        for (_, subpkgname, subispkg) in pkgutil.walk_packages(pkg.__path__):
            # do not import __main__ module to avoid side effects!
            if subpkgname == "__main__" or subpkgname.startswith("tests") or subpkgname.startswith("cli"):
                continue
            # import the submodule and add it to the tests
            module = importlib.import_module(".".join([pkg.__name__, subpkgname]))
            globs = dict(pyrodigal=pyrodigal, json=json, gzip=gzip, Bio=Bio, **module.__dict__)
            tests.addTests(
                doctest.DocTestSuite(
                    module,
                    globs=globs,
                    setUp=setUp,
                    tearDown=tearDown,
                    optionflags=+doctest.ELLIPSIS,
                )
            )
            # if the submodule is a package, we need to process its submodules
            # as well, so we add it to the package queue
            if subispkg and subpkgname != "tests":
                packages.append(module)

    return tests
