import collections.abc
import functools
import gzip
import os
import sys
import unittest
import random
import typing

from .. import TrainingInfo, Nodes, Sequence, METAGENOMIC_BINS, lib
from ..lib import ConnectionScorer
from . import data


@functools.lru_cache()
def extract_nodes(record):
    seq = Sequence(record.seq)
    tinf = METAGENOMIC_BINS[0].training_info
    nodes = Nodes()
    nodes.extract(seq, translation_table=tinf.translation_table)
    nodes.sort()
    return nodes


@functools.lru_cache()
def scored_nodes(record, final=True, backend=None):
    # extract nodes from the record
    tinf = METAGENOMIC_BINS[0].training_info
    scorer = ConnectionScorer(backend=backend)
    nodes = extract_nodes(record).copy()
    scorer.index(nodes)
    # compute scores for all positions
    scorer.score_connections(nodes, tinf, final=final)
    return nodes


class _TestConnectionScorerBase:
    backend = None  # type: typing.Optional[str]

    def assertNodeEqual(self, n1, n2):
        self.assertEqual(n1.index, n2.index)
        self.assertEqual(n1.strand, n2.strand)
        self.assertEqual(n1.type, n2.type)
        self.assertEqual(n1.edge, n2.edge)
        self.assertEqual(n1.gc_bias, n2.gc_bias)
        self.assertEqual(n1.gc_cont, n2.gc_cont)
        self.assertAlmostEqual(n1.score, n2.score)
        self.assertAlmostEqual(n1.cscore, n2.cscore)
        self.assertAlmostEqual(n1.rscore, n2.rscore)
        self.assertAlmostEqual(n1.sscore, n2.sscore)
        self.assertAlmostEqual(n1.tscore, n2.tscore)

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_score_connections_final(self):
        record = data.load_record("MIIJ01000039.fna.gz")
        nodes_expected = scored_nodes(record, final=True, backend=None)
        nodes_actual = scored_nodes(record, final=True, backend=self.backend)
        for n1, n2 in zip(nodes_expected, nodes_actual):
            self.assertNodeEqual(n1, n2)

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_score_connections_train(self):
        record = data.load_record("GCF_001457455.1_NCTC11397_genomic.fna.gz")
        nodes_expected = scored_nodes(record, final=False, backend=None)
        nodes_actual = scored_nodes(record, final=False, backend=self.backend)
        for n1, n2 in zip(nodes_expected, nodes_actual):
            self.assertNodeEqual(n1, n2)


class TestConnectionScorerGeneric(_TestConnectionScorerBase, unittest.TestCase):
    backend = "generic"


@unittest.skipUnless(lib._MMX_BUILD_SUPPORT, "extension compiled without MMX support")
@unittest.skipUnless(lib._MMX_RUNTIME_SUPPORT, "requires machine with MMX support")
class TestConnectionScorerMMX(_TestConnectionScorerBase, unittest.TestCase):
    backend = "mmx"


@unittest.skipUnless(lib._SSE2_BUILD_SUPPORT, "extension compiled without SSE2 support")
@unittest.skipUnless(lib._SSE2_RUNTIME_SUPPORT, "requires machine with SSE2 support")
class TestConnectionScorerSSE(_TestConnectionScorerBase, unittest.TestCase):
    backend = "sse"


@unittest.skipUnless(lib._AVX2_BUILD_SUPPORT, "extension compiled without AVX2 support")
@unittest.skipUnless(lib._AVX2_RUNTIME_SUPPORT, "requires machine with AVX2 support")
class TestConnectionScorerAVX(_TestConnectionScorerBase, unittest.TestCase):
    backend = "avx"


@unittest.skipUnless(lib._AVX512_BUILD_SUPPORT, "extension compiled without AVX512 support")
@unittest.skipUnless(lib._AVX512_RUNTIME_SUPPORT, "requires machine with AVX512 support")
class TestConnectionScorerAVX512(_TestConnectionScorerBase, unittest.TestCase):
    backend = "avx512"


@unittest.skipUnless(lib._NEON_BUILD_SUPPORT, "extension compiled without NEON support")
@unittest.skipUnless(lib._NEON_RUNTIME_SUPPORT, "requires machine with NEON support")
class TestConnectionScorerNEON(_TestConnectionScorerBase, unittest.TestCase):
    backend = "neon"
