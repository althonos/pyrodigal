import collections.abc
import functools
import gzip
import os
import sys
import unittest
import random

from .. import TrainingInfo, Nodes, Sequence, _pyrodigal
from .._pyrodigal import METAGENOMIC_BINS, ConnectionScorer
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
    # compute scores only for some pseudo-randomly selected positions
    rng = random.Random(42)
    for i in rng.sample(range(len(nodes)), 10):
        # compute boundary (MAX_NODE_DIST = 500)
        j = 0 if i < 500 else i - 500
        # score connections without fast-indexing skippable nodes
        scorer.compute_skippable(j, i)
        scorer.score_connections(nodes, j, i, tinf, final=final)
    return nodes


class _TestConnectionScorerBase:
    backend = None

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

    @unittest.skipUnless(data.resources, "importlib.resources not available")
    def test_score_connections_final(self):
        record = data.load_record("MIIJ01000039.fna.gz")
        nodes_expected = scored_nodes(record, final=True, backend=None)
        nodes_actual = scored_nodes(record, final=True, backend=self.backend)
        for n1, n2 in zip(nodes_expected, nodes_actual):
            self.assertNodeEqual(n1, n2)

    @unittest.skipUnless(data.resources, "importlib.resources not available")
    def test_score_connections_train(self):
        record = data.load_record("GCF_001457455.1_NCTC11397_genomic.fna.gz")
        nodes_expected = scored_nodes(record, final=False, backend=None)
        nodes_actual = scored_nodes(record, final=False, backend=self.backend)
        for n1, n2 in zip(nodes_expected, nodes_actual):
            self.assertNodeEqual(n1, n2)


class TestConnectionScorerGeneric(_TestConnectionScorerBase, unittest.TestCase):
    backend = "generic"


@unittest.skipUnless(_pyrodigal._MMX_BUILD_SUPPORT, "extension compiled without MMX support")
@unittest.skipUnless(_pyrodigal._MMX_RUNTIME_SUPPORT, "requires machine with MMX support")
class TestConnectionScorerMMX(_TestConnectionScorerBase, unittest.TestCase):
    backend = "mmx"


@unittest.skipUnless(_pyrodigal._SSE2_BUILD_SUPPORT, "extension compiled without SSE2 support")
@unittest.skipUnless(_pyrodigal._SSE2_RUNTIME_SUPPORT, "requires machine with SSE2 support")
class TestConnectionScorerSSE(_TestConnectionScorerBase, unittest.TestCase):
    backend = "sse"


@unittest.skipUnless(_pyrodigal._AVX2_BUILD_SUPPORT, "extension compiled without AVX2 support")
@unittest.skipUnless(_pyrodigal._AVX2_RUNTIME_SUPPORT, "requires machine with AVX2 support")
class TestConnectionScorerAVX(_TestConnectionScorerBase, unittest.TestCase):
    backend = "avx"


@unittest.skipUnless(_pyrodigal._NEON_BUILD_SUPPORT, "extension compiled without NEON support")
@unittest.skipUnless(_pyrodigal._NEON_RUNTIME_SUPPORT, "requires machine with NEON support")
class TestConnectionScorerNEON(_TestConnectionScorerBase, unittest.TestCase):
    backend = "neon"
