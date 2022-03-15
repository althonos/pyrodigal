import collections.abc
import gzip
import os
import sys
import unittest

from .. import Nodes, Sequence, _pyrodigal
from .._pyrodigal import METAGENOMIC_BINS, ConnectionScorer

from .fasta import parse


class TestConnectionScorer(unittest.TestCase):

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

    @classmethod
    def setUpClass(cls):
        data = os.path.realpath(os.path.join(__file__, "..", "data"))
        fna = os.path.join(data, "MIIJ01000039.fna.gz")
        with gzip.open(fna, "rt") as f:
            cls.record = next(parse(f))

    @unittest.skipUnless(_pyrodigal._TARGET_CPU == "x86", "requires x86 CPU")
    @unittest.skipUnless(_pyrodigal._SSE2_BUILD_SUPPORT, "requires extension compiled with SSE2 support")
    @unittest.skipUnless(_pyrodigal._SSE2_RUNTIME_SUPPORT, "requires machine with SSE2 support")
    def test_score_connections_sse(self):
        # setup
        seq = Sequence.from_string(self.record.seq)
        tinf = METAGENOMIC_BINS[0].training_info
        scorer_sse = ConnectionScorer(backend="sse")
        scorer_none = ConnectionScorer(backend=None)
        # add nodes from the sequence
        nodes = Nodes()
        nodes.extract(seq, translation_table=tinf.translation_table)
        nodes.sort()
        # index nodes for the scorers
        scorer_sse.index(nodes)
        scorer_none.index(nodes)
        # use copies to compute both scores
        nodes_sse = nodes.copy()
        nodes_none = nodes.copy()
        for i in range(500, len(nodes)):
            # compute boundary
            j = 0 if i < 500 else i - 500
            # score connections without fast-indexing skippable nodes
            scorer_none.compute_skippable(j, i)
            scorer_none.score_connections(nodes_none, j, i, tinf, final=True)
            # compute skippable nodes with SSE and score connections with
            scorer_sse.compute_skippable(j, i)
            scorer_sse.score_connections(nodes_sse, j, i, tinf, final=True)
        # check that both methods scored the same
        for n1, n2 in zip(nodes_sse, nodes_none):
            self.assertNodeEqual(n1, n2)

    @unittest.skipUnless(_pyrodigal._TARGET_CPU == "x86", "requires x86 CPU")
    @unittest.skipUnless(_pyrodigal._AVX2_BUILD_SUPPORT, "requires extension compiled with AVX2 support")
    @unittest.skipUnless(_pyrodigal._AVX2_RUNTIME_SUPPORT, "requires machine with AVX2 support")
    def test_score_connections_avx(self):
        # setup
        seq = Sequence.from_string(self.record.seq)
        tinf = METAGENOMIC_BINS[0].training_info
        scorer_avx = ConnectionScorer(backend="avx")
        scorer_none = ConnectionScorer(backend=None)
        # add nodes from the sequence
        nodes = Nodes()
        nodes.extract(seq, translation_table=tinf.translation_table)
        nodes.sort()
        # index nodes for the scorers
        scorer_avx.index(nodes)
        scorer_none.index(nodes)
        # use copies to compute both scores
        nodes_avx = nodes.copy()
        nodes_none = nodes.copy()
        for i in range(500, len(nodes)):
            # compute boundary
            j = 0 if i < 500 else i - 500
            # score connections without fast-indexing skippable nodes
            scorer_none.compute_skippable(j, i)
            scorer_none.score_connections(nodes_none, j, i, tinf, final=True)
            # compute skippable nodes with SSE and score connections with
            scorer_avx.compute_skippable(j, i)
            scorer_avx.score_connections(nodes_avx, j, i, tinf, final=True)
        # check that both methods scored the same
        for n1, n2 in zip(nodes_avx, nodes_none):
            self.assertNodeEqual(n1, n2)

    @unittest.skipUnless(_pyrodigal._TARGET_CPU in ("arm", "aarch64"), "requires ARM CPU")
    @unittest.skipUnless(_pyrodigal._NEON_BUILD_SUPPORT, "requires extension compiled with NEON support")
    @unittest.skipUnless(_pyrodigal._NEON_RUNTIME_SUPPORT, "requires machine with NEON support")
    def test_score_connections_neon(self):
        # setup
        seq = Sequence.from_string(self.record.seq)
        tinf = METAGENOMIC_BINS[0].training_info
        scorer_avx = ConnectionScorer(backend="neon")
        scorer_none = ConnectionScorer(backend=None)
        # add nodes from the sequence
        nodes = Nodes()
        nodes.extract(seq, translation_table=tinf.translation_table)
        nodes.sort()
        # index nodes for the scorers
        scorer_avx.index(nodes)
        scorer_none.index(nodes)
        # use copies to compute both scores
        nodes_avx = nodes.copy()
        nodes_none = nodes.copy()
        for i in range(500, len(nodes)):
            # compute boundary
            j = 0 if i < 500 else i - 500
            # score connections without fast-indexing skippable nodes
            scorer_none.compute_skippable(j, i)
            scorer_none.score_connections(nodes_none, j, i, tinf, final=True)
            # compute skippable nodes with SSE and score connections with
            scorer_avx.compute_skippable(j, i)
            scorer_avx.score_connections(nodes_avx, j, i, tinf, final=True)
        # check that both methods scored the same
        for n1, n2 in zip(nodes_avx, nodes_none):
            self.assertNodeEqual(n1, n2)

    def test_score_connections_generic(self):
        # setup
        seq = Sequence.from_string(self.record.seq)
        tinf = METAGENOMIC_BINS[0].training_info
        scorer_generic = ConnectionScorer(backend="generic")
        scorer_none = ConnectionScorer(backend=None)
        # add nodes from the sequence
        nodes = Nodes()
        nodes.extract(seq, translation_table=tinf.translation_table)
        nodes.sort()
        # index nodes for the scorers
        scorer_generic.index(nodes)
        scorer_none.index(nodes)
        # use copies to compute both scores
        nodes_generic = nodes.copy()
        nodes_none = nodes.copy()
        for i in range(500, len(nodes)):
            # compute boundary
            j = 0 if i < 500 else i - 500
            # score connections without fast-indexing skippable nodes
            scorer_none.compute_skippable(j, i)
            scorer_none.score_connections(nodes_none, j, i, tinf, final=True)
            # compute skippable nodes with generic filter and score connections
            scorer_generic.compute_skippable(j, i)
            scorer_generic.score_connections(nodes_generic, j, i, tinf, final=True)
        # check that both methods scored the same
        for n1, n2 in zip(nodes_generic, nodes_none):
            self.assertNodeEqual(n1, n2)
