import collections.abc
import gzip
import os
import sys
import unittest

from .. import Nodes, Sequence, _pyrodigal
from .._pyrodigal import METAGENOMIC_BINS, ConnectionScorer

from .fasta import parse


class TestConnectionScorer(unittest.TestCase):

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
        scorer_generic = ConnectionScorer(backend=None)
        # add nodes from the sequence
        nodes = Nodes()
        nodes.extract(seq, translation_table=tinf.translation_table)
        nodes.sort()
        # index nodes for the scorers
        scorer_sse.index(nodes)
        scorer_generic.index(nodes)
        # use copies to compute both scores
        nodes_sse = nodes.copy()
        nodes_generic = nodes.copy()
        for i in range(500, len(nodes)):
            # compute boundary
            j = 0 if i < 500 else i - 500
            # score connections without fast-indexing skippable nodes
            scorer_generic.compute_skippable(j, i)
            scorer_generic.score_connections(nodes_generic, j, i, tinf, final=True)
            # compute skippable nodes with SSE and score connections with
            scorer_sse.compute_skippable(j, i)
            scorer_sse.score_connections(nodes_sse, j, i, tinf, final=True)
        # check that both methods scored the same
        for n_sse, n_generic in zip(nodes_sse, nodes_generic):
            self.assertEqual(n_sse.score, n_generic.score)

    @unittest.skipUnless(_pyrodigal._TARGET_CPU == "x86", "requires x86 CPU")
    @unittest.skipUnless(_pyrodigal._AVX2_BUILD_SUPPORT, "requires extension compiled with AVX2 support")
    @unittest.skipUnless(_pyrodigal._AVX2_RUNTIME_SUPPORT, "requires machine with AVX2 support")
    def test_score_connections_avx(self):
        # setup
        seq = Sequence.from_string(self.record.seq)
        tinf = METAGENOMIC_BINS[0].training_info
        scorer_avx = ConnectionScorer(backend="avx")
        scorer_generic = ConnectionScorer(backend=None)
        # add nodes from the sequence
        nodes = Nodes()
        nodes.extract(seq, translation_table=tinf.translation_table)
        nodes.sort()
        # index nodes for the scorers
        scorer_avx.index(nodes)
        scorer_generic.index(nodes)
        # use copies to compute both scores
        nodes_avx = nodes.copy()
        nodes_generic = nodes.copy()
        for i in range(500, len(nodes)):
            # compute boundary
            j = 0 if i < 500 else i - 500
            # score connections without fast-indexing skippable nodes
            scorer_generic.compute_skippable(j, i)
            scorer_generic.score_connections(nodes_generic, j, i, tinf, final=True)
            # compute skippable nodes with SSE and score connections with
            scorer_avx.compute_skippable(j, i)
            scorer_avx.score_connections(nodes_avx, j, i, tinf, final=True)
        # check that both methods scored the same
        for n_avx, n_generic in zip(nodes_avx, nodes_generic):
            self.assertEqual(n_avx.score, n_generic.score)

    @unittest.skipUnless(_pyrodigal._TARGET_CPU in ("arm", "aarch64"), "requires ARM CPU")
    @unittest.skipUnless(_pyrodigal._NEON_BUILD_SUPPORT, "requires extension compiled with NEON support")
    @unittest.skipUnless(_pyrodigal._NEON_RUNTIME_SUPPORT, "requires machine with NEON support")
    def test_score_connections_neon(self):
        # setup
        seq = Sequence.from_string(self.record.seq)
        tinf = METAGENOMIC_BINS[0].training_info
        scorer_avx = ConnectionScorer(backend="neon")
        scorer_generic = ConnectionScorer(backend=None)
        # add nodes from the sequence
        nodes = Nodes()
        nodes.extract(seq, translation_table=tinf.translation_table)
        nodes.sort()
        # index nodes for the scorers
        scorer_avx.index(nodes)
        scorer_generic.index(nodes)
        # use copies to compute both scores
        nodes_avx = nodes.copy()
        nodes_generic = nodes.copy()
        for i in range(500, len(nodes)):
            # compute boundary
            j = 0 if i < 500 else i - 500
            # score connections without fast-indexing skippable nodes
            scorer_generic.compute_skippable(j, i)
            scorer_generic.score_connections(nodes_generic, j, i, tinf, final=True)
            # compute skippable nodes with SSE and score connections with
            scorer_avx.compute_skippable(j, i)
            scorer_avx.score_connections(nodes_avx, j, i, tinf, final=True)
        # check that both methods scored the same
        for n_avx, n_generic in zip(nodes_avx, nodes_generic):
            self.assertEqual(n_avx.score, n_generic.score)
