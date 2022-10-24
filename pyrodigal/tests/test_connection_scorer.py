import collections.abc
import gzip
import os
import sys
import unittest

from .. import TrainingInfo, Nodes, Sequence, _pyrodigal
from .._pyrodigal import METAGENOMIC_BINS, ConnectionScorer

from .fasta import parse


class TestConnectionScorer(unittest.TestCase):
    backend = "generic"

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
        fna_train = os.path.join(data, "GCF_000009045.1_ASM904v1_genomic.fna.gz")
        with gzip.open(fna_train, "rt") as f:
            cls.record_train = next(parse(f))

    def test_score_connections_final(self):
        # setup
        seq = Sequence.from_string(self.record.seq)
        tinf = METAGENOMIC_BINS[0].training_info
        scorer_simd = ConnectionScorer(backend=self.backend)
        scorer_none = ConnectionScorer(backend=None)
        # add nodes from the sequence
        nodes = Nodes()
        nodes.extract(seq, translation_table=tinf.translation_table)
        nodes.sort()
        # index nodes for the scorers
        scorer_simd.index(nodes)
        scorer_none.index(nodes)
        # use copies to compute both scores
        nodes_simd = nodes.copy()
        nodes_none = nodes.copy()
        for i in range(500, len(nodes)):
            # compute boundary (MAX_NODE_DIST = 500)
            j = 0 if i < 500 else i - 500
            # score connections without fast-indexing skippable nodes
            scorer_none.compute_skippable(j, i)
            scorer_none.score_connections(nodes_none, j, i, tinf, final=True)
            # compute skippable nodes with SIMD and then score connections
            scorer_simd.compute_skippable(j, i)
            scorer_simd.score_connections(nodes_simd, j, i, tinf, final=True)
        # check that both methods scored the same
        for n1, n2 in zip(nodes_none, nodes_simd):
            self.assertNodeEqual(n1, n2)

    def test_score_connections_training(self):
        # setup
        seq = Sequence.from_string(self.record_train.seq)
        tinf = TrainingInfo(seq.gc, translation_table=11)
        scorer_simd = ConnectionScorer(backend=self.backend)
        scorer_none = ConnectionScorer(backend=None)
        # add nodes from the sequence
        nodes = Nodes()
        nodes.extract(seq, closed=True)
        nodes.sort()
        self.assertEqual(len(nodes), 196059)
        # index nodes for the scorers
        scorer_simd.index(nodes)
        scorer_none.index(nodes)
        # use copies to compute both scores
        nodes_simd = nodes.copy()
        nodes_none = nodes.copy()
        for i in range(500, len(nodes)):
            # compute boundary (MAX_NODE_DIST = 500)
            j = 0 if i < 500 else i - 500
            # score connections without fast-indexing skippable nodes
            scorer_none.compute_skippable(j, i)
            scorer_none.score_connections(nodes_none, j, i, tinf, final=False)
            # compute skippable nodes with SIMD and then score connections
            scorer_simd.compute_skippable(j, i)
            scorer_simd.score_connections(nodes_simd, j, i, tinf, final=False)
        # check that both methods scored the same
        for n1, n2 in zip(nodes_none, nodes_simd):
            self.assertNodeEqual(n1, n2)


@unittest.skipUnless(_pyrodigal._MMX_BUILD_SUPPORT, "extension compiled without MMX support")
@unittest.skipUnless(_pyrodigal._MMX_RUNTIME_SUPPORT, "requires machine with MMX support")
class TestConnectionScorerMMX(TestConnectionScorer):
    backend = "mmx"


@unittest.skipUnless(_pyrodigal._SSE2_BUILD_SUPPORT, "extension compiled without SSE2 support")
@unittest.skipUnless(_pyrodigal._SSE2_RUNTIME_SUPPORT, "requires machine with SSE2 support")
class TestConnectionScorerSSE(TestConnectionScorer):
    backend = "sse"


@unittest.skipUnless(_pyrodigal._AVX2_BUILD_SUPPORT, "extension compiled without AVX2 support")
@unittest.skipUnless(_pyrodigal._AVX2_RUNTIME_SUPPORT, "requires machine with AVX2 support")
class TestConnectionScorerAVX(TestConnectionScorer):
    backend = "avx"


@unittest.skipUnless(_pyrodigal._NEON_BUILD_SUPPORT, "extension compiled without NEON support")
@unittest.skipUnless(_pyrodigal._NEON_RUNTIME_SUPPORT, "requires machine with NEON support")
class TestConnectionScorerNEON(TestConnectionScorer):
    backend = "neon"
