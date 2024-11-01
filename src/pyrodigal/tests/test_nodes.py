import collections.abc
import gzip
import os
import sys
import pickle
import unittest

from .. import Nodes, Sequence, METAGENOMIC_BINS
from . import data
from .fasta import parse


class TestNodes(unittest.TestCase):
    def assertNodeEqual(self, n1, n2):
        self.assertEqual(n1.index, n2.index, "indices differ")
        self.assertEqual(n1.strand, n2.strand, "strands differ")
        self.assertEqual(n1.type, n2.type, "types differ")
        self.assertEqual(n1.edge, n2.edge, "edge differ")
        self.assertEqual(n1.gc_bias, n2.gc_bias, "GC biases differ")
        self.assertEqual(n1.cscore, n2.cscore, "cscores differ")
        self.assertEqual(n1.gc_cont, n2.gc_cont, "GC contents differ")
        self.assertEqual(n1.score, n2.score, "GC contents differ")
        self.assertEqual(n1.rscore, n2.rscore, "rscores differ")
        self.assertEqual(n1.sscore, n2.sscore, "sscores differ")
        self.assertEqual(n1.tscore, n2.tscore, "tscores differ")

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_add_nodes_metagenomic_bins(self):
        record = data.load_record("SRR492066.fna.gz")
        seq = Sequence(record.seq)
        nodes = Nodes()
        # nodes should start empty
        self.assertEqual(len(nodes), 0)
        # numbers below obtained directly in Prodigal by `printf`-ing the
        # node numbers on a normal run
        for bin, expected in [(0, 2970), (2, 2970), (11, 2293), (24, 2293)]:
            tt = METAGENOMIC_BINS[bin].training_info.translation_table
            self.assertEqual(nodes.extract(seq, translation_table=tt), expected)
            self.assertEqual(len(nodes), expected)
            nodes.clear()

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_copy(self):
        record = data.load_record("SRR492066.fna.gz")
        tt = METAGENOMIC_BINS[0].training_info.translation_table
        seq = Sequence(record.seq)
        nodes1 = Nodes()
        nodes1.extract(seq, translation_table=tt)
        nodes2 = nodes1.copy()
        for n1, n2 in zip(nodes1, nodes2):
            self.assertNodeEqual(n1, n2)

    def test_copy_empty(self):
        nodes = Nodes()
        copy = nodes.copy()
        self.assertEqual(len(nodes), 0)
        self.assertEqual(len(copy), 0)

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_pickle(self):
        record = data.load_record("SRR492066.fna.gz")
        tt = METAGENOMIC_BINS[0].training_info.translation_table
        seq = Sequence(record.seq)
        nodes1 = Nodes()
        nodes1.extract(seq, translation_table=tt)
        nodes2 = pickle.loads(pickle.dumps(nodes1))
        self.assertEqual(len(nodes1), len(nodes2), "lengths differ")
        for n1, n2 in zip(nodes1, nodes2):
            self.assertNodeEqual(n1, n2)

    def test_pickle_empty(self):
        nodes1 = Nodes()
        nodes2 = pickle.loads(pickle.dumps(nodes1))
        self.assertEqual(len(nodes1), 0)
        self.assertEqual(len(nodes2), 0)

    def test_extract_edge_start(self):
        # make sure that start nodes on edges are not extracted twice
        # when in open genome mode as it was a bug at some point (#22)
        nodes = Nodes()
        seq = Sequence(
            "ATGGTTAACGCTTCCGGCGACCCCGTAATCGAGGCCGCCC"   # "ATG" at index 0
            "ACATCTGGTCAGACACGCTGACGGTGCTCAAACACAGCGC"
            "TTCGCTCAGCCCACGAGAAAAAGGCTGGTTGGAAGGCGTT"
            "GTTCCTGAAGGCGTCTTCGGTTCGACCATCGTGCTGTGTG"
            "TGGACAACAACGACACGCTTCAAGCCATTCAGGGTGATTT"
            "GAACGATTCCCTGCTTCAGGCATTGCGTACGGTAACCGGC"
            "GAAAATATGTTTCCCGCGTTCAAGGTCGTGCCGAAAACCG"
        )
        nodes.extract(seq, closed=False)
        nodes.sort()
        # check the first node is *not* an edge not despite being on
        # the edge and the node extraction running in open genome mode
        self.assertEqual(nodes[0].index, 0)
        self.assertFalse(nodes[0].edge)
        self.assertEqual(nodes[0].strand, 1)
        self.assertEqual(nodes[0].type, "ATG")
        # check the second node is *not* at index 0
        self.assertNotEqual(nodes[1].index, 0)
