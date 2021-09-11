import collections.abc
import gzip
import os
import sys
import unittest

from .. import Nodes, Sequence
from .._pyrodigal import METAGENOMIC_BINS, add_nodes

from .fasta import parse


class TestNodes(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        data = os.path.realpath(os.path.join(__file__, "..", "data"))
        fna = os.path.join(data, "SRR492066.fna.gz")
        with gzip.open(fna, "rt") as f:
            cls.record = next(parse(f))

    def test_add_nodes_metagenomic_bins(self):
        seq = Sequence.from_string(self.record.seq)
        nodes = Nodes()

        self.assertEqual(len(nodes), 0)

        self.assertEqual(add_nodes(nodes, seq, METAGENOMIC_BINS[0].training_info), 2970)
        self.assertEqual(len(nodes), 2970)
        nodes.clear()

        self.assertEqual(add_nodes(nodes, seq, METAGENOMIC_BINS[2].training_info), 2970)
        self.assertEqual(len(nodes), 2970)
        nodes.clear()

        self.assertEqual(add_nodes(nodes, seq, METAGENOMIC_BINS[11].training_info), 2293)
        self.assertEqual(len(nodes), 2293)
        nodes.clear()

        self.assertEqual(add_nodes(nodes, seq, METAGENOMIC_BINS[24].training_info), 2293)
        self.assertEqual(len(nodes), 2293)
        nodes.clear()
