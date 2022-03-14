import collections.abc
import gzip
import os
import sys
import unittest

from .. import Nodes, Sequence
from .._pyrodigal import METAGENOMIC_BINS

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
        # nodes should start empty
        self.assertEqual(len(nodes), 0)
        # numbers below obtained directly in Prodigal by `printf`-ing the
        # node numbers on a normal run
        for bin, expected in [(0, 2970), (2, 2970), (11, 2293), (24, 2293)]:
            tt = METAGENOMIC_BINS[bin].training_info.translation_table
            self.assertEqual(nodes.extract(seq, translation_table=tt), expected)
            self.assertEqual(len(nodes), expected)
            nodes.clear()

    def test_copy(self):
        tt = METAGENOMIC_BINS[0].training_info.translation_table
        seq = Sequence.from_string(self.record.seq)
        nodes1 = Nodes()
        nodes1.extract(seq, translation_table=tt)
        nodes2 = nodes1.copy()
        for n1, n2 in zip(nodes1, nodes2):
            self.assertEqual(n1.type, n2.type)

    def test_copy_empty(self):
        nodes = Nodes()
        copy = nodes.copy()
        self.assertEqual(len(nodes), 0)
        self.assertEqual(len(copy), 0)
