import collections.abc
import gzip
import os
import sys
import unittest

import Bio.SeqIO
from pyrodigal import Pyrodigal



class TestGenes(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        data = os.path.realpath(os.path.join(__file__, "..", "data"))
        fna = os.path.join(data, "SRR492066.fna.gz")
        meta_faa = os.path.join(data, "SRR492066.meta.faa.gz")

        with gzip.open(fna, "rt") as f:
            cls.record = next(Bio.SeqIO.parse(f, "fasta"))

        cls.p = Pyrodigal(meta=True)
        cls.genes = cls.p.find_genes(str(cls.record.seq))

    def test_indexing(self):
        length = len(self.genes)
        self.assertEqual(self.genes[0]._data, self.genes[-length]._data)
        with self.assertRaises(IndexError):
            self.genes[length]
        with self.assertRaises(IndexError):
            self.genes[-length-1]

    def test_iter(self):
        for i, gene in zip(range(len(self.genes)), self.genes):
            self.assertEqual(gene._data, self.genes[i]._data)

    def test_reversed(self):
        for i, gene in zip(range(1, len(self.genes)+1), reversed(self.genes)):
            self.assertEqual(gene._data, self.genes[-i]._data)

    def test_bool(self):
        self.assertTrue(bool(self.genes))
        self.assertFalse(bool(self.p.find_genes("TTT")))

    def test_translation_table(self):
        valid = set((*range(1, 7), *range(9, 17), *range(21, 26)))
        for gene in self.genes:
            self.assertIn(gene.translation_table, valid)

    @unittest.skipIf(sys.implementation.name != 'cpython', 'can panic with PyPy')
    def test_collection_abc_subclass(self):
        self.assertIsInstance(self.genes, collections.abc.Sequence)
        self.assertIsInstance(self.genes, collections.abc.Sized)
        self.assertIsInstance(self.genes, collections.abc.Container)
        self.assertIsInstance(self.genes, collections.abc.Iterable)
        if sys.version_info >= (3, 6):
            self.assertIsInstance(self.genes, collections.abc.Reversible)
