import abc
import gzip
import os
import textwrap
import unittest
import warnings

from .. import OrfFinder
from . import data


@unittest.skipUnless(data.resources, "importlib.resources not available")
class TestGene(unittest.TestCase):
    @classmethod
    def find_genes(cls, seq):
        p = OrfFinder(meta=True)
        return p.find_genes(seq)

    @classmethod
    def setUpClass(cls):
        cls.record = data.load_record("SRR492066.fna.gz")
        cls.preds = cls.find_genes(str(cls.record.seq))
        cls.proteins = [
            record 
            for record in data.load_records("SRR492066.meta.faa.gz")
            if record.id.startswith("{}_".format(cls.record.id)) 
        ]
        cls.genes = [
            record 
            for record in data.load_records("SRR492066.meta.fna.gz")
            if record.id.startswith("{}_".format(cls.record.id)) 
        ]

    def test_translate(self):
        gene = self.preds[0]
        prot = self.proteins[0].seq
        self.assertEqual(self.preds[0].translation_table, 4)
        self.assertEqual(gene.translate(), prot)
        self.assertEqual(gene.translate(translation_table=4), prot)
        self.assertNotEqual(gene.translate(translation_table=2), prot)
        self.assertEqual(gene.translate(), prot)

    def test_translate_error(self):
        gene = self.preds[0]
        self.assertRaises(ValueError, gene.translate, 7)

    def test_sequence(self):
        for gene in self.preds:
            if gene.strand == 1:
                self.assertEqual(gene.sequence(), self.record.seq[gene.begin-1:gene.end])