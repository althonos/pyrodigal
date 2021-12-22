import abc
import gzip
import os
import textwrap
import unittest
import warnings

from .. import OrfFinder
from .fasta import parse


class TestGene(unittest.TestCase):

    @classmethod
    def find_genes(cls, seq):
        p = OrfFinder(meta=True)
        return p.find_genes(seq)

    @classmethod
    def setUpClass(cls):
        data = os.path.realpath(os.path.join(__file__, "..", "data"))
        fna = os.path.join(data, "SRR492066.fna.gz")
        meta_fna = os.path.join(data, "SRR492066.meta.fna.gz")
        meta_faa = os.path.join(data, "SRR492066.meta.faa.gz")

        with gzip.open(fna, "rt") as f:
            cls.record = next(parse(f))
        with gzip.open(meta_faa, "rt") as f:
            cls.proteins = [
                record
                for record in parse(f)
                if record.id.startswith("{}_".format(cls.record.id))
            ]
        with gzip.open(meta_fna, "rt") as f:
            cls.genes = [
                record
                for record in parse(f)
                if record.id.startswith("{}_".format(cls.record.id))
            ]

        cls.preds = cls.find_genes(str(cls.record.seq))

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
