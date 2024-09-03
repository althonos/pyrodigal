import abc
import gzip
import os
import textwrap
import unittest
import warnings

from .. import GeneFinder
from . import data


@unittest.skipUnless(data.files, "importlib.resources not available")
class TestGene(unittest.TestCase):
    @classmethod
    def find_genes(cls, seq):
        p = GeneFinder(meta=True)
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
        for gene, prot in zip(self.preds, self.proteins):
            self.assertEqual(self.preds[0].translation_table, 4)
            self.assertEqual(gene.translate(), prot.seq)
            self.assertEqual(gene.translate(translation_table=4), prot.seq)
            self.assertNotEqual(gene.translate(translation_table=13), prot.seq)
            self.assertEqual(gene.translate(), prot.seq)

    def test_translate_warning(self):
        gene = self.preds[0]
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.assertRaises(
                UserWarning, 
                gene.translate, 
                translation_table=2
            )

    def test_translate_strictness(self):
        finder = GeneFinder(meta=True, closed=True)
        seq = "ATG" + ("CCN" + "CGN" + "CTN" + "ACN" + "GTN" + "GCN" + "GGN" +"TCN") * 10 + "TGA"
        gene = finder.find_genes(seq)[0]
        
        strict = gene.translate(strict=True, translation_table=11)
        non_strict = gene.translate(strict=False, translation_table=11)

        self.assertEqual(strict[1:9], "XXXXXXXX")
        self.assertEqual(non_strict[1:9], "PRLTVAGS")

        strict = gene.translate(strict=True, translation_table=12)
        non_strict = gene.translate(strict=False, translation_table=12)

        self.assertEqual(strict[1:9], "XXXXXXXX")
        self.assertEqual(non_strict[1:9], "PRXTVAGS")

    def test_translate_no_include_stop(self):
        for gene, prot in zip(self.preds, self.proteins):
            self.assertEqual(gene.translate(include_stop=False), prot.seq.rstrip("*"))

    def test_translate_error(self):
        gene = self.preds[0]
        self.assertRaises(ValueError, gene.translate, 7)

    def test_sequence(self):
        for gene in self.preds:
            if gene.strand == 1:
                self.assertEqual(gene.sequence(), self.record.seq[gene.begin-1:gene.end])