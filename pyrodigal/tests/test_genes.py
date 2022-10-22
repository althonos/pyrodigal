import collections.abc
import csv
import gzip
import io
import os
import platform
import sys
import unittest
import pickle

from .. import OrfFinder, METAGENOMIC_BINS
from .fasta import parse


class TestGenes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data = os.path.realpath(os.path.join(__file__, "..", "data"))
        fna = os.path.join(data, "SRR492066.fna.gz")
        meta_faa = os.path.join(data, "SRR492066.meta.faa.gz")

        with gzip.open(fna, "rt") as f:
            cls.record = next(parse(f))

        cls.p = OrfFinder(meta=True)
        cls.genes = cls.p.find_genes(str(cls.record.seq))

    def test_indexing(self):
        length = len(self.genes)
        self.assertEqual(self.genes[0]._gene_data(1), self.genes[-length]._gene_data(1))
        with self.assertRaises(IndexError):
            self.genes[length]
        with self.assertRaises(IndexError):
            self.genes[-length - 1]

    def test_iter(self):
        for i, gene in zip(range(len(self.genes)), self.genes):
            self.assertEqual(gene._gene_data(1), self.genes[i]._gene_data(1))

    def test_reversed(self):
        for i, gene in zip(range(1, len(self.genes) + 1), reversed(self.genes)):
            self.assertEqual(gene._gene_data(1), self.genes[-i]._gene_data(1))

    def test_bool(self):
        self.assertTrue(bool(self.genes))
        self.assertFalse(bool(self.p.find_genes("TTT")))

    def test_translation_table(self):
        valid = set((*range(1, 7), *range(9, 17), *range(21, 26)))
        for gene in self.genes:
            self.assertIn(gene.translation_table, valid)

    @unittest.skipIf(sys.implementation.name != "cpython", "can panic with PyPy")
    def test_collection_abc_subclass(self):
        self.assertIsInstance(self.genes, collections.abc.Sequence)
        self.assertIsInstance(self.genes, collections.abc.Sized)
        self.assertIsInstance(self.genes, collections.abc.Container)
        self.assertIsInstance(self.genes, collections.abc.Iterable)
        if sys.version_info >= (3, 6):
            self.assertIsInstance(self.genes, collections.abc.Reversible)

    def test_write_scores(self):
        buffer = io.StringIO()
        self.genes.write_scores(buffer, self.record.id)
        actual = [
            line.strip()
            for line in buffer.getvalue().splitlines()
            if not line.startswith("#") and line.strip()
        ]

        data = os.path.realpath(os.path.join(__file__, "..", "data"))
        tsv = os.path.join(data, "SRR492066.meta.tsv")
        with open(tsv) as f:
            expected = [
                line.strip() for line in f if not line.startswith("#") and line.strip()
            ]

        r1 = csv.reader(actual, dialect="excel-tab")
        r2 = csv.reader(expected, dialect="excel-tab")
        for row1, row2 in zip(r1, r2):
            self.assertEqual(row1[0], row2[0], "begin coordinates differ")
            self.assertEqual(row1[1], row2[1], "end coordinates differ")
            self.assertEqual(row1[2], row2[2], "strands differ")
            self.assertEqual(row1[6], row2[6], "codons differ")
            self.assertEqual(row1[7], row2[7], "RBS motifs differ")
            self.assertEqual(row1[8], row2[8], "RBS spacers differ")

        # NOTE(@althonos): This is rather convoluted but unfortunately we cannot
        #                  compare lines directly because of rounding issues that
        #                  can occur between different platforms: for instance,
        #                  depending on the compilation flags, the score
        #                  `-32.2550` may be rounded as `-32.25` or `-32.26`.

    def test_pickle(self):
        genes = pickle.loads(pickle.dumps(self.genes))
        mb = self.genes.training_info.metagenomic_bin
        # gene data should be the same before/after pickling
        for gene1, gene2 in zip(self.genes, genes):
            self.assertEqual(gene1._gene_data(1), gene2._gene_data(1))
        # training info should be preserved, in particular if coming
        # from a metagenomic bin the references should be the same
        self.assertIs(genes.training_info, mb.training_info)
        self.assertIs(genes.training_info.metagenomic_bin, mb)

    def test_write_translations(self):
        buffer = io.StringIO()
        self.genes.write_translations(buffer, self.record.id)
        actual = buffer.getvalue()
        
        data = os.path.realpath(os.path.join(__file__, "..", "data"))
        faa = os.path.join(data, "SRR492066.meta.faa.gz")
        with gzip.open(faa) as f:
            expected = f.read().decode()

        self.assertEqual(actual, expected)

    def test_write_genes(self):
        buffer = io.StringIO()
        self.genes.write_genes(buffer, self.record.id)
        actual = buffer.getvalue()
        
        data = os.path.realpath(os.path.join(__file__, "..", "data"))
        fna = os.path.join(data, "SRR492066.meta.fna.gz")
        with gzip.open(fna) as f:
            expected = f.read().decode()

        self.assertEqual(actual, expected)

    def test_write_gff(self):
        buffer = io.StringIO()
        self.genes.write_gff(buffer, self.record.id)
        actual = buffer.getvalue().splitlines()

        data = os.path.realpath(os.path.join(__file__, "..", "data"))
        gff = os.path.join(data, "SRR492066.meta.gff")
        with open(gff) as f:
            expected = f.read().splitlines()

        self.assertEqual(actual[0], expected[0])
        self.assertEqual(actual[1], expected[1])
        self.assertEqual(actual[2].split(";")[1:], expected[2].split(";")[1:])

        for line_actual, line_expected in zip(actual[3:], expected[3:]):
            row_actual = line_actual.split("\t", maxsplit=8)
            row_expected = line_expected.split("\t", maxsplit=8)
            self.maxDiff = None
            self.assertEqual(row_actual[0], row_expected[0])
            self.assertEqual(row_actual[2:8], row_expected[2:8])
            # NOTE: Don't compare all attributes because:
            #       - IDs are differents, see #18
            #       - some discrepancies in the score data, see #19
            attributes_actual = dict(x.split("=", 1) for x in row_actual[8].split(";") if x)
            attributes_expected = dict(x.split("=", 1) for x in row_expected[8].split(";") if x)
            self.assertEqual(attributes_actual['partial'], attributes_expected['partial'])
            self.assertEqual(attributes_actual['rbs_motif'], attributes_expected['rbs_motif'])
            self.assertEqual(attributes_actual['rbs_spacer'], attributes_expected['rbs_spacer'])
            self.assertEqual(attributes_actual['start_type'], attributes_expected['start_type'])
