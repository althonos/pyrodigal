import collections.abc
import csv
import gzip
import io
import os
import platform
import sys
import unittest
import pickle

from .. import GeneFinder, Genes, METAGENOMIC_BINS
from . import data


@unittest.skipUnless(data.files, "importlib.resources not available")
class TestGenes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = data.load_record("SRR492066.fna.gz")
        cls.p = GeneFinder(meta=True)
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

    def test_pickle(self):
        genes = pickle.loads(pickle.dumps(self.genes))
        mb = self.genes.metagenomic_bin
        # gene data should be the same before/after pickling
        for gene1, gene2 in zip(self.genes, genes):
            self.assertEqual(gene1._gene_data(1), gene2._gene_data(1))
        # training info should be preserved
        self.assertIs(genes.training_info, genes.metagenomic_bin.training_info)
        self.assertEqual(genes.metagenomic_bin.description, mb.description)
        self.assertEqual(genes.training_info.gc, mb.training_info.gc)



@unittest.skipUnless(data.files, "importlib.resources not available")
class _TestWrite(object):

    @classmethod
    def setUpClass(cls):
        cls.record = data.load_record("SRR492066.fna.gz")
        cls.p = GeneFinder(meta=True)
        cls.genes = cls.p.find_genes(str(cls.record.seq))

    def _write(self, genes, buffer, **kwargs):
        raise NotImplementedError

    def test_empty(self):
        genes = self.p.find_genes("")
        buffer = io.StringIO()
        self._write(genes, buffer)

    def test_reported_bytes(self):
        buffer = io.StringIO()
        n = self._write(self.genes, buffer)
        self.assertEqual(n, buffer.tell())


class TestWriteGFF(_TestWrite, unittest.TestCase):

    def _write(self, genes, buffer, **kwargs):
        return genes.write_gff(buffer, self.record.id, **kwargs)

    def test_example(self):
        buffer = io.StringIO()
        self._write(self.genes, buffer)
        actual = buffer.getvalue().splitlines()
        expected = data.load_text("SRR492066.meta.gff").splitlines()

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

    def test_header(self):
        buffer = io.StringIO()
        self._write(self.genes, buffer, header=True)
        actual = buffer.getvalue().splitlines()
        self.assertEqual(actual[0], "##gff-version  3")
        buffer = io.StringIO()
        self._write(self.genes, buffer, header=False)
        actual = buffer.getvalue().splitlines()
        self.assertNotEqual(actual[0], "##gff-version  3")


class TestWriteTranslation(_TestWrite, unittest.TestCase):

    def _write(self, genes, buffer, **kwargs):
        return genes.write_translations(buffer, self.record.id, **kwargs)

    def test_example(self):
        buffer = io.StringIO()
        self._write(self.genes, buffer)
        actual = buffer.getvalue()
        expected = data.load_text("SRR492066.meta.faa.gz")
        self.assertEqual(actual, expected)


class TestWriteGenes(_TestWrite, unittest.TestCase):

    def _write(self, genes, buffer, **kwargs):
        return genes.write_genes(buffer, self.record.id, **kwargs)

    def test_example(self):
        buffer = io.StringIO()
        self._write(self.genes, buffer)
        actual = buffer.getvalue()
        expected = data.load_text("SRR492066.meta.fna.gz")
        self.assertEqual(actual, expected)


class TestWriteGenbank(_TestWrite, unittest.TestCase):

    def _write(self, genes, buffer, **kwargs):
        return genes.write_genbank(buffer, self.record.id, **kwargs)


class TestWriteScores(_TestWrite, unittest.TestCase):

    def _write(self, genes, buffer, **kwargs):
        return genes.write_scores(buffer, self.record.id, **kwargs)

    def test_example(self):
        buffer = io.StringIO()
        self._write(self.genes, buffer)
        actual = [
            line.strip()
            for line in buffer.getvalue().splitlines()
            if not line.startswith("#") and line.strip()
        ]
        expected = [
            line.strip()
            for line in data.load_text("SRR492066.meta.tsv").splitlines()
            if not line.startswith("#") and line.strip()
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
