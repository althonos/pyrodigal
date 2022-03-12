import abc
import gzip
import os
import textwrap
import unittest
import warnings

from .. import OrfFinder
from .utils import load_record, load_proteins, load_genes


class _OrfFinderTestCase(object):

    def assertTranslationsEqual(self, predictions, proteins):
        self.assertEqual(len(predictions), len(proteins))
        for pred, protein in zip(predictions, proteins):
            t = pred.translate()
            self.assertEqual(len(t), len(protein.seq))
            self.assertSequenceEqual(t, str(protein.seq))

    def assertCoordinatesEqual(self, predictions, proteins):
        self.assertEqual(len(predictions), len(proteins))
        for gene, protein in zip(predictions, proteins):
            id_, start, end, strand, *_ = protein.description.split(" # ")
            self.assertEqual(gene.begin, int(start))
            self.assertEqual(gene.end, int(end))
            self.assertEqual(gene.strand, int(strand))

    def assertRbsMotifsEqual(self, predictions, proteins):
        self.assertEqual(len(predictions), len(predictions))
        for gene, protein in zip(predictions, proteins):
            *_, raw_data = protein.description.split(" # ")
            data = dict(keyval.split("=") for keyval in raw_data.split(";"))
            if data["rbs_motif"] != "None":
                self.assertEqual(gene.rbs_motif, data["rbs_motif"])
            else:
                self.assertIs(gene.rbs_motif, None)

    def assertStartTypesEqual(self, predictions, proteins):
        self.assertEqual(len(predictions), len(proteins))
        for gene, protein in zip(predictions, proteins):
            *_, raw_data = protein.description.split(" # ")
            data = dict(keyval.split("=") for keyval in raw_data.split(";"))
            self.assertEqual(gene.start_type, data["start_type"])

    def assertRbsSpacersEqual(self, predictions, proteins):
        self.assertEqual(len(predictions), len(proteins))
        for gene, protein in zip(predictions, proteins):
            *_, raw_data = protein.description.split(" # ")
            data = dict(keyval.split("=") for keyval in raw_data.split(";"))
            if data["rbs_spacer"] != "None":
                self.assertEqual(gene.rbs_spacer, data["rbs_spacer"])
            else:
                self.assertIs(gene.rbs_spacer, None)

    def assertGCContentEqual(self, predictions, proteins):
        self.assertEqual(len(predictions), len(proteins))
        for gene, protein in zip(predictions, proteins):
            *_, raw_data = protein.description.split(" # ")
            data = dict(keyval.split("=") for keyval in raw_data.split(";"))
            self.assertAlmostEqual(gene.gc_cont, float(data["gc_cont"]), places=2)

    def assertGenesEqual(self, predictions, genes):
        self.assertEqual(len(predictions), len(genes))
        for pred, gene in zip(predictions, genes):
            self.assertSequenceEqual(pred.sequence(), str(gene.seq))

    def assertGeneDataEqual(self, predictions, proteins):
        self.assertEqual(len(predictions), len(proteins))
        for gene, protein in zip(predictions, proteins):
            *_, gene_data = protein.description.split(" # ")
            self.assertEqual(gene._gene_data, gene_data.strip())

    def assertPredictionsEqual(self, predictions, proteins):
        self.assertTranslationsEqual(predictions, proteins)
        self.assertCoordinatesEqual(predictions, proteins)
        self.assertRbsMotifsEqual(predictions, proteins)
        self.assertStartTypesEqual(predictions, proteins)
        self.assertRbsSpacersEqual(predictions, proteins)
        self.assertGCContentEqual(predictions, proteins)
        self.assertGeneDataEqual(predictions, proteins)


class _TestMode(_OrfFinderTestCase):

    def test_find_genes_KK037166(self):
        record = load_record("KK037166")
        proteins = load_proteins("KK037166", self.mode)
        genes = load_genes("KK037166", self.mode)

        preds = self.find_genes(self.get_sequence(record))
        self.assertGenesEqual(preds, genes)
        self.assertPredictionsEqual(preds, proteins)

    def test_find_genes_SRR492066(self):
        record = load_record("SRR492066")
        proteins = load_proteins("SRR492066", self.mode)
        genes = load_genes("SRR492066", self.mode)

        preds = self.find_genes(self.get_sequence(record))
        self.assertGenesEqual(preds, genes)
        self.assertPredictionsEqual(preds, proteins)

    def test_find_genes_MIIJ01000039(self):
        record = load_record("MIIJ01000039")
        proteins = load_proteins("MIIJ01000039", self.mode)
        genes = load_genes("MIIJ01000039", self.mode)

        preds = self.find_genes(self.get_sequence(record))
        self.assertGenesEqual(preds, genes)
        self.assertPredictionsEqual(preds, proteins)


class _TestBin(object):
    @classmethod
    def get_sequence(cls, r):
        return r.seq.encode('ascii')


class _TestTxt(object):
    @classmethod
    def get_sequence(cls, r):
        return r.seq


class _TestSingle(object):
    mode = "single"
    @classmethod
    def find_genes(cls, seq):
        p = OrfFinder(meta=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.train(seq)
        return p.find_genes(seq)


class _TestMeta(object):
    mode = "meta"
    @classmethod
    def find_genes(cls, seq):
        p = OrfFinder(meta=True)
        return p.find_genes(seq)


class TestMetaTxt(_TestMeta, _TestTxt, _TestMode, unittest.TestCase):
    pass


class TestMetaBin(_TestMeta, _TestBin, _TestMode, unittest.TestCase):
    pass


class TestSingleTxt(_TestSingle, _TestTxt, _TestMode, unittest.TestCase):
    pass


class TestSingleBin(_TestSingle, _TestBin, _TestMode, unittest.TestCase):
    pass


class TestOrfFinder(_OrfFinderTestCase, unittest.TestCase):

    def test_invalid_overlap(self):
        self.assertRaises(ValueError, OrfFinder, min_gene=10, max_overlap=100)
        self.assertRaises(ValueError, OrfFinder, max_overlap=-1)

    def test_invalid_min_gene(self):
        self.assertRaises(ValueError, OrfFinder, min_gene=-1)


class TestMeta(_OrfFinderTestCase, unittest.TestCase):

    def test_train(self):
        record = load_record("SRR492066")
        p = OrfFinder(meta=True)
        self.assertRaises(RuntimeError, p.train, str(record.seq))

    def test_overflow(self):
        # > 180195.SAMN03785337.LFLS01000089
        seq = """
        AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAAC
        AGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATC
        GCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACG
        GTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTG
        ACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACC
        AACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGC
        AATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGC
        GTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCT
        GGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACC
        ACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTG
        ACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACC
        GTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAAC
        CTGAGC
        """
        p = OrfFinder(meta=True, closed=False)
        genes = p.find_genes(textwrap.dedent(seq).replace("\n", ""))
        self.assertEqual(len(genes), 1)
        self.assertEqual(genes[0].start_type, "Edge")
        self.assertTrue(genes[0].partial_begin)
        self.assertTrue(genes[0].partial_end)

    def test_short_sequences(self):
        seq = "AATGTAGGAAAAACAGCATTTTCATTTCGCCATTTT"
        p = OrfFinder(meta=True)
        for i in range(1, len(seq)):
            genes = p.find_genes(seq[:i])
            self.assertEqual(len(genes), 0)
            self.assertRaises(StopIteration, next, iter(genes))

    def test_empty_sequence(self):
        p = OrfFinder(meta=True)
        genes = p.find_genes("")
        self.assertEqual(len(genes), 0)
        self.assertRaises(StopIteration, next, iter(genes))

    def test_find_genes_masked_MIIJ01000039(self):
        record = load_record("MIIJ01000039")
        proteins = load_proteins("MIIJ01000039", "meta+mask")
        genes = load_genes("MIIJ01000039", "meta+mask")

        orf_finder = OrfFinder(meta=True, mask=True)
        preds = orf_finder.find_genes(record.seq)
        self.assertEqual(len(preds.sequence.masks), 1)
        self.assertGenesEqual(preds, genes)
        self.assertPredictionsEqual(preds, proteins)

    def test_find_genes_large_minsize(self):
        record = load_record("KK037166")
        genes = load_genes("KK037166", "meta+mask")
        large_genes = [gene for gene in genes if len(gene.seq) >= 200]

        orf_finder = OrfFinder(meta=True, min_gene=200, min_edge_gene=200, mask=True)
        preds = orf_finder.find_genes(record.seq)
        self.assertEqual(len(preds), len(large_genes))

    def test_find_genes_small_minsize(self):
        record = load_record("KK037166")
        genes = load_genes("KK037166", "meta+mask")

        orf_finder = OrfFinder(meta=True, min_gene=30, min_edge_gene=20, max_overlap=20, mask=True)
        preds = orf_finder.find_genes(record.seq)
        self.assertGreaterEqual(len(preds), len(genes))


class TestSingle(_OrfFinderTestCase, unittest.TestCase):

    def test_train_info(self):
        record = load_record("SRR492066")
        p = OrfFinder(meta=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            info = p.train(record.seq)

        self.assertEqual(info.translation_table, 11)
        self.assertEqual(info.gc, 0.3010045159434068)
        self.assertEqual(info.start_weight, 4.35)
        self.assertEqual(info.bias[0], 2.6770525781861187)
        self.assertEqual(info.bias[1], 0.17260535063729165)
        self.assertEqual(info.bias[2], 0.1503420711765898)
        self.assertEqual(info.type_weights[0], 0.71796361273324)
        self.assertEqual(info.type_weights[1], -1.3722361344058844)
        self.assertEqual(info.type_weights[2], -2.136731395763296)
        self.assertTrue(info.uses_sd)

    def test_train_not_called(self):
        record = load_record("SRR492066")
        p = OrfFinder(meta=False)
        self.assertRaises(RuntimeError, p.find_genes, str(record.seq))

    def test_training_info_deallocation(self):
        record = load_record("SRR492066")
        proteins = load_proteins("SRR492066", "single")
        p = OrfFinder(meta=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.train(str(record.seq))
        genes = p.find_genes(str(record.seq))
        del p # normally should not deallocate training info since it's RC
        self.assertEqual(genes[0].translate(), str(proteins[0].seq))

    def test_short_sequences(self):
        record = load_record("SRR492066")
        seq = "AATGTAGGAAAAACAGCATTTTCATTTCGCCATTTT"
        p = OrfFinder(meta=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.train(str(record.seq[:20000]))
        for i in range(1, len(seq)):
            genes = p.find_genes(seq[:i])
            self.assertEqual(len(genes), 0)
            self.assertRaises(StopIteration, next, iter(genes))

    def test_empty_sequence(self):
        record = load_record("SRR492066")
        p = OrfFinder(meta=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.train(str(record.seq[:20000]))
        genes = p.find_genes("")
        self.assertEqual(len(genes), 0)
        self.assertRaises(StopIteration, next, iter(genes))
