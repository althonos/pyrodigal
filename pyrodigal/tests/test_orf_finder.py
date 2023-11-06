import abc
import gzip
import os
import pickle
import textwrap
import unittest
import warnings

from .. import GeneFinder, MetagenomicBins, METAGENOMIC_BINS
from . import data


class _GeneFinderTestCase(object):
    def assertGeneEqual(self, gene1, gene2):
        self.assertEqual(gene1.begin, gene2.begin)
        self.assertEqual(gene1.end, gene2.end)
        self.assertEqual(gene1.strand, gene2.strand)
        self.assertEqual(gene1.partial_begin, gene2.partial_begin)
        self.assertEqual(gene1.partial_end, gene2.partial_end)
        self.assertEqual(gene1.start_type, gene2.start_type)
        self.assertEqual(gene1.rbs_spacer, gene2.rbs_spacer)
        self.assertEqual(gene1.gc_cont, gene2.gc_cont)
        self.assertEqual(gene1.translation_table, gene2.translation_table)
        self.assertEqual(gene1.cscore, gene2.cscore)
        self.assertEqual(gene1.rscore, gene2.rscore)
        self.assertEqual(gene1.sscore, gene2.sscore)
        self.assertEqual(gene1.tscore, gene2.tscore)
        self.assertEqual(gene1.uscore, gene2.uscore)
        self.assertEqual(gene1.score, gene2.score)

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
            self.assertEqual(gene._gene_data(1), gene_data.strip())

    def assertPredictionsEqual(self, predictions, proteins):
        self.assertTranslationsEqual(predictions, proteins)
        self.assertCoordinatesEqual(predictions, proteins)
        self.assertRbsMotifsEqual(predictions, proteins)
        self.assertStartTypesEqual(predictions, proteins)
        self.assertRbsSpacersEqual(predictions, proteins)
        self.assertGCContentEqual(predictions, proteins)
        self.assertGeneDataEqual(predictions, proteins)


class _TestMode(_GeneFinderTestCase):
    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_find_genes_KK037166(self):
        record = data.load_record("KK037166.fna.gz")
        proteins = data.load_records("KK037166.{}.faa.gz".format(self.mode))
        genes = data.load_records("KK037166.{}.fna.gz".format(self.mode))

        preds = self.find_genes(self.get_sequence(record))
        self.assertGenesEqual(preds, genes)
        self.assertPredictionsEqual(preds, proteins)

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_find_genes_SRR492066(self):
        record = data.load_record("SRR492066.fna.gz")
        proteins = data.load_records("SRR492066.{}.faa.gz".format(self.mode))
        genes = data.load_records("SRR492066.{}.fna.gz".format(self.mode))

        preds = self.find_genes(self.get_sequence(record))
        self.assertGenesEqual(preds, genes)
        self.assertPredictionsEqual(preds, proteins)

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_find_genes_MIIJ01000039(self):
        record = data.load_record("MIIJ01000039.fna.gz")
        proteins = data.load_records("MIIJ01000039.{}.faa.gz".format(self.mode))
        genes = data.load_records("MIIJ01000039.{}.fna.gz".format(self.mode))

        preds = self.find_genes(self.get_sequence(record))
        self.assertGenesEqual(preds, genes)
        self.assertPredictionsEqual(preds, proteins)


class _TestBin(object):
    @classmethod
    def get_sequence(cls, r):
        return r.seq.encode("ascii")


class _TestTxt(object):
    @classmethod
    def get_sequence(cls, r):
        return r.seq


class _TestSingle(object):
    mode = "single"

    @classmethod
    def find_genes(cls, seq):
        p = GeneFinder(meta=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.train(seq)
        return p.find_genes(seq)


class _TestMeta(object):
    mode = "meta"

    @classmethod
    def find_genes(cls, seq):
        p = GeneFinder(meta=True)
        return p.find_genes(seq)


class TestMetaTxt(_TestMeta, _TestTxt, _TestMode, unittest.TestCase):
    pass


class TestMetaBin(_TestMeta, _TestBin, _TestMode, unittest.TestCase):
    pass


class TestSingleTxt(_TestSingle, _TestTxt, _TestMode, unittest.TestCase):
    pass


class TestSingleBin(_TestSingle, _TestBin, _TestMode, unittest.TestCase):
    pass


class TestGeneFinder(_GeneFinderTestCase, unittest.TestCase):
    def test_invalid_overlap(self):
        self.assertRaises(ValueError, GeneFinder, min_gene=10, max_overlap=100)
        self.assertRaises(ValueError, GeneFinder, max_overlap=-1)

    def test_invalid_min_gene(self):
        self.assertRaises(ValueError, GeneFinder, min_gene=-1)


class TestMeta(_GeneFinderTestCase, unittest.TestCase):
    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_train(self):
        record = data.load_record("SRR492066.fna.gz")
        p = GeneFinder(meta=True)
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
        p = GeneFinder(meta=True, closed=False)
        genes = p.find_genes(textwrap.dedent(seq).replace("\n", ""))
        self.assertEqual(len(genes), 1)
        self.assertEqual(genes[0].start_type, "Edge")
        self.assertTrue(genes[0].partial_begin)
        self.assertTrue(genes[0].partial_end)

    def test_short_sequences(self):
        seq = "AATGTAGGAAAAACAGCATTTTCATTTCGCCATTTT"
        p = GeneFinder(meta=True)
        for i in range(1, len(seq)):
            genes = p.find_genes(seq[:i])
            self.assertEqual(len(genes), 0)
            self.assertRaises(StopIteration, next, iter(genes))

    def test_empty_sequence(self):
        p = GeneFinder(meta=True)
        genes = p.find_genes("")
        self.assertEqual(len(genes), 0)
        self.assertRaises(StopIteration, next, iter(genes))

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_find_genes_masked_MIIJ01000039(self):
        record = data.load_record("MIIJ01000039.fna.gz")
        proteins = data.load_records("MIIJ01000039.{}.faa.gz".format("meta+mask"))
        genes = data.load_records("MIIJ01000039.{}.fna.gz".format("meta+mask"))

        orf_finder = GeneFinder(meta=True, mask=True)
        preds = orf_finder.find_genes(record.seq)
        self.assertEqual(len(preds.sequence.masks), 1)
        self.assertGenesEqual(preds, genes)
        self.assertPredictionsEqual(preds, proteins)

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_find_genes_large_minsize(self):
        record = data.load_record("KK037166.fna.gz")
        genes = data.load_records("KK037166.{}.fna.gz".format("meta+mask"))
        large_genes = [gene for gene in genes if len(gene.seq) >= 200]

        orf_finder = GeneFinder(meta=True, min_gene=200, min_edge_gene=200, mask=True)
        preds = orf_finder.find_genes(record.seq)
        self.assertEqual(len(preds), len(large_genes))

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_find_genes_small_minsize(self):
        record = data.load_record("KK037166.fna.gz")
        genes = data.load_records("KK037166.{}.fna.gz".format("meta+mask"))

        orf_finder = GeneFinder(
            meta=True, min_gene=30, min_edge_gene=20, max_overlap=20, mask=True
        )
        preds = orf_finder.find_genes(record.seq)
        self.assertGreaterEqual(len(preds), len(genes))

    def test_find_small_genes_consistency(self):
        # reported in issue #13
        seq = """
        TTCGTCAGTCGTTCTGTTTCATTCAATACGATAGTAATGTATTTTTCGTGCATTTCCGGT
        GGAATCGTGCCGTCCAGCATAGCCTCCAGATATCCCCTTATAGAGGTCAGAGGGGAACGG
        AAATCGTGGGATACATTGGCTACAAACTTTTTCTGATCATCCTCGGAACGGGCAATTTCG
        CTTGCCATATAATTCAGACAGGAAGCCAGATAACCGATTTCATCCTCACTATCGACCTGA
        AATTCATAATGCATATTACCGGCAGCATACTGCTCTGTGGCATGAGTGATCTTCCTCAGA
        GGAATATATACGATCTCAGTGAAAAAGATCAGAATGATCAGGGATAGCAGGAACAGGATT
        GCCAGGGTGATATAGGAAATATTCAGCAGGTTGTTACAGGATTTCTGAATATCATTCATA
        TCAGTATGGATGACTACATAGCCTTTTACCTTGTAGTTGGAGGTAATGGGAGCAAATACA
        GTAAGTACATCCGAATCAAAATTACCGAAGAAATCACCAACAATGTAATAGGAGCCGCTG
        GTTACGGTCGAATCAAAATTCTCAATGACAACCACATTCTCCACATCTAAGGGACTATTG
        GTATCCAGTACCAGTCGTCCGGAGGGATTGATGATGCGAATCTCGGAATTCAGGTAGACC
        GCCAGGGAGTCCAGCTGCATTTTAACGGTCTCCAAAGTTGTTTCACTGGTGTACAATCCG
        CCGGCATAGGTTCCGGCGATCAGGGTTGCTTCGGAATAGAGACTTTCTGCCTTTTCCCGG
        ATCAGATGTTCTTTGGTCATATTGGGAACAAAAGTTGTAACAATGATGAAACCAAATACA
        CCAAAAATAAAATATGCGAGTATAAATTTTAGATAAAGTGTTTTTTTCATAACAAATCCT
        GCTTTTGGTATGACTTAATTACGTACTTCGAATTTATAGCCGATGCCCCAGATGGTGCTG
        ATCTTCCAGTTGGCATGATCCTTGATCTTCTC
        """
        p = GeneFinder(meta=True, closed=True, min_gene=33, max_overlap=0)
        for _ in range(10):
            genes = p.find_genes(textwrap.dedent(seq).replace("\n", ""))
            self.assertEqual(len(genes), 2)
            self.assertEqual(genes[0].start_type, "GTG")
            self.assertEqual(genes[0].begin, 48)
            self.assertEqual(genes[0].end, 347)
            self.assertEqual(genes[1].start_type, "ATG")
            self.assertEqual(genes[1].begin, 426)
            self.assertEqual(genes[1].end, 590)

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_custom_metagenomic_bins(self):
        record = data.load_record("SRR492066.fna.gz")
        
        orf_finder = GeneFinder(meta=True)
        preds = orf_finder.find_genes(record.seq)
        self.assertIsNot(preds.metagenomic_bin, None)
        self.assertEqual(preds.metagenomic_bin.description, METAGENOMIC_BINS[0].description)

        metagenomic_bins = MetagenomicBins(( METAGENOMIC_BINS[0], METAGENOMIC_BINS[4] ))
        orf_finder = GeneFinder(meta=True, metagenomic_bins=metagenomic_bins)
        preds = orf_finder.find_genes(record.seq)
        self.assertIsNot(preds.metagenomic_bin, None)
        self.assertEqual(preds.metagenomic_bin.description, METAGENOMIC_BINS[0].description)

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_empty_metagenomic_bins(self):
        record = data.load_record("SRR492066.fna.gz")
        metagenomic_bins = MetagenomicBins()
        orf_finder = GeneFinder(meta=True, metagenomic_bins=metagenomic_bins)
        preds = orf_finder.find_genes(record.seq)
        self.assertEqual(len(preds), 0)
        self.assertIs(preds.metagenomic_bin, None)
        self.assertIs(preds.training_info, None)


class TestSingle(_GeneFinderTestCase, unittest.TestCase):
    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_train_info(self):
        record = data.load_record("SRR492066.fna.gz")
        p = GeneFinder(meta=False)
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

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_train_not_called(self):
        record = data.load_record("SRR492066.fna.gz")
        p = GeneFinder(meta=False)
        self.assertRaises(RuntimeError, p.find_genes, str(record.seq))

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_training_info_deallocation(self):
        record = data.load_record("SRR492066.fna.gz")
        proteins = data.load_records("SRR492066.{}.faa.gz".format("single"))
        p = GeneFinder(meta=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.train(str(record.seq))
        genes = p.find_genes(str(record.seq))
        del p  # normally should not deallocate training info since it's RC
        self.assertEqual(genes[0].translate(), str(proteins[0].seq))

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_short_sequences(self):
        record = data.load_record("SRR492066.fna.gz")
        seq = "AATGTAGGAAAAACAGCATTTTCATTTCGCCATTTT"
        p = GeneFinder(meta=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.train(str(record.seq[:20000]))
        for i in range(1, len(seq)):
            genes = p.find_genes(seq[:i])
            self.assertEqual(len(genes), 0)
            self.assertRaises(StopIteration, next, iter(genes))

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_empty_sequence(self):
        record = data.load_record("SRR492066.fna.gz")
        p = GeneFinder(meta=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.train(str(record.seq[:20000]))
        genes = p.find_genes("")
        self.assertEqual(len(genes), 0)
        self.assertRaises(StopIteration, next, iter(genes))

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_pickle(self):
        record = data.load_record("SRR492066.fna.gz")
        # train separately
        p1 = GeneFinder(meta=False, min_gene=60)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p1.train(str(record.seq[:20000]))
        # pickle/unpickle the GeneFinder
        p2 = pickle.loads(pickle.dumps(p1))
        # make sure the same genes are found
        g1 = p1.find_genes(record.seq)
        g2 = p2.find_genes(record.seq)
        # make sure genes are the same
        self.assertEqual(len(g1), len(g2))
        for gene1, gene2 in zip(g1, g2):
            self.assertGeneEqual(gene1, gene2)

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_training_info_pickle(self):
        record = data.load_record("SRR492066.fna.gz")
        # train separately
        p1 = GeneFinder(meta=False, min_gene=60)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p1.train(str(record.seq[:20000]))
        # pickle/unpickle the TrainingInfo
        ti = pickle.loads(pickle.dumps(p1.training_info))
        p2 = GeneFinder(meta=False, training_info=ti, min_gene=60)
        # make sure the same genes are found
        g1 = p1.find_genes(record.seq)
        g2 = p2.find_genes(record.seq)
        # make sure genes are the same
        self.assertEqual(len(g1), len(g2))
        for gene1, gene2 in zip(g1, g2):
            self.assertGeneEqual(gene1, gene2)
