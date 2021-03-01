import abc
import gzip
import os
import textwrap
import unittest
import warnings

import Bio.SeqIO
from pyrodigal import Pyrodigal


class _TestPyrodigalMode(object):

    mode = None

    @classmethod
    @abc.abstractmethod
    def find_genes(s):
        return NotImplemented

    @classmethod
    def setUpClass(cls):
        data = os.path.realpath(os.path.join(__file__, "..", "data"))
        fna = os.path.join(data, "SRR492066.fna.gz")
        meta_fna = os.path.join(data, "SRR492066.{}.fna.gz".format(cls.mode))
        meta_faa = os.path.join(data, "SRR492066.{}.faa.gz".format(cls.mode))

        with gzip.open(fna, "rt") as f:
            cls.record = next(Bio.SeqIO.parse(f, "fasta"))
        with gzip.open(meta_faa, "rt") as f:
            cls.proteins = [
                record
                for record in Bio.SeqIO.parse(f, "fasta")
                if record.id.startswith("{}_".format(cls.record.id))
            ]
        with gzip.open(meta_fna, "rt") as f:
            cls.genes = [
                record
                for record in Bio.SeqIO.parse(f, "fasta")
                if record.id.startswith("{}_".format(cls.record.id))
            ]

        cls.preds = cls.find_genes(str(cls.record.seq))

    def test_translate(self):
        self.assertEqual(len(self.preds), len(self.proteins))
        for gene, protein in zip(self.preds, self.proteins):
            self.assertEqual(gene.translate(), str(protein.seq))

    def test_coordinates(self):
        self.assertEqual(len(self.preds), len(self.proteins))
        for gene, protein in zip(self.preds, self.proteins):
            id_, start, end, strand, *_ = protein.description.split(" # ")
            self.assertEqual(gene.begin, int(start))
            self.assertEqual(gene.end, int(end))
            self.assertEqual(gene.strand, int(strand))

    def test_rbs_motif(self):
        self.assertEqual(len(self.preds), len(self.proteins))
        for gene, protein in zip(self.preds, self.proteins):
            *_, raw_data = protein.description.split(" # ")
            data = dict(keyval.split("=") for keyval in raw_data.split(";"))
            if data["rbs_motif"] != "None":
                self.assertEqual(gene.rbs_motif, data["rbs_motif"])
            else:
                self.assertIs(gene.rbs_motif, None)

    def test_rbs_spacer(self):
        self.assertEqual(len(self.preds), len(self.proteins))
        for gene, protein in zip(self.preds, self.proteins):
            *_, raw_data = protein.description.split(" # ")
            data = dict(keyval.split("=") for keyval in raw_data.split(";"))
            if data["rbs_spacer"] != "None":
                self.assertEqual(gene.rbs_spacer, data["rbs_spacer"])
            else:
                self.assertIs(gene.rbs_spacer, None)

    def test_start_type(self):
        self.assertEqual(len(self.preds), len(self.proteins))
        for gene, protein in zip(self.preds, self.proteins):
            *_, raw_data = protein.description.split(" # ")
            data = dict(keyval.split("=") for keyval in raw_data.split(";"))
            self.assertEqual(gene.start_type, data["start_type"])


class TestPyrodigalMeta(_TestPyrodigalMode, unittest.TestCase):
    mode = "meta"

    @classmethod
    def find_genes(cls, seq):
        p = Pyrodigal(meta=True)
        return p.find_genes(seq)

    def test_train(self):
        p = Pyrodigal(meta=True)
        self.assertRaises(RuntimeError, p.train, str(self.record.seq))

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
        p = Pyrodigal(meta=True, closed=False)
        genes = p.find_genes(textwrap.dedent(seq).replace("\n", ""))
        self.assertEqual(len(genes), 1)
        self.assertEqual(genes[0].start_type, "Edge")
        self.assertTrue(genes[0].partial_begin)
        self.assertTrue(genes[0].partial_end)

    def test_short_sequence(self):
        seq = "AATGTAGGAAAAACAGCATTTTCATTTCGCCATTTT"
        p = Pyrodigal(meta=True)
        genes = p.find_genes(seq)
        self.assertEqual(len(genes), 0)
        self.assertRaises(StopIteration, next, iter(genes))


class TestPyrodigalSingle(_TestPyrodigalMode, unittest.TestCase):
    mode = "single"

    @classmethod
    def find_genes(cls, seq):
        p = Pyrodigal(meta=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.train(seq)
        return p.find_genes(seq)

    def test_train_not_called(self):
        p = Pyrodigal(meta=False)
        self.assertRaises(RuntimeError, p.find_genes, str(self.record.seq))

    def test_training_info_deallocation(self):
        p = Pyrodigal(meta=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.train(str(self.record.seq))
        genes = p.find_genes(str(self.record.seq))
        del p # normally should not deallocate training info since it's RC
        self.assertEqual(genes[0].translate(), str(self.proteins[0].seq))
