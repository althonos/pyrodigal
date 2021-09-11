import abc
import gzip
import os
import textwrap
import unittest
import warnings

from .. import Pyrodigal
from .fasta import parse


def _load_record(name):
    data = os.path.realpath(os.path.join(__file__, "..", "data"))
    fna = os.path.join(data, "{name}.fna.gz".format(name=name))
    with gzip.open(fna, "rt") as f:
        return next(parse(f))

def _load_proteins(name, mode):
    data = os.path.realpath(os.path.join(__file__, "..", "data"))
    faa = os.path.join(data, "{name}.{mode}.faa.gz".format(name=name, mode=mode))
    with gzip.open(faa, "rt") as f:
        return list(parse(f))

def _load_genes(name, mode):
    data = os.path.realpath(os.path.join(__file__, "..", "data"))
    fna = os.path.join(data, "{name}.{mode}.fna.gz".format(name=name, mode=mode))
    with gzip.open(fna, "rt") as f:
        return list(parse(f))


class _TestPyrodigalMode(object):

    mode = None

    @classmethod
    @abc.abstractmethod
    def find_genes(s):
        return NotImplemented

    def assertTranslationsEqual(self, predictions, proteins):
        self.assertEqual(len(predictions), len(proteins))
        for pred, protein in zip(predictions, proteins):
            self.assertEqual(pred.translate(), str(protein.seq))

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

    def test_find_genes_KK037166(self):
        record = _load_record("KK037166")
        proteins = _load_proteins("KK037166", self.mode)

        preds = self.find_genes(record.seq)
        self.assertTranslationsEqual(preds, proteins)
        self.assertCoordinatesEqual(preds, proteins)
        self.assertRbsMotifsEqual(preds, proteins)
        self.assertStartTypesEqual(preds, proteins)
        self.assertRbsSpacersEqual(preds, proteins)

        preds_bin = self.find_genes(record.seq.encode("ascii"))
        self.assertTranslationsEqual(preds_bin, proteins)
        self.assertCoordinatesEqual(preds_bin, proteins)
        self.assertRbsMotifsEqual(preds_bin, proteins)
        self.assertStartTypesEqual(preds_bin, proteins)
        self.assertRbsSpacersEqual(preds_bin, proteins)

    def test_find_genes_SRR492066(self):
        record = _load_record("SRR492066")
        proteins = _load_proteins("SRR492066", self.mode)

        preds = self.find_genes(record.seq)
        self.assertTranslationsEqual(preds, proteins)
        self.assertCoordinatesEqual(preds, proteins)
        self.assertRbsMotifsEqual(preds, proteins)
        self.assertStartTypesEqual(preds, proteins)
        self.assertRbsSpacersEqual(preds, proteins)

        preds_bin = self.find_genes(record.seq.encode("ascii"))
        self.assertTranslationsEqual(preds_bin, proteins)
        self.assertCoordinatesEqual(preds_bin, proteins)
        self.assertRbsMotifsEqual(preds_bin, proteins)
        self.assertStartTypesEqual(preds_bin, proteins)
        self.assertRbsSpacersEqual(preds_bin, proteins)


class TestPyrodigalMeta(_TestPyrodigalMode, unittest.TestCase):
    mode = "meta"

    @classmethod
    def find_genes(cls, seq):
        p = Pyrodigal(meta=True)
        return p.find_genes(seq)

    def test_train(self):
        record = _load_record("SRR492066")
        p = Pyrodigal(meta=True)
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
        p = Pyrodigal(meta=True, closed=False)
        genes = p.find_genes(textwrap.dedent(seq).replace("\n", ""))
        self.assertEqual(len(genes), 1)
        self.assertEqual(genes[0].start_type, "Edge")
        self.assertTrue(genes[0].partial_begin)
        self.assertTrue(genes[0].partial_end)

    def test_short_sequences(self):
        seq = "AATGTAGGAAAAACAGCATTTTCATTTCGCCATTTT"
        p = Pyrodigal(meta=True)
        for i in range(1, len(seq)):
            genes = p.find_genes(seq[:i])
            self.assertEqual(len(genes), 0)
            self.assertRaises(StopIteration, next, iter(genes))

    def test_empty_sequence(self):
        p = Pyrodigal(meta=True)
        genes = p.find_genes("")
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

    def test_train_info(self):
        record = _load_record("SRR492066")
        p = Pyrodigal(meta=False)
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
        record = _load_record("SRR492066")
        p = Pyrodigal(meta=False)
        self.assertRaises(RuntimeError, p.find_genes, str(record.seq))

    def test_training_info_deallocation(self):
        record = _load_record("SRR492066")
        proteins = _load_proteins("SRR492066", self.mode)
        p = Pyrodigal(meta=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.train(str(record.seq))
        genes = p.find_genes(str(record.seq))
        del p # normally should not deallocate training info since it's RC
        self.assertEqual(genes[0].translate(), str(proteins[0].seq))

    def test_short_sequences(self):
        record = _load_record("SRR492066")
        seq = "AATGTAGGAAAAACAGCATTTTCATTTCGCCATTTT"
        p = Pyrodigal(meta=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.train(str(record.seq[:20000]))
        for i in range(1, len(seq)):
            genes = p.find_genes(seq[:i])
            self.assertEqual(len(genes), 0)
            self.assertRaises(StopIteration, next, iter(genes))

    def test_empty_sequence(self):
        record = _load_record("SRR492066")
        p = Pyrodigal(meta=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p.train(str(record.seq[:20000]))
        genes = p.find_genes("")
        self.assertEqual(len(genes), 0)
        self.assertRaises(StopIteration, next, iter(genes))
