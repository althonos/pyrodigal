import gzip
import os
import unittest

import Bio.SeqIO
from pyrodigal import Pyrodigal



class TestPyrodigalMeta(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        data = os.path.realpath(os.path.join(__file__, "..", "data"))
        fna = os.path.join(data, "SRR492066.fna.gz")
        meta_fna = os.path.join(data, "SRR492066.meta.fna.gz")
        meta_faa = os.path.join(data, "SRR492066.meta.faa.gz")

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

    def test_translate(self):
        p = Pyrodigal(meta=True)
        genes = p.find_genes(str(self.record.seq))
        self.assertEqual(len(genes), len(self.proteins))

        for gene, protein in zip(genes, self.proteins):
            self.assertEqual(gene.translate(), str(protein.seq))
