import abc
import gzip
import os
import tempfile
import textwrap
import unittest
import pickle
import warnings

from .. import OrfFinder, TrainingInfo
from .._pyrodigal import METAGENOMIC_BINS
from .fasta import parse
from .utils import load_record


class TestTrainingInfo(unittest.TestCase):
    def assertTrainingInfoEqual(self, t1, t2):
        self.assertEqual(t1.translation_table, t2.translation_table)
        self.assertEqual(t1.gc, t2.gc)
        self.assertEqual(t1.bias, t2.bias)
        self.assertEqual(t1.type_weights, t2.type_weights)
        self.assertEqual(t1.uses_sd, t2.uses_sd)
        self.assertEqual(t1.start_weight, t2.start_weight)

    @classmethod
    def setUpClass(cls):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cls.record = load_record("SRR492066")
            cls.training_info = OrfFinder().train(cls.record.seq)

    def test_roundtrip(self):
        try:
            fd, filename = tempfile.mkstemp()
            with os.fdopen(fd, "wb") as dst:
                self.training_info.dump(dst)
            with open(filename, "rb") as src:
                tinf = TrainingInfo.load(src)
        finally:
            if os.path.exists(filename):
                os.remove(filename)
        self.assertTrainingInfoEqual(tinf, self.training_info)

    def test_load_error(self):
        try:
            fd, filename = tempfile.mkstemp()
            with os.fdopen(fd, "wb") as dst:
                dst.write(b"not ok\n")
            with open(filename, "rb") as src:
                self.assertRaises(EOFError, TrainingInfo.load, src)
        finally:
            if os.path.exists(filename):
                os.remove(filename)

    def test_pickle(self):
        t1 = METAGENOMIC_BINS[0].training_info
        t2 = pickle.loads(pickle.dumps(t1))
        self.assertTrainingInfoEqual(t1, t2)
