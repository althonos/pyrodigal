import abc
import gzip
import os
import tempfile
import textwrap
import unittest
import warnings

from .. import OrfFinder, TrainingInfo
from .fasta import parse
from .utils import load_record


class TestTrainingInfo(unittest.TestCase):

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

        self.assertEqual(tinf.translation_table, self.training_info.translation_table)
        self.assertEqual(tinf.gc, self.training_info.gc)
        self.assertEqual(tinf.bias, self.training_info.bias)
        self.assertEqual(tinf.type_weights, self.training_info.type_weights)
        self.assertEqual(tinf.uses_sd, self.training_info.uses_sd)
        self.assertEqual(tinf.start_weight, self.training_info.start_weight)

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
