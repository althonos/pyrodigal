import abc
import gzip
import os
import sys
import tempfile
import textwrap
import unittest
import pickle
import platform
import warnings

from .. import GeneFinder, TrainingInfo, METAGENOMIC_BINS
from . import data


class TestTrainingInfo(unittest.TestCase):
    def assertTrainingInfoEqual(self, t1, t2):
        self.assertEqual(t1.translation_table, t2.translation_table)
        self.assertEqual(t1.gc, t2.gc)
        self.assertEqual(t1.bias, t2.bias)
        self.assertEqual(t1.type_weights, t2.type_weights)
        self.assertEqual(t1.uses_sd, t2.uses_sd)
        self.assertEqual(t1.start_weight, t2.start_weight)
        self.assertEqual(t1.upstream_compositions, t2.upstream_compositions)
        self.assertEqual(t1.motif_weights, t2.motif_weights)
        self.assertEqual(t1.rbs_weights, t2.rbs_weights)

    def test_roundtrip(self):
        try:
            tinf = METAGENOMIC_BINS[0].training_info
            fd, filename = tempfile.mkstemp()
            with os.fdopen(fd, "wb") as dst:
                tinf.dump(dst)
            with open(filename, "rb") as src:
                tinf2 = TrainingInfo.load(src)
        finally:
            if os.path.exists(filename):
                os.remove(filename)
        self.assertTrainingInfoEqual(tinf, tinf2)

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

    @unittest.skipUnless(data.files, "importlib.resources not available")
    def test_pickle(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            record = data.load_record("SRR492066.fna.gz")
            t1 = GeneFinder().train(record.seq)
        t2 = pickle.loads(pickle.dumps(t1))
        self.assertTrainingInfoEqual(t1, t2)

    def test_pickle_metagenomic(self):
        t1 = METAGENOMIC_BINS[0].training_info
        t2 = pickle.loads(pickle.dumps(t1))
        self.assertTrainingInfoEqual(t1, t2)

    @unittest.skipUnless(data.files, "importlib.resources not available")
    @unittest.skipUnless(platform.machine() == "x86_64", "Reference training file was created on x86-64")
    @unittest.skipUnless(sys.platform == "linux", "Reference training file was created on Linux")
    def test_train_closed(self):
        records = data.load_records("GCF_001457455.1_NCTC11397_genomic.fna.gz")
        with data.load("GCF_001457455.1_NCTC11397_genomic.tinf_closed.bin", "rb") as f:
            expected = TrainingInfo.load(f)
        orf_finder = GeneFinder(closed=True)
        actual = orf_finder.train(*(str(r.seq) for r in records))
        self.assertTrainingInfoEqual(actual, expected)
