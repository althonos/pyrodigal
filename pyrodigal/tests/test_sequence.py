import unittest

from .. import Sequence
from .._pyrodigal import METAGENOMIC_BINS


class TestSequence(unittest.TestCase):

    def test_str(self):
        s = "ATGCNNNNNNNNNNATGCNNNNNNNNTGC"
        seq = Sequence.from_string(s, mask=False)
        self.assertEqual(str(seq), s)

    def test_no_region_masking(self):
        seq = Sequence.from_string("ATGCNNNNNNNNNNATGCNNNNNNNNTGC", mask=False)
        self.assertEqual(len(seq.masks), 0)

    def test_region_masking(self):
        seq = Sequence.from_string("ATGCNNNNNNNNNNATGCNNNNNNNNTGC", mask=True)
        self.assertEqual(len(seq.masks), 2)
        self.assertEqual(seq.masks[0].begin, 4)
        self.assertEqual(seq.masks[0].end, 13)
        self.assertEqual(seq.masks[1].begin, 18)
        self.assertEqual(seq.masks[1].end, 25)

    def test_shine_dalgarno_exact(self):
        tinf = METAGENOMIC_BINS[0].training_info
        seq = Sequence.from_string("AGGAGGTTAGCAAATATG")
        for i in range(10):
            # AGGAGG if i == 0 else AGG if i == 1 else None
            expected = 24 if i == 0 else 13 if i == 3 else 0
            self.assertEqual(seq.shine_dalgarno(i, 15, tinf), expected, i)
        seq = Sequence.from_string("AGGTGGTTAGCAAATATG")
        for i in range(10):
            # AGG if i == 0 else AGG if i == 1 else None
            expected = 6 if i == 0 else 0
            self.assertEqual(seq.shine_dalgarno(i, 15, tinf), expected, i)

    def test_shine_dalgarno_mismatch(self):
        tinf = METAGENOMIC_BINS[0].training_info
        seq = Sequence.from_string("AGGAGGTTAGCAAATATG")
        for i in range(10):
            expected = 0
            self.assertEqual(seq.shine_dalgarno(i, 15, tinf, exact=False), expected, i)
        seq = Sequence.from_string("AGGTGGTTAGCAAATATG")
        for i in range(10):
            # AGGxGG if i == 0 else None
            expected = 19 if i == 0 else 0
            self.assertEqual(seq.shine_dalgarno(i, 15, tinf, exact=False), expected, i)
