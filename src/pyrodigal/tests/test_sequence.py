import unittest
import pickle

from .. import Sequence, METAGENOMIC_BINS


class TestSequence(unittest.TestCase):
    def test_pickle(self):
        s1 = Sequence("ATGCNNNNNNNNNNATGCNNNNNNNNTGC", mask=True)
        s2 = pickle.loads(pickle.dumps(s1))
        self.assertEqual(len(s1), len(s2))
        self.assertEqual(str(s1), str(s2))
        self.assertEqual(s1.gc, s2.gc)
        self.assertEqual(len(s1.masks), len(s2.masks))
        for m1, m2 in zip(s1.masks, s2.masks):
            self.assertEqual(m1.begin, m2.begin)
            self.assertEqual(m1.end, m2.end)

    def test_unknown(self):
        raw_seq = "ATGCNNNNNNNNNNATGCNNNNNNNNTGC"
        s = Sequence(raw_seq)
        self.assertEqual(s.unknown, raw_seq.count("N"))

    def test_gc_known(self):
        raw_seq = "ATGCNNNNNNNNNNATGCNNNNNNNNTGC"
        s = Sequence(raw_seq)
        gc = raw_seq.count("G") + raw_seq.count("C")
        n_known = len(raw_seq) - raw_seq.count("N")
        self.assertEqual(s.gc_known, gc / n_known)

    def test_str(self):
        s = "ATGCNNNNNNNNNNATGCNNNNNNNNTGC"
        seq = Sequence(s, mask=False)
        self.assertEqual(str(seq), s)

    def test_no_region_masking(self):
        seq = Sequence("ATGCNNNNNNNNNNATGCNNNNNNNNTGC", mask=False)
        self.assertEqual(len(seq.masks), 0)

    def test_region_masking(self):
        seq = Sequence("ATGCNNNNNNNNNNATGCNNNNNNNNTGC", mask=True, mask_size=0)
        self.assertEqual(len(seq.masks), 2)
        self.assertEqual(seq.masks[0].begin, 4)
        self.assertEqual(seq.masks[0].end, 14)
        self.assertEqual(seq.masks[1].begin, 18)
        self.assertEqual(seq.masks[1].end, 26)
        seq = Sequence("ATGCNNNNNNNNNNATGCNNNNNNNNTGC", mask=True, mask_size=10)
        self.assertEqual(len(seq.masks), 1)
        self.assertEqual(seq.masks[0].begin, 4)
        self.assertEqual(seq.masks[0].end, 14)

    def test_shine_dalgarno_exact(self):
        tinf = METAGENOMIC_BINS[0].training_info
        seq = Sequence("AGGAGGTTAGCAAATATG")
        for i in range(10):
            # AGGAGG if i == 0 else AGG if i == 1 else None
            expected = 24 if i == 0 else 13 if i == 3 else 0
            self.assertEqual(seq.shine_dalgarno(i, 15, tinf), expected, i)
        seq = Sequence("AGGTGGTTAGCAAATATG")
        for i in range(10):
            # AGG if i == 0 else AGG if i == 1 else None
            expected = 6 if i == 0 else 0
            self.assertEqual(seq.shine_dalgarno(i, 15, tinf), expected, i)

    def test_shine_dalgarno_mismatch(self):
        tinf = METAGENOMIC_BINS[0].training_info
        seq = Sequence("AGGAGGTTAGCAAATATG")
        for i in range(10):
            expected = 0
            self.assertEqual(seq.shine_dalgarno(i, 15, tinf, exact=False), expected, i)
        seq = Sequence("AGGTGGTTAGCAAATATG")
        for i in range(10):
            # AGGxGG if i == 0 else None
            expected = 19 if i == 0 else 0
            self.assertEqual(seq.shine_dalgarno(i, 15, tinf, exact=False), expected, i)

    def test_mask_trailing(self):
        seq = Sequence("AGCGGGCTACTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", mask=True, mask_size=10)
        self.assertEqual(len(seq.masks), 1)