import unittest

from .. import Sequence


class TestSequence(unittest.TestCase):

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
