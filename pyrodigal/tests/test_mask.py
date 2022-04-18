import unittest

from .. import Mask


class TestMask(unittest.TestCase):
    def test_repr(self):
        mask = Mask(1, 2)
        self.assertEqual(repr(mask), "<pyrodigal._pyrodigal.Mask begin=1 end=2>")

    def test_eq(self):
        m1 = Mask(1, 2)
        m2 = Mask(1, 2)
        self.assertIsNot(m1, m2)
        self.assertEqual(m1, m2)
        self.assertNotEqual(m1, object())
        self.assertNotEqual(m1, None)
