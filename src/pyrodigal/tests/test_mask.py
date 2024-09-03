import unittest

from .. import Mask


class TestMask(unittest.TestCase):
    def test_repr(self):
        mask = Mask(1, 2)
        self.assertEqual(repr(mask), "<pyrodigal.lib.Mask begin=1 end=2>")

    def test_eq(self):
        m1 = Mask(1, 2)
        m2 = Mask(1, 2)
        self.assertIsNot(m1, m2)
        self.assertEqual(m1, m2)
        self.assertNotEqual(m1, object())
        self.assertNotEqual(m1, None)

    def test_intersect(self):
        self.assertTrue(Mask(2, 4).intersects(1, 3))
        self.assertTrue(Mask(2, 4).intersects(1, 5))
        self.assertTrue(Mask(2, 4).intersects(2, 4))
        self.assertFalse(Mask(2, 4).intersects(4, 6)) # range exclusive