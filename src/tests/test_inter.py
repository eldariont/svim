import unittest

from svim.SVIM_inter import is_similar

class TestSVIMInter(unittest.TestCase):

    def test_is_similar(self):
        self.assertFalse(is_similar("chrI", 0, 100, "chrII", 0, 100))
        self.assertTrue(is_similar("chrI", 0, 100, "chrI", 0, 100))
        self.assertTrue(is_similar("chrI", 0, 100, "chrI", 10, 90))
        self.assertFalse(is_similar("chrI", 0, 100, "chrI", 21, 100))

if __name__ == '__main__':
    unittest.main()
