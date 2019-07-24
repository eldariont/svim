import unittest

from svim.SVIM_intra import analyze_cigar_indel

class TestSVIMIntra(unittest.TestCase):

    def test_analyze_cigar_indel(self):
        tuples = [(5,10), (4,20), (0,10), (7,10), (8,5), (0,5), (1,50), (0,30), (4,25), (5,15)]
        indels = [(30, 50, 50, "INS")]
        self.assertEqual(analyze_cigar_indel(tuples, 30), indels)

        tuples = [(5,10), (4,20), (0,30), (2,50), (0,30), (4,25), (5,15)]
        indels = [(30, 50, 50, "DEL")]
        self.assertEqual(analyze_cigar_indel(tuples, 30), indels)

        tuples = [(5,10), (4,20), (0,30), (2,40), (1,50), (0,30), (4,25), (5,15)]
        indels = [(30, 50, 40, "DEL"), (70, 50, 50, "INS")]
        self.assertEqual(analyze_cigar_indel(tuples, 30), indels)

        tuples = [(5,10), (4,20), (0,30), (1,40), (2,50), (0,30), (4,25), (5,15)]
        indels = [(30, 50, 40, "INS"), (30, 90, 50, "DEL")]
        self.assertEqual(analyze_cigar_indel(tuples, 30), indels)

if __name__ == '__main__':
    unittest.main()
