import unittest

from svim.SVSignature import SignatureDeletion,SignatureInsertion

class TestSVSignature(unittest.TestCase):

    def test_accessors(self):
        deletion = SignatureDeletion("chr1", 100, 300, "cigar", "read1")

        self.assertEqual(deletion.get_source(), ("chr1", 100, 300))
        self.assertEqual(deletion.get_key(), ("del", "chr1", 200))
    
    def test_mean_distance_to(self):
        deletion1 = SignatureDeletion("chr1", 100, 300, "cigar", "read1")
        deletion2 = SignatureDeletion("chr1", 150, 200, "cigar", "read2")
        deletion3 = SignatureDeletion("chr2", 150, 200, "cigar", "read2")
        insertion = SignatureInsertion("chr1", 150, 200, "cigar", "read2")

        self.assertEqual(deletion1.mean_distance_to(deletion2), 25)
        self.assertEqual(deletion1.mean_distance_to(deletion3), float("Inf"))
        self.assertEqual(deletion1.mean_distance_to(insertion), float("Inf"))
    
    def test_as_string(self):
        deletion1 = SignatureDeletion("chr1", 100, 300, "cigar", "read1")

        self.assertEqual(deletion1.as_string(), "chr1\t100\t300\tdel;cigar\tread1")
        self.assertEqual(deletion1.as_string(":"), "chr1:100:300:del;cigar:read1")


if __name__ == '__main__':
    unittest.main()
