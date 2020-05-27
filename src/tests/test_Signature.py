import unittest

from svim.SVSignature import SignatureDeletion,SignatureInsertion

class TestSVSignature(unittest.TestCase):

    def test_accessors(self):
        deletion = SignatureDeletion("chr1", 100, 300, "cigar", "read1")

        self.assertEqual(deletion.get_source(), ("chr1", 100, 300))
        self.assertEqual(deletion.get_key(), ("DEL", "chr1", 300))
    
    def test_downstream_distance_to(self):
        deletion1 = SignatureDeletion("chr1", 100, 300, "cigar", "read1")
        deletion2 = SignatureDeletion("chr1", 450, 500, "cigar", "read2")
        deletion3 = SignatureDeletion("chr1", 150, 200, "cigar", "read3")
        deletion4 = SignatureDeletion("chr2", 350, 400, "cigar", "read3")
        insertion = SignatureInsertion("chr1", 150, 200, "cigar", "read2", "ACGTAGTAGCTAGCTTTGCTAGCATTAGCGACTGCTTACGCAGCTCCCTA")

        self.assertEqual(deletion1.downstream_distance_to(deletion2), 150)
        self.assertEqual(deletion1.downstream_distance_to(deletion3), 0)
        self.assertEqual(deletion1.downstream_distance_to(deletion4), float("Inf"))
        self.assertEqual(deletion1.downstream_distance_to(insertion), float("Inf"))
    
    def test_as_string(self):
        deletion1 = SignatureDeletion("chr1", 100, 300, "cigar", "read1")

        self.assertEqual(deletion1.as_string(), "chr1\t100\t300\tDEL;cigar\tread1")
        self.assertEqual(deletion1.as_string(":"), "chr1:100:300:DEL;cigar:read1")


if __name__ == '__main__':
    unittest.main()
