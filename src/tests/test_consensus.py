import unittest
from unittest.mock import Mock
from random import choices

from svim.SVIM_COMBINE import generate_insertion_consensus
from svim.SVSignature import SignatureInsertion, SignatureClusterUniLocal

class TestSVIMConsensus(unittest.TestCase):
    def setUp(self):
        self.nucleotides = ["A", "C", "G", "T"]
        #Prepare reference
        self.genome = "A"*100 + "C"*100
        self.reference = Mock()        
        self.reference.fetch = lambda contig, start, end: self.genome[start:end]
        
    def test_skipping(self):
        insertion_signatures = []
        for i in range(10):
            insertion_signatures.append(SignatureInsertion("chr1", 100, 100100, "suppl", "read"+str(i), "".join(choices(self.nucleotides, k=100000))))
        cluster_long = SignatureClusterUniLocal("chr1", 100, 100100, 10, 10, insertion_signatures, "INS", 0, 0)
        status_code, _ = generate_insertion_consensus(cluster_long, self.reference, maximum_haplotype_length = 10000)
        self.assertEqual(status_code, 1)
    
    def test_simple(self):
        #Prepare cluster
        insertion_sequence = "".join(choices(self.nucleotides, k=100))
        insertion_signatures = []
        for i in range(10):
            insertion_signatures.append(SignatureInsertion("chr1", 100, 200, "cigar", "read"+str(i), insertion_sequence))
        cluster_simple = SignatureClusterUniLocal("chr1", 100, 200, 10, 10, insertion_signatures, "INS", 0, 0)
        status_code, result = generate_insertion_consensus(cluster_simple, self.reference)
        self.assertEqual(status_code, 0)
        realigned_insertion_start, realigned_insertion_size, insertion_consensus = result
        self.assertEqual(realigned_insertion_start, 100)
        self.assertEqual(realigned_insertion_size, 100)
        self.assertEqual(insertion_consensus, insertion_sequence)

if __name__ == '__main__':
    unittest.main()
