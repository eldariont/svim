import unittest

from random import uniform

from svim.SVIM_clustering import partition_and_cluster
from svim.SVSignature import SignatureDeletion
from svim.SVIM_input_parsing import parse_arguments

class TestSVIMClustering(unittest.TestCase):
    def setUp(self):
        self.signatures = []
        for i in range(10):
            center = 100000 + uniform(-100, 100)
            half_span = 1000 + uniform(-100, 100)
            new_sig = SignatureDeletion("chr1", center - half_span, center + half_span, "cigar", str(i))
            self.signatures.append(new_sig)
        for i in range(10, 20):
            center = 200000 + uniform(-100, 100)
            half_span = 1000 + uniform(-100, 100)
            new_sig = SignatureDeletion("chr1", center - half_span, center + half_span, "cigar", str(i))
            self.signatures.append(new_sig)
        for i in range(20, 30):
            center = 100000 + uniform(-100, 100)
            half_span = 2000 + uniform(-100, 100)
            new_sig = SignatureDeletion("chr1", center - half_span, center + half_span, "cigar", str(i))
            self.signatures.append(new_sig)
        self.options = parse_arguments('1.2.0', ['alignment', 'myworkdir', 'mybamfile', 'mygenome'])

    def test_clustering(self):
        clusters = partition_and_cluster(self.signatures, options=self.options, type="deleted regions")
        self.assertEqual(len(clusters), 3)
        for cluster in clusters:
            self.assertEqual(len(set([int(member.read) // 10 for member in cluster.members])), 1)
    
    def test_scores(self):
        clusters = partition_and_cluster(self.signatures, options=self.options, type="deleted regions")
        for cluster in clusters:
            self.assertGreaterEqual(cluster.score, 10)
            self.assertLessEqual(cluster.score, 10 + 20/8)

if __name__ == '__main__':
    unittest.main()
