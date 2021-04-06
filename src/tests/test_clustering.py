import unittest
import tempfile

from random import uniform, choices

from svim.SVIM_clustering import form_partitions, partition_and_cluster
from svim.SVSignature import SignatureDeletion
from svim.SVIM_input_parsing import parse_arguments

class TestSVIMClustering(unittest.TestCase):
    def setUp(self):
        self.signatures = []
        #group 0
        for i in range(10):
            center = 100000 + uniform(-100, 100)
            half_span = 1000 + uniform(-100, 100)
            new_sig = SignatureDeletion("chr1", center - half_span, center + half_span, "cigar", str(i))
            self.signatures.append(new_sig)
        #group 1
        for i in range(10, 20):
            center = 200000 + uniform(-100, 100)
            half_span = 1000 + uniform(-100, 100)
            new_sig = SignatureDeletion("chr1", center - half_span, center + half_span, "cigar", str(i))
            self.signatures.append(new_sig)
        #group 2
        for i in range(20, 30):
            center = 100000 + uniform(-100, 100)
            half_span = 2000 + uniform(-100, 100)
            new_sig = SignatureDeletion("chr1", center - half_span, center + half_span, "cigar", str(i))
            self.signatures.append(new_sig)
        self.temp_genome = tempfile.NamedTemporaryFile(mode="w")
        self.temp_genome.write('>chr1\n')
        self.temp_genome.write("".join(choices(["A", "C", "G", "T"], k=300000)))
        self.options = parse_arguments('1.5.0', ['alignment', 'myworkdir', 'mybamfile', self.temp_genome.name])

    def tearDown(self):
        self.temp_genome.close()

    def test_partitioning(self):
        partitions = form_partitions(self.signatures, 100)
        self.assertEqual(len(partitions), 2)
        for partition in partitions:
            groups_in_partition = set([int(member.read) // 10 for member in partition])
            self.assertTrue(groups_in_partition in [set([0,2]), set([1])])

    def test_partitioning_large_distance(self):
        partitions = form_partitions(self.signatures, 100000)
        self.assertEqual(len(partitions), 1)
        groups_in_partition = set([int(member.read) // 10 for member in partitions[0]])
        self.assertEqual(groups_in_partition, set([0,1,2]))

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
