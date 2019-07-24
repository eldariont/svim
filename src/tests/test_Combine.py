import unittest
import tempfile

from random import uniform

from svim.SVIM_input_parsing import parse_arguments
from svim.SVSignature import SignatureDeletion
from svim.SVIM_CLUSTER import cluster_sv_signatures, write_signature_clusters_bed, write_signature_clusters_vcf
from svim.SVIM_COMBINE import combine_clusters

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

    def test_combine(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            options = parse_arguments('1.2.0', ['alignment', tmpdirname, 'mybamfile', 'mygenome'])
            signature_clusters = cluster_sv_signatures(self.signatures, options)

            # Write SV signature clusters
            write_signature_clusters_bed(options.working_dir, signature_clusters)
            write_signature_clusters_vcf(options.working_dir, signature_clusters, '1.2.0')

            combine_clusters(signature_clusters, options)

if __name__ == '__main__':
    unittest.main()
