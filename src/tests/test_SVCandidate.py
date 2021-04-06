import unittest
import unittest.mock

from svim.SVCandidate import CandidateDeletion, CandidateInversion
from svim.SVSignature import SignatureDeletion

class DeletionCandidateTestCase(unittest.TestCase):
    def setUp(self):
        self.contig = "chr1"
        self.start = 1000
        self.end = 2000
        self.members = [SignatureDeletion(self.contig, self.start, self.end, "cigar", "read1")]
        self.score = 2
        self.std_span = 10.2346
        self.std_pos = 21.3453
        self.deletion = CandidateDeletion(self.contig,
                                          self.start,
                                          self.end,
                                          self.members,
                                          self.score,
                                          self.std_span,
                                          self.std_pos)
        self.distance = 1000
        self.size = 1000
        self.other_contig = 'chr9'
        self.deletion2 = CandidateDeletion(self.contig,
                                           self.end+self.distance,
                                           self.end+self.distance+self.size,
                                           self.members,
                                           self.score,
                                           self.std_span,
                                           self.std_pos)
        self.inversion = CandidateInversion(self.contig,
                                            self.end+self.distance,
                                            self.end+self.distance+self.size,
                                            self.members,
                                            self.score,
                                            self.std_span,
                                            self.std_pos)
        self.deletion3 = CandidateDeletion(self.other_contig,
                                           self.end+self.distance,
                                           self.end+self.distance+self.size,
                                           self.members,
                                           self.score,
                                           self.std_span,
                                           self.std_pos)

    def test_get_key(self):
        self.assertEqual(self.deletion.get_key(), ("DEL", self.contig, self.end),
                         'incorrect key')

    def test_get_source(self):
        self.assertEqual(self.deletion.get_source(), (self.contig, self.start, self.end),
                         'incorrect source')

    def test_get_std_span(self):
        self.assertEqual(self.deletion.get_std_span(), round(self.std_span, 2))
        self.assertEqual(self.deletion.get_std_span(3), round(self.std_span, 3))

    def test_get_std_pos(self):
        self.assertEqual(self.deletion.get_std_pos(), round(self.std_pos, 2))
        self.assertEqual(self.deletion.get_std_pos(3), round(self.std_pos, 3))

    def test_downstream_distance_to(self):
        self.assertEqual(self.deletion.downstream_distance_to(self.deletion2), self.distance,
                         'incorrect distance')
        self.assertEqual(self.deletion2.downstream_distance_to(self.deletion), 0,
                         'distance should be 0 for upstream variant')
        self.assertEqual(self.deletion.downstream_distance_to(self.inversion), float('Inf'),
                         'distance should be Inf for other variant types')
        self.assertEqual(self.deletion.downstream_distance_to(self.deletion3), float('Inf'),
                         'distance should be Inf for variant on other contig')

    def test_get_vcf_entry(self):
        vcf_string_normal = "\t".join(["chr1",
                                       "1000",
                                       "PLACEHOLDERFORID",
                                       "N",
                                       "<DEL>",
                                       "2",
                                       "PASS",
                                       "SVTYPE=DEL;END=2000;SVLEN=-1000;SUPPORT=1;STD_SPAN=10.23;STD_POS=21.35",
                                       "GT:DP:AD",
                                       "./.:.:.,."])
        self.assertEqual(self.deletion.get_vcf_entry(), vcf_string_normal,
                         'incorrect VCF record')

        ref_allele = "ACGTCGGATCGCAT"
        alt_allele = "A"
        reference = unittest.mock.Mock()
        reference.fetch.side_effect = [ref_allele, alt_allele]
        vcf_string_sequence_allele = "\t".join(["chr1",
                                                "1000",
                                                "PLACEHOLDERFORID",
                                                "ACGTCGGATCGCAT",
                                                "A",
                                                "2",
                                                "PASS",
                                                "SVTYPE=DEL;END=2000;SVLEN=-1000;SUPPORT=1;STD_SPAN=10.23;STD_POS=21.35",
                                                "GT:DP:AD",
                                                "./.:.:.,."])
        self.assertEqual(self.deletion.get_vcf_entry(sequence_alleles=True, reference=reference),
                         vcf_string_sequence_allele,
                         'incorrect VCF record with sequence alleles')
