import unittest
import pysam
import os

from svim.SVIM_COLLECT import retrieve_other_alignments

class TestSAExtraction(unittest.TestCase):
    def setUp(self):
        TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'chimeric_read.bam')
        self.samfile = pysam.AlignmentFile(TESTDATA_FILENAME, "rb")
        self.alignments = list(self.samfile.fetch(until_eof=True))

    def test_fetch(self):
        self.assertEqual(len(self.alignments), 4)
    
    def test_satag_length(self):
        primary = self.alignments[0]
        supplementary_alns = retrieve_other_alignments(primary, self.samfile)
        self.assertEqual(len(supplementary_alns), 3)
    
    def test_satag_extraction_complete(self):
        primary = self.alignments[0]
        supplementary_alns = retrieve_other_alignments(primary, self.samfile)
        for index, aln in enumerate(supplementary_alns):
            self.assertEqual(aln.cigarstring, self.alignments[index+1].cigarstring)
            self.assertEqual(aln.reference_id, self.alignments[index+1].reference_id)
            self.assertEqual(aln.reference_start, self.alignments[index+1].reference_start)
            self.assertEqual(aln.reference_end, self.alignments[index+1].reference_end)
            self.assertEqual(aln.flag, self.alignments[index+1].flag)
            self.assertEqual(aln.mapping_quality, self.alignments[index+1].mapping_quality)
            self.assertEqual(aln.query_sequence, self.alignments[index+1].query_sequence)
            self.assertEqual(aln.query_name, self.alignments[index+1].query_name)
            self.assertEqual(aln.query_alignment_start, self.alignments[index+1].query_alignment_start)
            self.assertEqual(aln.query_alignment_end, self.alignments[index+1].query_alignment_end)


class TestSAExtractionError(unittest.TestCase):
    def setUp(self):
        TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'chimeric_read_errors.bam')
        self.samfile = pysam.AlignmentFile(TESTDATA_FILENAME, "rb")
        self.alignments = list(self.samfile.fetch(until_eof=True))
    
    def test_satag_too_many_fields(self):
        primary = self.alignments[0]
        supplementary_alns = retrieve_other_alignments(primary, self.samfile)
        #first SA entry has too many fields
        self.assertEqual(len(supplementary_alns), 2)

    def test_satag_negative_mapq(self):
        primary = self.alignments[1]
        supplementary_alns = retrieve_other_alignments(primary, self.samfile)
        #negative mapping quality is interpreted as 0
        self.assertEqual(len(supplementary_alns), 1)
        self.assertEqual(supplementary_alns[0].mapping_quality, 0)

if __name__ == '__main__':
    unittest.main()