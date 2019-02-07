import unittest

from svim.SVIM_input_parsing import guess_file_type

class TestSVSignature(unittest.TestCase):

    def test_guess_file_type(self):
        fasta_paths = ["/test/path/test.file.fa", "/test/path/test.file.fasta", "/test/path/test.file.FA"]
        fastq_paths = ["/test/path/test.file.fq", "/test/path/test.file.FQ", "/test/path/test.file.fastq"]

        for p in fasta_paths:
            self.assertEqual(guess_file_type(p), "fasta")
            zipped_path1 = p + ".gz"
            zipped_path2 = p + ".gzip"
            self.assertEqual(guess_file_type(zipped_path1), "fasta_gzip")
            self.assertEqual(guess_file_type(zipped_path2), "fasta_gzip")
            list_path = p + ".fn"
            self.assertEqual(guess_file_type(list_path), "list")
        
        for p in fastq_paths:
            self.assertEqual(guess_file_type(p), "fastq")
            zipped_path1 = p + ".gz"
            zipped_path2 = p + ".gzip"
            self.assertEqual(guess_file_type(zipped_path1), "fastq_gzip")
            self.assertEqual(guess_file_type(zipped_path2), "fastq_gzip")
            list_path = p + ".fn"
            self.assertEqual(guess_file_type(list_path), "list")