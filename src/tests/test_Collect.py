import unittest
import pysam
import tempfile

from svim.SVIM_COLLECT import bam_iterator, analyze_alignment_file_querysorted, retrieve_other_alignments
from svim.SVIM_input_parsing import parse_arguments
from random import choice, triangular, uniform

class TestCollect(unittest.TestCase):

    def generate_random_sequence(self, length):
        sequence=""
        for i in range(length):
            sequence+=choice("ACGT")
        return sequence

    def generate_random_cigar_string(self, readlength):
        """Generate random cigar string for a read of a given length. Simulate small mismatches and indels but nothing larger than 10bp."""
        softclip_left = round(triangular(0, readlength, min(1000, readlength * 0.5)))
        non_clipped = readlength - softclip_left
        softclip_right = round(triangular(0, non_clipped, min(1000, non_clipped * 0.5)))
        non_clipped = readlength - softclip_left - softclip_right
        sequence = ""
        read_bases_consumed = 0
        while read_bases_consumed < non_clipped:
            #choose next operation
            if len(sequence) == 0 or sequence[-1] == "I" or sequence[-1] == "D":
                next_operation = "M"
                next_length = round(triangular(1, non_clipped - read_bases_consumed, min(30, non_clipped - read_bases_consumed)))
                read_bases_consumed += next_length
            else:
                next_operation = choice("ID")
                if next_operation == "I":
                    next_length = round(triangular(1, min(10, non_clipped - read_bases_consumed), 1))
                    read_bases_consumed += next_length
                else:
                    next_length = round(triangular(1, 10, 1))
            sequence += str(next_length) + next_operation
        return "{0}S{1}{2}S".format(softclip_left, sequence, softclip_right)


    def generate_random_cigar_string_hardclipped(self, readlength):
        """Generate random cigar string for a read of a given length.
        Simulate small mismatches and indels but nothing larger than 10bp. Simulate hard-clipping and return tuple (left-clipped, right-clipped, cigar)"""
        hardclip_left = round(triangular(0, readlength, min(1000, readlength * 0.5)))
        non_clipped = readlength - hardclip_left
        hardclip_right = round(triangular(0, non_clipped, min(1000, non_clipped * 0.5)))
        non_clipped = readlength - hardclip_left - hardclip_right
        sequence = ""
        read_bases_consumed = 0
        while read_bases_consumed < non_clipped:
            #choose next operation
            if len(sequence) == 0 or sequence[-1] == "I" or sequence[-1] == "D":
                next_operation = "M"
                next_length = round(triangular(1, non_clipped - read_bases_consumed, min(30, non_clipped - read_bases_consumed)))
                read_bases_consumed += next_length
            else:
                next_operation = choice("ID")
                if next_operation == "I":
                    next_length = round(triangular(1, min(10, non_clipped - read_bases_consumed), 1))
                    read_bases_consumed += next_length
                else:
                    next_length = round(triangular(1, 10, 1))
            sequence += str(next_length) + next_operation
        return (hardclip_left, hardclip_right, "{0}H{1}{2}H".format(hardclip_left, sequence, hardclip_right))

    def generate_read(self, qname, flag):
        rname = "chr1"
        pos = int(uniform(1,249250620))
        mapq = int(triangular(0, 60, 50))
        length = int(triangular(100, 20000, 15000))
        cigar = self.generate_random_cigar_string(length)
        seq = self.generate_random_sequence(length)

        read_info = (qname, flag, rname, pos, mapq, cigar, "*", 0, 0, seq, "*", "")

        return read_info


    def generate_split_read_with_sa_tags(self, qname, flag):
        length = int(triangular(100, 20000, 15000))
        seq = self.generate_random_sequence(length)

        suppl_rname = "chr1"
        suppl_pos = int(uniform(1,249250620))
        suppl_mapq = int(triangular(0, 60, 50))
        suppl_hardclipped_left, suppl_hardclipped_right, suppl_cigar = self.generate_random_cigar_string_hardclipped(length)

        prim_rname = "chr1"
        prim_pos = int(uniform(1,249250620))
        prim_mapq = int(triangular(0, 60, 50))
        prim_cigar = self.generate_random_cigar_string(length)

        supplementary_read_info = ( qname,
                                    flag + 2048,
                                    suppl_rname,
                                    suppl_pos,
                                    suppl_mapq,
                                    suppl_cigar,
                                    "*",
                                    0,
                                    0,
                                    seq[suppl_hardclipped_left:-suppl_hardclipped_right],
                                    "*",
                                    "SA:Z:{rname},{pos},{strand},{cigar},{mapq},{nm};".format(rname=prim_rname,
                                                                                              pos=prim_pos,
                                                                                              strand=("-" if flag & 16 else "+"),
                                                                                              cigar=prim_cigar,
                                                                                              mapq=prim_mapq,
                                                                                              nm=0))
        primary_read_info = (   qname,
                                flag,
                                prim_rname,
                                prim_pos,
                                prim_mapq,
                                prim_cigar,
                                "*",
                                0,
                                0,
                                seq,
                                "*",
                                "SA:Z:{rname},{pos},{strand},{cigar},{mapq},{nm};".format(rname=suppl_rname,
                                                                                          pos=suppl_pos,
                                                                                          strand=("-" if flag & 16 else "+"),
                                                                                          cigar=suppl_cigar.replace("H", "S"),
                                                                                          mapq=suppl_mapq,
                                                                                          nm=0))
        return (primary_read_info, supplementary_read_info)

    def setUp(self):
        self.bam_file = tempfile.NamedTemporaryFile()
        header = """@HD	VN:1.0	SO:queryname
@SQ	SN:chr1	LN:249250621
@SQ	SN:chr2	LN:243199373
@SQ	SN:chr3	LN:198022430
@SQ	SN:chr4	LN:191154276
@SQ	SN:chr5	LN:180915260
@SQ	SN:chr6	LN:171115067
@SQ	SN:chr7	LN:159138663
@SQ	SN:chr8	LN:146364022
@SQ	SN:chr9	LN:141213431
@SQ	SN:chr10	LN:135534747
@SQ	SN:chr11	LN:135006516
@SQ	SN:chr12	LN:133851895
@SQ	SN:chr13	LN:115169878
@SQ	SN:chr14	LN:107349540
@SQ	SN:chr15	LN:102531392
@SQ	SN:chr16	LN:90354753
@SQ	SN:chr17	LN:81195210
@SQ	SN:chr18	LN:78077248
@SQ	SN:chr19	LN:59128983
@SQ	SN:chr20	LN:63025520
@SQ	SN:chr21	LN:48129895
@SQ	SN:chr22	LN:51304566
@SQ	SN:chrX	LN:155270560
@SQ	SN:chrY	LN:59373566
@SQ	SN:chr6_ssto_hap7	LN:4928567
@SQ	SN:chr6_mcf_hap5	LN:4833398
@SQ	SN:chr6_cox_hap2	LN:4795371
@SQ	SN:chr6_mann_hap4	LN:4683263
@SQ	SN:chr6_apd_hap1	LN:4622290
@SQ	SN:chr6_qbl_hap6	LN:4611984
@SQ	SN:chr6_dbb_hap3	LN:4610396
@SQ	SN:chr17_ctg5_hap1	LN:1680828
@SQ	SN:chr4_ctg9_hap1	LN:590426
@SQ	SN:chr1_gl000192_random	LN:547496
@SQ	SN:chrUn_gl000225	LN:211173
@SQ	SN:chr4_gl000194_random	LN:191469
@SQ	SN:chr4_gl000193_random	LN:189789
@SQ	SN:chr9_gl000200_random	LN:187035
@SQ	SN:chrUn_gl000222	LN:186861
@SQ	SN:chrUn_gl000212	LN:186858
@SQ	SN:chr7_gl000195_random	LN:182896
@SQ	SN:chrUn_gl000223	LN:180455
@SQ	SN:chrUn_gl000224	LN:179693
@SQ	SN:chrUn_gl000219	LN:179198
@SQ	SN:chr17_gl000205_random	LN:174588
@SQ	SN:chrUn_gl000215	LN:172545
@SQ	SN:chrUn_gl000216	LN:172294
@SQ	SN:chrUn_gl000217	LN:172149
@SQ	SN:chr9_gl000199_random	LN:169874
@SQ	SN:chrUn_gl000211	LN:166566
@SQ	SN:chrUn_gl000213	LN:164239
@SQ	SN:chrUn_gl000220	LN:161802
@SQ	SN:chrUn_gl000218	LN:161147
@SQ	SN:chr19_gl000209_random	LN:159169
@SQ	SN:chrUn_gl000221	LN:155397
@SQ	SN:chrUn_gl000214	LN:137718
@SQ	SN:chrUn_gl000228	LN:129120
@SQ	SN:chrUn_gl000227	LN:128374
@SQ	SN:chr1_gl000191_random	LN:106433
@SQ	SN:chr19_gl000208_random	LN:92689
@SQ	SN:chr9_gl000198_random	LN:90085
@SQ	SN:chr17_gl000204_random	LN:81310
@SQ	SN:chrUn_gl000233	LN:45941
@SQ	SN:chrUn_gl000237	LN:45867
@SQ	SN:chrUn_gl000230	LN:43691
@SQ	SN:chrUn_gl000242	LN:43523
@SQ	SN:chrUn_gl000243	LN:43341
@SQ	SN:chrUn_gl000241	LN:42152
@SQ	SN:chrUn_gl000236	LN:41934
@SQ	SN:chrUn_gl000240	LN:41933
@SQ	SN:chr17_gl000206_random	LN:41001
@SQ	SN:chrUn_gl000232	LN:40652
@SQ	SN:chrUn_gl000234	LN:40531
@SQ	SN:chr11_gl000202_random	LN:40103
@SQ	SN:chrUn_gl000238	LN:39939
@SQ	SN:chrUn_gl000244	LN:39929
@SQ	SN:chrUn_gl000248	LN:39786
@SQ	SN:chr8_gl000196_random	LN:38914
@SQ	SN:chrUn_gl000249	LN:38502
@SQ	SN:chrUn_gl000246	LN:38154
@SQ	SN:chr17_gl000203_random	LN:37498
@SQ	SN:chr8_gl000197_random	LN:37175
@SQ	SN:chrUn_gl000245	LN:36651
@SQ	SN:chrUn_gl000247	LN:36422
@SQ	SN:chr9_gl000201_random	LN:36148
@SQ	SN:chrUn_gl000235	LN:34474
@SQ	SN:chrUn_gl000239	LN:33824
@SQ	SN:chr21_gl000210_random	LN:27682
@SQ	SN:chrUn_gl000231	LN:27386
@SQ	SN:chrUn_gl000229	LN:19913
@SQ	SN:chrM	LN:16571
@SQ	SN:chrUn_gl000226	LN:15008
@SQ	SN:chr18_gl000207_random	LN:4262
@PG	ID:ngmlr	PN:nextgenmap-lr	VN:0.2.7	CL:ngmlr -t 10 -r hg19.fa -q reads.fa -o reads.ngmlr.hg19.bam"""
        self.bam_file.write(header.encode('utf-8'))

        self.read_infos = []
        #10 reads with only primary alignment
        for index in range(10):
            read_info = self.generate_read("read{}".format(index+1), 0)
            self.read_infos.append(read_info)
            sam_entry = "\n" + "\t".join([str(el) for el in read_info])
            self.bam_file.write(sam_entry.encode('utf-8'))

        #10 reads with primary and supplementary alignment
        for index in range(10, 20):
            primary_read_info, supplementary_read_info = self.generate_split_read_with_sa_tags("read{}".format(index+1), 0)
            #primary with SA tag
            self.read_infos.append(primary_read_info)
            sam_entry = "\n" + "\t".join([str(el) for el in primary_read_info])
            self.bam_file.write(sam_entry.encode('utf-8'))
            #supplementary
            self.read_infos.append(supplementary_read_info)
            sam_entry = "\n" + "\t".join([str(el) for el in supplementary_read_info])
            self.bam_file.write(sam_entry.encode('utf-8'))            

        self.bam_file.seek(0)
        self.alignment_file = pysam.AlignmentFile(self.bam_file.name, "rb")

    def test_bam_iterator(self):
        bam_it = bam_iterator(self.alignment_file)

        num_primary_only = 0
        num_primary_supplementary = 0
        num_total = 0
        for prim, suppl, sec in bam_it:
            if len(prim) == 1 and len(suppl) == 0 and len(sec) == 0:
                num_primary_only += 1                
            if len(prim) == 1 and len(suppl) == 1 and len(sec) == 0:
                num_primary_supplementary += 1
            num_total += 1
        
        self.assertEqual(num_total, 20)
        self.assertEqual(num_primary_only, 10)
        self.assertEqual(num_primary_supplementary, 10)

    def test_analyze_alignment_file_querysorted(self):
        arguments = ['alignment', 'myworkdir', 'mybamfile', 'mygenome']
        options = parse_arguments('1.2.0', arguments)
        signatures, translocation_signatures_all_bnds = analyze_alignment_file_querysorted(self.alignment_file, options)
        self.assertEqual(len([sig for sig in signatures if sig.signature == "cigar"]), 0)
    
    def test_retrieve_supplementary_alignment_from_primary(self):
        alignment_it = self.alignment_file.fetch(until_eof=True)
        alignments = list(alignment_it)
        for i in range(10,30,2):
            primary = alignments[i]
            supplementary = alignments[i+1]
            retrieved_supplementary_alns = retrieve_other_alignments(primary, self.alignment_file)
            self.assertEqual(len(retrieved_supplementary_alns), 1)
            self.assertEqual(retrieved_supplementary_alns[0].cigarstring, supplementary.cigarstring.replace("H", "S"))
            self.assertEqual(retrieved_supplementary_alns[0].reference_id, supplementary.reference_id)
            self.assertEqual(retrieved_supplementary_alns[0].reference_start, supplementary.reference_start)
            self.assertEqual(retrieved_supplementary_alns[0].reference_end, supplementary.reference_end)
            self.assertEqual(retrieved_supplementary_alns[0].flag, supplementary.flag)
            self.assertEqual(retrieved_supplementary_alns[0].mapping_quality, supplementary.mapping_quality)
            self.assertEqual(retrieved_supplementary_alns[0].query_sequence[retrieved_supplementary_alns[0].query_alignment_start:retrieved_supplementary_alns[0].query_alignment_end], supplementary.query_sequence)
            self.assertEqual(retrieved_supplementary_alns[0].query_name, supplementary.query_name)

    def test_retrieve_primary_alignment_from_supplementary(self):
        alignment_it = self.alignment_file.fetch(until_eof=True)
        alignments = list(alignment_it)
        for i in range(10,30,2):
            primary = alignments[i]
            supplementary = alignments[i+1]
            retrieved_primary_alns = retrieve_other_alignments(supplementary, self.alignment_file)
            self.assertEqual(len(retrieved_primary_alns), 0)