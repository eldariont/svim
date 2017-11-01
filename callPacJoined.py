from __future__ import print_function

import sys
import argparse
import os
import re
from subprocess import call, Popen, PIPE

import pysam
from Bio import SeqIO
import networkx as nx

from callPacParams import callPacParams
from kmerCounting import find_svs
from callPacIndels import parse_sam_file, find_indels_in_cigar_tuples


class SVEvidence:
    def __init__(self, contig, start, end, type, evidence, read):
        self.contig = contig
        self.start = start
        self.end = end
        self.type = type
        self.evidence = evidence
        self.read = read

    def as_tuple(self):
        return (self.contig, self.start, self.end, self.type, self.evidence, self.read)


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""""")
    parser.add_argument('fasta', type=argparse.FileType('r'))
    parser.add_argument('genome', default="/scratch/cluster/heller_d/genomes/hg19/hg19.fa", type=str)
    parser.add_argument('temp_dir', type=str, help='temp directory')

    parser.add_argument('--tail_span', type=int, default=1000, help='length of read tails')
    parser.add_argument('--tail_min_mapq', type=int, default=30, help='minimum mapping quality')
    parser.add_argument('--tail_min_deviation', type=float, default=-0.1, help='minimum deviation')
    parser.add_argument('--tail_max_deviation', type=float, default=0.2, help='maximum deviation')
    parser.add_argument('--count_win_size', type=int, default=50, help='window size for k-mer counting')
    parser.add_argument('--count_k', type=int, default=7, help='k for k-mer counting')
    parser.add_argument('--count_band', type=float, default=0.5, help='band width')
    parser.add_argument('--stretch_threshold', type=int, default=7, help='z-score threshold')
    parser.add_argument('--stretch_tolerance', type=int, default=2, help='tolerance for stretch finding')
    parser.add_argument('--stretch_min_length', type=int, default=3, help='minimum stretch length')
    parser.add_argument('--path_constant_gap_cost', type=int, default=0, help='constant gap cost for path finding')
    parser.add_argument('--path_convex_gap_cost', type=int, default=3, help='convex gap cost for path finding')
    parser.add_argument('--path_tolerance', type=int, default=2, help='tolerance for overlapping segments')
    parser.add_argument('--align_costs_match', type=int, default=3, help='match cost for alignment')
    parser.add_argument('--align_costs_mismatch', type=int, default=-12, help='mismatch cost for alignment')
    parser.add_argument('--align_costs_gap', type=int, default=-12, help='gap cost for alignment')
    return parser.parse_args()


def create_temp_files(temp_dir, fasta, span):
    """Create temporary FASTA files if they do not exist."""
    if not os.path.exists(temp_dir):
        print("ERROR: Given temp directory does not exist", file=sys.stderr)
        sys.exit()

    if not os.path.exists(temp_dir + '/left.fa'):
        left_file = open(temp_dir + '/left.fa', 'w')
        write_left = True
    else:
        write_left = False
        print("WARNING: Temp file for left sequences exists. Skip", file=sys.stderr)

    if not os.path.exists(temp_dir + '/right.fa'):
        right_file = open(temp_dir + '/right.fa', 'w')
        write_right = True
    else:
        write_right = False
        print("WARNING: Temp file for right sequences exists. Skip", file=sys.stderr)

    if write_left or write_right:
        for line in fasta:
            if line.startswith('>'):
                read_name = line.strip()[1:]
            else:
                sequence = line.strip()
                if len(sequence) < 2000:
                    continue

                if write_left:
                    prefix = sequence[:span]
                    print(">" + read_name, file=left_file)
                    print(prefix, file=left_file)

                if write_right:
                    suffix = sequence[-span:]
                    print(">" + read_name, file=right_file)
                    print(suffix, file=right_file)

    fasta.close()
    if write_left:
        left_file.close()
    if write_right:
        right_file.close()

    print("INFO: Temporary files written", file=sys.stderr)


def run_alignments(temp_dir, genome, fasta):
    """Align full reads and read tails with NGM-LR and BWA MEM, respectively."""
    if not os.path.exists(temp_dir + '/left_aln.rsorted.bam'):
        bwa = Popen(['/scratch/ngsvin/bin/bwa.kit/bwa',
                     'mem', '-x', 'pacbio', '-t', '30', genome, temp_dir + '/left.fa'], stdout=PIPE)
        view = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'view', '-b', '-@', '10'], stdin=bwa.stdout, stdout=PIPE)
        sort = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'sort', '-@', '10', '-n', '-o', temp_dir + '/left_aln.rsorted.bam'], stdin=view.stdout)
        sort.wait()
    else:
        print("WARNING: Alignment for left sequences exists. Skip", file=sys.stderr)

    if not os.path.exists(temp_dir + '/right_aln.rsorted.bam'):
        bwa = Popen(['/scratch/ngsvin/bin/bwa.kit/bwa',
                     'mem', '-x', 'pacbio', '-t', '30', genome, temp_dir + '/right.fa'], stdout=PIPE)
        view = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'view', '-b', '-@', '10'], stdin=bwa.stdout, stdout=PIPE)
        sort = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'sort', '-@', '10', '-n', '-o', temp_dir + '/right_aln.rsorted.bam'], stdin=view.stdout)
        sort.wait()
    else:
        print("WARNING: Alignment for right sequences exists. Skip", file=sys.stderr)

    # Align full reads with NGM-LR
    if not os.path.exists(temp_dir + '/full_aln.chained.rsorted.bam'):
        ngmlr = Popen(['/home/heller_d/bin/miniconda2/bin/ngmlr',
                       '-t', '30', '-r', genome, '-q', os.path.realpath(fasta.name), ], stdout=PIPE)
        view = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'view', '-b', '-@', '10'], stdin=ngmlr.stdout, stdout=PIPE)
        sort = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'sort', '-n', '-@', '10', '-o', temp_dir + '/full_aln.rsorted.bam'],
                     stdin=view.stdout)
        sort.wait()
        if call(['python', '/home/heller_d/bin/bamChain', temp_dir + '/full_aln.rsorted.bam',
                 temp_dir + '/full_aln.chained.bam', '--minmapq', '40']) != 0:
            print("ERROR: Calling bamchain on full sequences failed", file=sys.stderr)
        if call(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools', 'sort', '-n', '-@', '10',
                 temp_dir + '/full_aln.chained.bam', '-o',
                 temp_dir + '/full_aln.chained.rsorted.bam']) != 0:
            print("ERROR: Calling samtools sort on full sequences failed", file=sys.stderr)
    else:
        print("WARNING: Alignment for full sequences exists. Skip", file=sys.stderr)

    print("INFO: Alignment finished", file=sys.stderr)


def natural_representation(qname): 
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', qname)]


def bam_iterator(bam):
    """Returns an iterator for the given BAM file (must be query-sorted). 
    In each call, the alignments of a single read are yielded as a 4-tuple: (read_name, list of primary pysam.AlignedSegment, list of supplementary pysam.AlignedSegment, list of secondary pysam.AlignedSegment)."""
    alignments = bam.fetch(until_eof=True)
    current_aln = alignments.next()
    current_read_name = current_aln.query_name
    current_prim = []
    current_suppl = []
    current_sec = []
    if current_aln.is_secondary:
        current_sec.append(current_aln)
    elif current_aln.is_supplementary:
        current_suppl.append(current_aln)
    else:
        current_prim.append(current_aln)
    while True:
        try:
            next_aln = alignments.next()
            next_read_name = next_aln.query_name
            if next_read_name != current_read_name:
                yield (current_read_name, current_prim, current_suppl, current_sec)
                current_read_name = next_read_name
                current_prim = []
                current_suppl = []
                current_sec = []
            if next_aln.is_secondary:
                current_sec.append(next_aln)
            elif next_aln.is_supplementary:
                current_suppl.append(next_aln)
            else:
                current_prim.append(next_aln)
        except StopIteration:
            break
    yield (current_read_name, current_prim, current_suppl, current_sec)


def check_indel_candidate_minus(left_tail, right_tail, contig, full_read, reference, parameters):
    left_ref_start = left_tail.reference_start
    left_q_start = left_tail.query_alignment_start
    right_ref_end = right_tail.reference_end
    right_q_end = right_tail.query_alignment_end

    read_snippet = str(full_read[parameters.tail_span - left_q_start : len(full_read) - right_q_end].upper())
    ref_snippet = str(reference[contig].seq[right_ref_end:left_ref_start].upper().reverse_complement())
    sv_results = find_svs(ref_snippet, read_snippet, parameters, debug = False)

    sv_evidences = []
    for typ, start, end in sv_results:
        if typ == "del":
            print("Deletion detected: {0}:{1}-{2} (length {3})".format(contig, left_ref_start - end , left_ref_start - start, end - start), file=sys.stdout)
            sv_evidences.append(SVEvidence(contig, left_ref_start - end, left_ref_start - start, typ, "kmer", left_tail.query_name))
        if typ == "ins":
            print("Insertion detected: {0}:{1}-{2} (length {3})".format(contig, left_ref_start - start, left_ref_start - start + (end - start), end - start), file=sys.stdout)
            sv_evidences.append(SVEvidence(contig, left_ref_start - start, left_ref_start - start + (end - start), typ, "kmer", left_tail.query_name))
    return sv_evidences


def check_indel_candidate_plus(left_tail, right_tail, contig, full_read, reference, parameters):
    left_ref_end = left_tail.reference_end
    left_q_end = left_tail.query_alignment_end
    right_ref_start = right_tail.reference_start
    right_q_start = right_tail.query_alignment_start
    
    read_snippet = str(full_read[left_q_end:len(full_read) - parameters.tail_span + right_q_start].upper())
    ref_snippet = str(reference[contig].seq[left_ref_end:right_ref_start].upper())
    sv_results = find_svs(ref_snippet, read_snippet, parameters, debug = False)

    sv_evidences = []
    for typ, start, end in sv_results:
        if typ == "del":
            print("Deletion detected: {0}:{1}-{2} (length {3})".format(contig, left_ref_end + start, left_ref_end + end, end - start), file=sys.stdout)
            sv_evidences.append(SVEvidence(contig, left_ref_end + start, left_ref_end + end, typ, "kmer", left_tail.query_name))
        if typ == "ins":
            print("Insertion detected: {0}:{1}-{2} (length {3})".format(contig, left_ref_end + start, left_ref_end + end, end - start), file=sys.stdout)
            sv_evidences.append(SVEvidence(contig, left_ref_end + start, left_ref_end + end, typ, "kmer", left_tail.query_name))
    return sv_evidences


def check_inv_1(left_tail, right_tail, contig, full_read, reference, parameters):
    left_ref_end = left_tail.reference_end
    left_q_end = left_tail.query_alignment_end
    right_ref_end = right_tail.reference_end
    right_q_end = right_tail.query_alignment_start
    
    read_snippet = str(full_read[left_q_end:len(full_read) - right_q_end].upper())
    ref_snippet_1 = str(reference[contig].seq[left_ref_end:left_ref_end+len(read_snippet)].upper())
    ref_snippet_2 = str(reference[contig].seq[right_ref_end:right_ref_end+len(read_snippet)].upper().reverse_complement())
    sv_results = find_svs(ref_snippet_1 + ref_snippet_2, read_snippet, parameters, debug = False)

    sv_evidences = []
    for typ, start, end in sv_results:
        if typ == "del":
            inv_start = left_ref_end + start
            inv_end = right_ref_end + (len(ref_snippet_1 + ref_snippet_2) - end)
            print("Inversion detected: {0}:{1}-{2} (length {3})".format(contig, inv_start, inv_end, inv_end - inv_start), file=sys.stdout)
            sv_evidences.append(SVEvidence(contig, inv_start, inv_end, "inv", "kmer", left_tail.query_name))
    return sv_evidences


def check_inv_2(left_tail, right_tail, contig, full_read, reference, parameters):
    left_ref_start = left_tail.reference_start
    left_q_start = left_tail.query_alignment_start
    right_ref_start = right_tail.reference_start
    right_q_start = right_tail.query_alignment_start
    
    read_snippet = str(full_read[parameters.tail_span - left_q_start : len(full_read) - parameters.tail_span + right_q_start].upper())
    ref_snippet_1 = str(reference[contig].seq[left_ref_start - len(read_snippet):left_ref_start].upper().reverse_complement())
    ref_snippet_2 = str(reference[contig].seq[right_ref_start - len(read_snippet):right_ref_start].upper())
    sv_results = find_svs(ref_snippet_1 + ref_snippet_2, read_snippet, parameters, debug = False)

    sv_evidences = []
    for typ, start, end in sv_results:
        if typ == "del":
            inv_start = left_ref_start - start
            inv_end = right_ref_start - (len(ref_snippet_1 + ref_snippet_2) - end)
            print("Inversion detected: {0}:{1}-{2} (length {3})".format(contig, inv_start, inv_end, inv_end - inv_start), file=sys.stdout)
            sv_evidences.append(SVEvidence(contig, inv_start, inv_end, "inv", "kmer", left_tail.query_name))
    return sv_evidences


def check_inv_3(left_tail, right_tail, contig, full_read, reference, parameters):
    left_ref_end = left_tail.reference_end
    left_q_end = left_tail.query_alignment_end
    right_ref_end = right_tail.reference_end
    right_q_end = right_tail.query_alignment_start
    
    read_snippet = str(full_read[left_q_end:len(full_read) - right_q_end].upper())
    ref_snippet_1 = str(reference[contig].seq[left_ref_end:left_ref_end+len(read_snippet)].upper())
    ref_snippet_2 = str(reference[contig].seq[right_ref_end:right_ref_end+len(read_snippet)].upper().reverse_complement())
    sv_results = find_svs(ref_snippet_1 + ref_snippet_2, read_snippet, parameters, debug = False)

    sv_evidences = []
    for typ, start, end in sv_results:
        if typ == "del":
            inv_start = right_ref_end + (len(ref_snippet_1 + ref_snippet_2) - end)
            inv_end = left_ref_end + start
            print("Inversion detected: {0}:{1}-{2} (length {3})".format(contig, inv_start, inv_end, inv_end - inv_start), file=sys.stdout)
            sv_evidences.append(SVEvidence(contig, inv_start, inv_end, "inv", "kmer", left_tail.query_name))
    return sv_evidences


def check_inv_4(left_tail, right_tail, contig, full_read, reference, parameters):
    left_ref_start = left_tail.reference_start
    left_q_start = left_tail.query_alignment_start
    right_ref_start = right_tail.reference_start
    right_q_end = right_tail.query_alignment_end
    
    read_snippet = str(full_read[parameters.tail_span - left_q_start : len(full_read) - parameters.tail_span + right_q_start].upper())
    ref_snippet_1 = str(reference[contig].seq[left_ref_start - len(read_snippet):left_ref_start].upper().reverse_complement())
    ref_snippet_2 = str(reference[contig].seq[right_ref_start - len(read_snippet):right_ref_start].upper())
    sv_results = find_svs(ref_snippet_1 + ref_snippet_2, read_snippet, parameters, debug = False)

    sv_evidences = []
    for typ, start, end in sv_results:
        if typ == "del":
            inv_start = right_ref_start - (len(ref_snippet_1 + ref_snippet_2) - end)
            inv_end = left_ref_start - start
            print("Inversion detected: {0}:{1}-{2} (length {3})".format(contig, inv_start, inv_end, inv_end - inv_start), file=sys.stdout)
            sv_evidences.append(SVEvidence(contig, inv_start, inv_end, "inv", "kmer", left_tail.query_name))
    return sv_evidences


def analyze_pair_of_read_tails(left_iterator_object, right_iterator_object, left_bam, right_bam, reads, reference, parameters):
    left_read_name, left_prim, left_suppl, left_sec = left_iterator_object
    right_read_name, right_prim, right_suppl, right_sec = right_iterator_object

    if len(left_prim) != 1 or left_prim[0].is_unmapped or left_prim[0].mapping_quality < parameters.tail_min_mapq:
        return None
    if len(right_prim) != 1 or right_prim[0].is_unmapped or right_prim[0].mapping_quality < parameters.tail_min_mapq:
        return None

    left_ref_chr = left_bam.getrname(left_prim[0].reference_id)
    left_ref_start = left_prim[0].reference_start
    left_ref_end = left_prim[0].reference_end
    left_q_start = left_prim[0].query_alignment_start
    left_q_end = left_prim[0].query_alignment_end

    right_ref_chr = right_bam.getrname(right_prim[0].reference_id)
    right_ref_start = right_prim[0].reference_start
    right_ref_end = right_prim[0].reference_end
    right_q_start = right_prim[0].query_alignment_start
    right_q_end = right_prim[0].query_alignment_end

    full_read = reads[left_read_name].seq
    read_length = len(full_read)

    if left_ref_chr == right_ref_chr:
        if left_prim[0].is_reverse and right_prim[0].is_reverse:
            reference_dist = left_ref_start - right_ref_end
            if reference_dist > 0:
                individual_dist = read_length - right_q_end - (parameters.tail_span - left_q_start)
                percent_shift = (individual_dist - reference_dist) / float(read_length)
                if percent_shift > parameters.tail_max_deviation or percent_shift < parameters.tail_min_deviation:
                    size_estimate = individual_dist - reference_dist - (0.04 * read_length)
                    if size_estimate > -10000:
                        #INDEL candidate, check with k-mer counting
                        #print("INDEL detected between {0} and {1} on {2} (estimated size: {3} bps, read {4})".format(right_ref_end, left_ref_start, left_ref_chr, size_estimate, left_read_name), file=sys.stdout)
                        return check_indel_candidate_minus(left_prim[0], right_prim[0], left_ref_chr, full_read, reference, parameters)
                    else:
                        #Either very large INS or TRANS
                        pass 
            else:
                #TRANS candidate
                pass
        elif not left_prim[0].is_reverse and not right_prim[0].is_reverse:
            reference_dist = right_ref_start - left_ref_end
            if reference_dist > 0:
                individual_dist = read_length - left_q_end - (parameters.tail_span - right_q_start)
                percent_shift = (individual_dist - reference_dist) / float(read_length)
                if percent_shift > parameters.tail_max_deviation or percent_shift < parameters.tail_min_deviation:
                    size_estimate = individual_dist - reference_dist - (0.04 * read_length)
                    if size_estimate > -10000:
                        #INDEL candidate, check with k-mer counting
                        #print("INDEL detected between {0} and {1} on {2} (estimated size: {3} bps, read {4})".format(left_ref_end, right_ref_start, left_ref_chr, size_estimate, left_read_name), file=sys.stdout)
                        return check_indel_candidate_plus(left_prim[0], right_prim[0], left_ref_chr, full_read, reference, parameters)
                    else:
                        #Either very large INS or TRANS
                        pass
            else:
                #TRANS candidate
                pass
        elif not left_prim[0].is_reverse and right_prim[0].is_reverse:
            reference_dist = right_ref_start - left_ref_end
            if reference_dist > 0:
                #INV candidate, right tail in inverted region
                print("INV detected between {0} and {1} on {2} (1, read {3})".format(left_ref_end, right_ref_start, left_ref_chr, left_read_name), file=sys.stdout)
                return check_inv_1(left_prim[0], right_prim[0], left_ref_chr, full_read, reference, parameters)
            else:
                #INV candidate, left tail in inverted region
                print("INV detected between {0} and {1} on {2} (3, read {3})".format(right_ref_end, left_ref_start, left_ref_chr, left_read_name), file=sys.stdout)
                return check_inv_3(left_prim[0], right_prim[0], left_ref_chr, full_read, reference, parameters)
        elif left_prim[0].is_reverse and not right_prim[0].is_reverse:
            reference_dist = right_ref_start - left_ref_end
            if reference_dist > 0:
                #INV candidate, left tail in inverted region
                print("INV detected between {0} and {1} on {2} (2, read {3})".format(left_ref_end, right_ref_start, left_ref_chr, left_read_name), file=sys.stdout)
                return check_inv_2(left_prim[0], right_prim[0], left_ref_chr, full_read, reference, parameters)
            else:
                #INV candidate, right tail in inverted region
                print("INV detected between {0} and {1} on {2} (4, read {3})".format(right_ref_end, left_ref_start, left_ref_chr, left_read_name), file=sys.stdout)
                return check_inv_4(left_prim[0], right_prim[0], left_ref_chr, full_read, reference, parameters)
    else:
        #TRANS candidate
        pass
    return None


def analyze_read_tails(temp_dir, genome, fasta, parameters):
    left_bam = pysam.AlignmentFile(temp_dir + '/left_aln.rsorted.bam')
    right_bam = pysam.AlignmentFile(temp_dir + '/right_aln.rsorted.bam')
    left_it = bam_iterator(left_bam)
    right_it = bam_iterator(right_bam)

    reads = SeqIO.index(fasta.name, "fasta")
    reference = SeqIO.index(genome, "fasta")
    print("INFO: Indexing reads and reference finished", file=sys.stderr)

    sv_evidences = []

    left_iterator_object = left_it.next()
    while True:
        try:
            right_iterator_object = right_it.next()             
            while natural_representation(left_iterator_object[0]) < natural_representation(right_iterator_object[0]):
                left_iterator_object = left_it.next()
            if left_iterator_object[0] == right_iterator_object[0]:
                if int(left_iterator_object[0].split("_")[1]) % 1000 == 0:
                    print("INFO: Processed read", left_iterator_object[0].split("_")[1], file=sys.stderr)
                result = analyze_pair_of_read_tails(left_iterator_object, right_iterator_object, left_bam, right_bam, reads, reference, parameters)
                if result != None:
                    sv_evidences.extend(result)
                left_iterator_object = left_it.next()
        except StopIteration:
            break
    return sv_evidences


def search_svs(temp_dir, genome, fasta, parameters):
    """Search for SVs using aligned read tails and the read and reference regions in-between."""
    left_it = bam_iterator(temp_dir + '/left_aln.rsorted.bam')
    right_it = bam_iterator(temp_dir + '/right_aln.rsorted.bam')
    full_bam = pysam.AlignmentFile(temp_dir + '/full_aln.chained.rsorted.bam')

    # reference = SeqIO.index(genome, "fasta")
    # reads = SeqIO.index(fasta.name, "fasta")
    # print("INFO: Opening output and reference finished", file=sys.stderr)

    sv_evidences = []

    try:
        for left_contig in left_contigs:
            print("INFO: Searching for SVs in", left_contig, "..", file=sys.stderr)
            left_aln_dict = parse_sam_file(left_bam, left_contig)
            right_aln_dict = parse_sam_file(right_bam, left_contig)
            full_aln_dict = parse_sam_file(full_bam, left_contig)

            for read in left_aln_dict.keys():
                read_id = read.split("_")[1]
                original_length = len(reads[read].seq)

                left_aln = left_aln_dict[read][0]
                try:
                    right_aln = right_aln_dict[read][0]
                except IndexError:
                    # print("Right tail is not mapped to same chromosome", file=sys.stderr)
                    continue

                if left_aln.is_unmapped or right_aln.is_unmapped:
                    # print("One or both of the tails is unmapped", file=sys.stderr)
                    continue

                if left_aln.mapping_quality < parameters.tail_min_mapq or right_aln.mapping_quality < parameters.tail_min_mapq:
                    # print("One or both of the tails are mapped with low quality", file=sys.stderr)
                    continue

                left_ref_start = left_aln.reference_start
                left_ref_end = left_aln.reference_end
                left_q_start = left_aln.query_alignment_start
                left_q_end = left_aln.query_alignment_end

                right_ref_start = right_aln.reference_start
                right_ref_end = right_aln.reference_end
                right_q_start = right_aln.query_alignment_start
                right_q_end = right_aln.query_alignment_end

                if left_aln.is_reverse and right_aln.is_reverse:
                    reference_dist = left_ref_start - right_ref_end
                    if reference_dist > 0:
                        individual_dist = original_length - right_q_end - (parameters.tail_span - left_q_start)
                        percent_shift = (individual_dist - reference_dist) / float(original_length)

                        if percent_shift > parameters.tail_max_deviation:
                            insertion_size_estimate = individual_dist - reference_dist - (0.04 * original_length)
                            print("Insertion detected between {0} and {1} on {2} (estimated size: {3} bps, read {4})".format(right_ref_end, left_ref_start, left_contig, insertion_size_estimate, read_id), file=sys.stdout)
                            if insertion_size_estimate < 10000:
                                read_snippet = str(reads[read].seq[parameters.tail_span - left_q_start : original_length - right_q_end].upper())
                                ref_snippet = str(reference[left_contig].seq[right_ref_end:left_ref_start].upper().reverse_complement())
                                sv_results = find_svs(ref_snippet, read_snippet, parameters, debug = False)
                                for typ, start, end in sv_results:
                                    if typ == "del":
                                        print("Deletion detected: {0}: {1} - {2} (length {3})".format(left_contig, left_ref_start - end , left_ref_start - start, end - start), file=sys.stdout)
                                        sv_evidences.append(SVEvidence(left_contig, left_ref_start - end, left_ref_start - start, typ, "kmer", read_id))
                                    if typ == "ins":
                                        print("Insertion detected: {0}: {1} - {2} (length {3})".format(left_contig, left_ref_start - start, left_ref_start - start + (end - start), end - start), file=sys.stdout)
                                        sv_evidences.append(SVEvidence(left_contig, left_ref_start - start, left_ref_start - start + (end - start), typ, "kmer", read_id))
                        elif percent_shift < parameters.tail_min_deviation:
                            deletion_size_estimate = reference_dist - individual_dist + (0.04 * original_length)
                            print("Deletion detected between {0} and {1} on {2} (estimated size: {3} bps, read {4})".format(right_ref_end, left_ref_start, left_contig, deletion_size_estimate, read_id), file=sys.stderr)
                            if deletion_size_estimate < 10000:
                                read_snippet = str(reads[read].seq[parameters.tail_span - left_q_start : original_length - right_q_end].upper())
                                ref_snippet = str(reference[left_contig].seq[right_ref_end:left_ref_start].upper().reverse_complement())
                                sv_results = find_svs(ref_snippet, read_snippet, parameters, debug = False)
                                for typ, start, end in sv_results:
                                    if typ == "del":
                                        print("Deletion detected: {0}: {1} - {2} (length {3})".format(left_contig, left_ref_start - end, left_ref_start - start, end - start), file=sys.stdout)
                                        sv_evidences.append(SVEvidence(left_contig, left_ref_start - end, left_ref_start - start, typ, "kmer", read_id))
                                    if typ == "ins":
                                        print("Insertion detected: {0}: {1} - {2} (length {3})".format(left_contig, left_ref_start - start, left_ref_start - start + (end - start), end - start), file=sys.stdout)
                                        sv_evidences.append(SVEvidence(left_contig, left_ref_start - start, left_ref_start - start + (end - start), typ, "kmer", read_id))

                        full_aln = full_aln_dict[read][0]
                        indels = find_indels_in_cigar_tuples(full_aln.cigartuples)
                        for pos, length, typ in indels:
                            sv_evidences.append(
                                SVEvidence(left_contig, full_aln.reference_start + pos, full_aln.reference_start + pos + length, typ, "cigar", read_id))

                    else:
                        pass # Todo
                elif not left_aln.is_reverse and not right_aln.is_reverse:
                    reference_dist = right_ref_start - left_ref_end
                    if reference_dist > 0:
                        individual_dist = original_length - left_q_end - (parameters.tail_span - right_q_start)
                        percent_shift = (individual_dist - reference_dist) / float(original_length)
                        #print("{0}\t{1}\t{2}\t{3}".format(left_contig, left_ref_end, right_ref_start, percent_shift), file=tail_diff_output);

                        if percent_shift > parameters.tail_max_deviation:
                            insertion_size_estimate = individual_dist - reference_dist - (0.04 * original_length)
                            print("Insertion detected between {0} and {1} on {2} (estimated size: {3} bps, read {4})".format(left_ref_end, right_ref_start, left_contig, insertion_size_estimate, read_id), file=sys.stderr)
                            if insertion_size_estimate < 10000:
                                read_snippet = str(reads[read].seq[left_q_end:original_length - parameters.tail_span + right_q_start].upper())
                                ref_snippet = str(reference[left_contig].seq[left_ref_end:right_ref_start].upper())
                                sv_results = find_svs(ref_snippet, read_snippet, parameters, debug = False)
                                for typ, start, end in sv_results:
                                    if typ == "del":
                                        #print("{0}\t{1}\t{2}\t{3}\tprecise".format(left_contig, left_ref_end + (start * 50), left_ref_end + (end * 50), (end - start) * 50), file=deletion_output)
                                        print("Deletion detected: {0}: {1} - {2} (length {3})".format(left_contig, left_ref_end + start, left_ref_end + end, end - start), file=sys.stdout)
                                        sv_evidences.append(SVEvidence(left_contig, left_ref_end + start, left_ref_end + end, typ, "kmer", read_id))
                                    if typ == "ins":
                                        #print("{0}\t{1}\t{2}\t{3}\tprecise".format(left_contig, left_ref_end + (start * 50), left_ref_end + (end * 50), (end - start) * 50), file=insertion_output)
                                        print("Insertion detected: {0}: {1} - {2} (length {3})".format(left_contig, left_ref_end + start, left_ref_end + end, end - start), file=sys.stdout)
                                        sv_evidences.append(SVEvidence(left_contig, left_ref_end + start, left_ref_end + end, typ, "kmer", read_id))
                        elif percent_shift < parameters.tail_min_deviation:
                            deletion_size_estimate = reference_dist - individual_dist + (0.04 * original_length)
                            print("Deletion detected between {0} and {1} on {2} (estimated size: {3} bps, read {4})".format(left_ref_end, right_ref_start, left_contig, deletion_size_estimate, read_id), file=sys.stderr)
                            if deletion_size_estimate < 10000:
                                read_snippet = str(reads[read].seq[left_q_end:original_length - parameters.tail_span + right_q_start].upper())
                                ref_snippet = str(reference[left_contig].seq[left_ref_end:right_ref_start].upper())
                                sv_results = find_svs(ref_snippet, read_snippet, parameters, debug = False)
                                for typ, start, end in sv_results:
                                    if typ == "del":
                                        #print("{0}\t{1}\t{2}\t{3}\tprecise".format(left_contig, left_ref_end + (start * 50), left_ref_end + (end * 50), (end - start) * 50), file=deletion_output)
                                        print("Deletion detected: {0}: {1} - {2} (length {3})".format(left_contig, left_ref_end + start, left_ref_end + end, end - start), file=sys.stdout)
                                        sv_evidences.append(SVEvidence(left_contig, left_ref_end + start, left_ref_end + end, typ, "kmer", read_id))
                                    if typ == "ins":
                                        #print("{0}\t{1}\t{2}\t{3}\tprecise".format(left_contig, left_ref_end + (start * 50), left_ref_end + (end * 50), (end - start) * 50), file=insertion_output)
                                        print("Insertion detected: {0}: {1} - {2} (length {3})".format(left_contig, left_ref_end + start, left_ref_end + end, end - start), file=sys.stdout)
                                        sv_evidences.append(SVEvidence(left_contig, left_ref_end + start, left_ref_end + end, typ, "kmer", read_id))

                        full_aln = full_aln_dict[read][0]
                        indels = find_indels_in_cigar_tuples(full_aln.cigartuples)
                        for pos, length, typ in indels:
                            sv_evidences.append(
                                SVEvidence(left_contig, full_aln.reference_start + pos, full_aln.reference_start + pos + length, typ, "cigar", read_id))

                    else:
                        pass # Todo
                elif left_aln.is_reverse and not right_aln.is_reverse:
                    print("Inversion detected on {0} between {1} and {2}".format(left_contig, left_ref_start, right_ref_end), file=sys.stdout)
                else:
                    print("Inversion detected on {0} between {1} and {2}".format(left_contig, left_ref_start, right_ref_end), file=sys.stdout)
    except KeyboardInterrupt:
        print("ERROR: Keyboard interrupt received. Interrupting SV detection..", file=sys.stderr)
    return sv_evidences


def mean_distance(evidence1, evidence2):
    """Return distance between means of two evidences."""
    if evidence1.contig == evidence2.contig and evidence1.type == evidence2.type:
        return abs(((evidence1.start + evidence1.end) / 2) - ((evidence2.start + evidence2.end) / 2))
    else:
        return float("inf")


def gowda_diday_distance(evidence1, evidence2, largest_indel_size):
    """Return Gowda-Diday distance between two evidences."""
    # different chromosomes
    if evidence1.contig != evidence2.contig:
        return float("inf")
    # different SV type
    if evidence1.type != evidence2.type:
        return float("inf")
    # non-intersecting
    if evidence1.end <= evidence2.start or evidence2.end <= evidence1.start:
        return float("inf")
    dist_pos = abs(evidence1.start - evidence2.start) / float(largest_indel_size)
    span1 = abs(evidence1.end - evidence1.start)
    span2 = abs(evidence2.end - evidence2.start)
    span_total = abs(max(evidence1.end, evidence2.end) - min(evidence1.start, evidence2.start))
    dist_span = abs(span1 - span2) / float(span_total)
    inter = min(evidence1.end, evidence2.end) - max(evidence1.start, evidence2.start)
    dist_content = (span1 + span2 - 2 * inter) / float(span_total)
    return dist_pos + dist_span + dist_content


def form_partitions(sv_evidences, max_delta=1000):
    """Form partitions of evidences using mean distance."""
    sorted_evidences = sorted(sv_evidences, key=lambda evi: (evi.type, evi.contig, (evi.start + evi.end) / 2))
    partitions = []
    current_partition = []
    for evidence in sorted_evidences:
        if len(current_partition) < 1:
            current_partition.append(evidence)
            continue
        if mean_distance(current_partition[0], evidence) > max_delta:
            partitions.append(current_partition[:])
            while len(current_partition) > 0 and mean_distance(current_partition[0], evidence) > max_delta:
                current_partition.pop(0)
        current_partition.append(evidence)
    partitions.append(current_partition[:])
    return partitions


def clusters_from_partitions(partitions, max_delta=1):
    """Form clusters in partitions using Gowda-Diday distance and clique finding in a distance graph."""
    clusters_full = []
    # Find clusters in each partition individually.
    for num, partition in enumerate(partitions):
        # print("Process partition", num, "of length", len(partition), file=sys.stderr)
        largest_evidence = sorted(partition, key=lambda evi: (evi.end - evi.start))[-1]
        largest_indel_size = largest_evidence.end - largest_evidence.start
        connection_graph = nx.Graph()
        connection_graph.add_nodes_from(range(len(partition)))
        for i1 in range(len(partition)):
            for i2 in range(len(partition)):
                if i1 == i2 or gowda_diday_distance(partition[i1], partition[i2], largest_indel_size) > max_delta:
                    pass
                else:
                    # Add edge in graph only if two indels are close to each other (distance <= max_delta)
                    connection_graph.add_edge(i1, i2)
        clusters_indices = nx.find_cliques(connection_graph)
        for cluster in clusters_indices:
            clusters_full.append([partition[index] for index in cluster])
    return clusters_full


def consolidate_clusters(clusters):
    """Consolidate clusters to a list of (type, contig, mean start, mean end, cluster size, members) tuples."""
    consolidated_clusters = []
    for cluster in clusters:
        length = len(cluster)
        cigar_evidences = [member for member in cluster if member.evidence == "cigar"]
        kmer_evidences = [member for member in cluster if member.evidence == "kmer"]
        score = len(cigar_evidences) + len(kmer_evidences)
        if len(cigar_evidences) > 0:
            score += 5
        if len(kmer_evidences) > 0:
            score += 5
        average_start = (2 * sum([member.start for member in cigar_evidences]) + sum([member.start for member in kmer_evidences])) / float(2*len(cigar_evidences) + len(kmer_evidences))
        average_end = (2 * sum([member.end for member in cigar_evidences]) + sum([member.end for member in kmer_evidences])) / float(2*len(cigar_evidences) + len(kmer_evidences))
        consolidated_clusters.append((cluster[0].type, cluster[0].contig,
                                      int(round(average_start)), int(round(average_end)),
                                      score, length, map(lambda x: x.as_tuple(), cluster)))
    return consolidated_clusters


def main():
    options = parse_arguments()
    parameters = callPacParams()
    parameters.set_with_options(options)

    # Run SV search only if raw SV results do not exist
    if not os.path.exists(options.temp_dir + '/sv_evidences.tsv'):
        create_temp_files(options.temp_dir, options.fasta, parameters.tail_span)
        run_alignments(options.temp_dir, options.genome, options.fasta)
        sv_evidences = analyze_read_tails(options.temp_dir, options.genome, options.fasta, parameters)

        raw_file = open(options.temp_dir + '/sv_evidences.tsv', 'w')
        for sv_evidence in sv_evidences:
            print("\t".join(map(str, sv_evidence.as_tuple())), file=raw_file)
        raw_file.close()
    else:
        print("WARNING: Result file with SV evidences already exists. Skip", file=sys.stderr)
        raw_file = open(options.temp_dir + '/sv_evidences.tsv', 'r')
        sv_evidences = []
        for line in raw_file:
            con, sta, end, typ, evi, sou = line.strip().split("\t")
            sv_evidences.append(SVEvidence(con, int(sta), int(end), typ, evi, sou))
        raw_file.close()

    # Cluster raw SVs
    partitions_raw = form_partitions(sv_evidences)
    print("Formed {0} partitions".format(len(partitions_raw)), file=sys.stderr)
    clusters_raw = clusters_from_partitions(partitions_raw)
    print("Subdivided partition into {0} clusters".format(len(clusters_raw)), file=sys.stderr)
    clusters_consolidated = consolidate_clusters(clusters_raw)
    partitions_consolidated = consolidate_clusters(partitions_raw)

    # Print output
    partition_output = open(options.temp_dir + '/partitions.bed', 'w')
    insertion_output = open(options.temp_dir + '/ins.bed', 'w')
    deletion_output = open(options.temp_dir + '/del.bed', 'w')
    inversion_output = open(options.temp_dir + '/inv.bed', 'w')

    for typ, contig, start, end, score, length, members in partitions_consolidated:
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(contig, start, end, score, length, members), file=partition_output)

    for typ, contig, start, end, score, length, members in clusters_consolidated:
        if typ == 'ins':
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(contig, start, end, score, length, members), file=insertion_output)
        if typ == 'del':
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(contig, start, end, score, length, members), file=deletion_output)
        if typ == 'inv':
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(contig, start, end, score, length, members), file=inversion_output)


if __name__ == "__main__":
    sys.exit(main())
