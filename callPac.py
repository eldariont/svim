import sys
import argparse
import os
import pysam
from subprocess import call, Popen
import subprocess
from collections import defaultdict

def parseArguments(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""""")
    parser.add_argument('fasta', type=argparse.FileType('r'))
    parser.add_argument('temp_dir', type=str, help='temp directory')
    return parser.parse_args()


def parseAlignmentFile(file_path, mode = 'r'):
    """Parses SAM and BAM files and returns two dictionaries: dict of reads (list of alignments for each read) and dict for mapping between chrID and ChrName"""
    sam = pysam.AlignmentFile(file_path, mode)

    alns = sam.fetch()
    aln_dict = defaultdict(list)
    chrIdToName = {}
    for aln in alns:
        aln_dict[aln.query_name].append(aln)
        if aln.reference_id >= 0 and not aln.reference_id in chrIdToName.keys():
            chrIdToName[aln.reference_id] = sam.get_reference_name(aln.reference_id)
    return aln_dict, chrIdToName


def parseCigarTuples(tuples, length, type):
    """Parses CIGAR tuples (op, len) and returns events of type 'type' that have a length similar to 'length'"""
    pos = 0
    candidates = []
    for op, l in tuples:
        if op == 0: #alignment match
            pos += l
        elif op == 1: #insertion
            if type == 'ins':
                if abs(length - l) < 300:
                    candidates.append((pos, l))
        elif op == 2: #deletion
            if type == 'del':
                if abs(length - l) < 300:
                    candidates.append((pos, l))
            pos += l
        elif op == 7 or op == 8: #match or mismatch
            pos += l
    return candidates


options = parseArguments(sys.argv)

if not os.path.exists(options.temp_dir):
    print("ERROR: Given temp directory does not exist", file=sys.stderr)
    sys.exit()

if not os.path.exists(options.temp_dir + '/full.fa'):
    full_file = open(options.temp_dir + '/full.fa', 'w')
    write_full = True
else:
    write_full = False
    print("WARNING: Temp file for full sequences exists. Skip", file=sys.stderr)

if not os.path.exists(options.temp_dir + '/left.fa'):
    left_file = open(options.temp_dir + '/left.fa', 'w')
    write_left = True
else:
    write_left = False
    print("WARNING: Temp file for left sequences exists. Skip", file=sys.stderr)

if not os.path.exists(options.temp_dir + '/right.fa'):
    right_file = open(options.temp_dir + '/right.fa', 'w')
    write_right = True
else:
    write_right = False
    print("WARNING: Temp file for right sequences exists. Skip", file=sys.stderr)

span = 1000

length_dict = {}
seq_nr = 0
for line in options.fasta:
    if line.startswith('>'):
        seq_nr += 1
        read_name = line.strip()[1:]
    else:
        sequence = line.strip()
        length = len(sequence)
        length_dict[read_name] = length
        if length < 2000:
            #print("WARNING: Sequence too short (Length {0})".format(length), file=sys.stderr)
            continue

        if write_full:
            print(">" + read_name, file=full_file)
            print(sequence, file=full_file)

        if write_left:
            prefix = sequence[:span]
            print(">" + read_name, file=left_file)
            print(prefix, file=left_file)

        if write_right:
            suffix = sequence[-span:]
            print(">" + read_name, file=right_file)
            print(suffix, file=right_file)


options.fasta.close()
if write_full:
    full_file.close()
if write_left:
    left_file.close()
if write_right:
    right_file.close()

print("INFO: Temporary files written", file=sys.stderr)

if not os.path.exists(options.temp_dir + '/full_aln.chained.sorted.bam'):
    ngmlr = Popen(['/home/heller_d/bin/miniconda2/bin/ngmlr', '-t', '30', '-r', '/scratch/cluster/heller_d/genomes/hg19/hg19.fa', '-q', options.temp_dir + '/full.fa', ], stdout=subprocess.PIPE)
    view = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools', 'view', '-b', '-@', '10'], stdin=ngmlr.stdout, stdout=subprocess.PIPE)
    sort = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools', 'sort', '-n', '-@', '10', '-o', options.temp_dir + '/full_aln.querysorted.bam'], stdin=view.stdout)
    sort.wait()
    if call(['python', '/home/heller_d/bin/bamChain', options.temp_dir + '/full_aln.querysorted.bam', options.temp_dir + '/full_aln.chained.bam', '--minmapq', '40']) != 0:
        print("ERROR: Calling bamchain on full sequences failed", file=sys.stderr)
    if call(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools', 'sort', '-@', '10', options.temp_dir + '/full_aln.chained.bam', '-o', options.temp_dir + '/full_aln.chained.sorted.bam']) != 0:
        print("ERROR: Calling samtools sort on full sequences failed", file=sys.stderr)
    if call(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools', 'index', options.temp_dir + '/full_aln.chained.sorted.bam']) != 0:
        print("ERROR: Calling samtools index on full sequences failed", file=sys.stderr)

else:
    print("WARNING: Alignment for full sequences exists. Skip", file=sys.stderr)

if not os.path.exists(options.temp_dir + '/left_aln.sam'):
    left_aln_file = open(options.temp_dir + '/left_aln.sam', 'w')
    if call(['/scratch/ngsvin/bin/bwa.kit/bwa', 'mem', '-x', 'pacbio', '-t', '30', '/scratch/cluster/heller_d/genomes/hg19/hg19.fa', options.temp_dir + '/left.fa'], stdout=left_aln_file) != 0:
        print("ERROR: Calling bwa on left sequences failed", file=sys.stderr)
    left_aln_file.close()
else:
    print("WARNING: Alignment for left sequences exists. Skip", file=sys.stderr)

if not os.path.exists(options.temp_dir + '/right_aln.sam'):
    right_aln_file = open(options.temp_dir + '/right_aln.sam', 'w')
    if call(['/scratch/ngsvin/bin/bwa.kit/bwa', 'mem', '-x', 'pacbio', '-t', '30', '/scratch/cluster/heller_d/genomes/hg19/hg19.fa', options.temp_dir + '/right.fa'], stdout=right_aln_file) != 0:
        print("ERROR: Calling bwa on right sequences failed", file=sys.stderr)
    right_aln_file.close()
else:
    print("WARNING: Alignment for right sequences exists. Skip", file=sys.stderr)

print("INFO: Alignment finished", file=sys.stderr)

full_aln_dict, full_name_table = parseAlignmentFile(options.temp_dir + '/full_aln.chained.sorted.bam', 'rb')
left_aln_dict, left_name_table = parseAlignmentFile(options.temp_dir + '/left_aln.sam')
right_aln_dict, right_name_table = parseAlignmentFile(options.temp_dir + '/right_aln.sam')

print("INFO: Alignment parsing finished", file=sys.stderr)

insertion_output = open(options.temp_dir + '/ins.bed', 'w')
deletion_output = open(options.temp_dir + '/del.bed', 'w')
inversion_output = open(options.temp_dir + '/inv.bed', 'w')
#tail_diff_output = open(options.temp_dir + '/tail_diff.tsv', 'w')

for read in left_aln_dict.keys():
    read_id = read.split("_")[1]
    original_length = length_dict[read]

    left_aln = left_aln_dict[read][0]
    try:
        right_aln = right_aln_dict[read][0]
    except IndexError:
        print("Right tail is unmapped", file=sys.stderr)
        continue

    if left_aln.is_unmapped or right_aln.is_unmapped:
        print("One or both of the tails is unmapped", file=sys.stderr)
        continue

    left_ref_start = left_aln.reference_start
    left_ref_end = left_aln.reference_end
    left_q_start = left_aln.query_alignment_start
    left_q_end = left_aln.query_alignment_end
    left_chrId = left_name_table[left_aln.reference_id]

    right_ref_start = right_aln.reference_start
    right_ref_end = right_aln.reference_end
    right_q_start = right_aln.query_alignment_start
    right_q_end = right_aln.query_alignment_end
    right_chrId = right_name_table[right_aln.reference_id]

    try:
        full_aln = full_aln_dict[read][0]
        full_ref_start = full_aln.reference_start
        full_ref_end = full_aln.reference_end
        full_q_start = full_aln.query_alignment_start
        full_q_end = full_aln.query_alignment_end
    except IndexError:
        full_aln = None

    if left_chrId != right_chrId:
        print("Tails map to different chromosomes ({0}, {1})".format(left_chrId, right_chrId), file=sys.stderr)
        continue

    if left_aln.is_reverse and right_aln.is_reverse:
        reference_dist = left_ref_start - right_ref_end
        if reference_dist > 0:
            individual_dist = original_length - right_q_end - (span - left_q_start)
            percent_shift = (individual_dist - reference_dist) / float(original_length)
            #print("{0}\t{1}\t{2}\t{3}".format(left_chrId, right_ref_end, left_ref_start, percent_shift), file=tail_diff_output);

            if percent_shift > 0.1:
                insertion_size_estimate = individual_dist - reference_dist - (0.04 * original_length)
                #print("Insertion detected between {0} and {1} on {2} (estimated size: {3} bps)".format(right_ref_end, left_ref_start, left_chrId, insertion_size_estimate), file=sys.stderr)
                if insertion_size_estimate < 10000:
                    if not full_aln is None:
                        full_aln_left_shift = abs((right_q_start - full_q_start) - (right_ref_start - full_ref_start))
                        full_aln_right_shift = abs((left_ref_end + (span - left_q_end)) - (full_ref_end + (original_length - full_q_end)))
                        #print("Read {0} ins1: left shift: {1} right shift: {2}".format(read_id, full_aln_left_shift, full_aln_right_shift), file=sys.stderr)
                        if full_aln_left_shift < 10 and full_aln_right_shift < 10:
                            #print("CIGAR string: {0}".format(full_aln.cigarstring), file=sys.stderr)
                            candidates = parseCigarTuples(full_aln.cigartuples, insertion_size_estimate, 'ins')
                            if len(candidates) == 1:
                                pos, length = candidates[0]
                                print("{0}\t{1}\t{2}\t{3}\tprecise".format(left_chrId, full_ref_start + pos, full_ref_start + pos + length, length), file=insertion_output)
                                print("Insertion detected: {0}: {1} - {2} (length {3})".format(left_chrId, full_ref_start + pos, full_ref_start + pos + length, length), file=sys.stdout)
                                continue
                    #print("{0}\t{1}\t{2}\t{3}\timprecise".format(left_chrId, right_ref_end, left_ref_start, insertion_size_estimate), file=insertion_output)

            if percent_shift < 0:
                deletion_size_estimate = reference_dist - individual_dist + (0.04 * original_length)
                #print("Deletion detected between {0} and {1} on {2} (estimated size: {3} bps)".format(right_ref_end, left_ref_start, left_chrId, deletion_size_estimate), file=sys.stderr)
                if deletion_size_estimate < 10000:
                    if not full_aln is None:
                        full_aln_left_shift = abs((right_q_start - full_q_start) - (right_ref_start - full_ref_start))
                        full_aln_right_shift = abs((left_ref_end + (span - left_q_end)) - (full_ref_end + (original_length - full_q_end)))
                        #print("Read {0} del1: left shift: {1} right shift: {2}".format(read_id, full_aln_left_shift, full_aln_right_shift), file=sys.stderr)
                        if full_aln_left_shift < 10 and full_aln_right_shift < 10:
                            #print("CIGAR string: {0}".format(full_aln.cigarstring), file=sys.stderr)
                            candidates = parseCigarTuples(full_aln.cigartuples, deletion_size_estimate, 'del')
                            if len(candidates) == 1:
                                pos, length = candidates[0]
                                print("{0}\t{1}\t{2}\t{3}\tprecise".format(left_chrId, full_ref_start + pos, full_ref_start + pos + length, length), file=deletion_output)
                                print("Deletion detected: {0}: {1} - {2} (length {3})".format(left_chrId, full_ref_start + pos, full_ref_start + pos + length, length), file=sys.stdout)
                                continue
                    #print("{0}\t{1}\t{2}\t{3}\timprecise".format(left_chrId, right_ref_end, left_ref_start, deletion_size_estimate), file=deletion_output)
        else:
            pass #Todo
    elif not left_aln.is_reverse and not right_aln.is_reverse:
        reference_dist = right_ref_start - left_ref_end
        if reference_dist > 0:
            individual_dist = original_length - left_q_end - (span - right_q_start)
            percent_shift = (individual_dist - reference_dist) / float(original_length)
            #print("{0}\t{1}\t{2}\t{3}".format(left_chrId, left_ref_end, right_ref_start, percent_shift), file=tail_diff_output);

            if percent_shift > 0.2:
                insertion_size_estimate = individual_dist - reference_dist - (0.04 * original_length)
                #print("Insertion detected between {0} and {1} on {2} (estimated size: {3} bps)".format(left_ref_end, right_ref_start, left_chrId, insertion_size_estimate), file=sys.stderr)
                if insertion_size_estimate < 10000:
                    if not full_aln is None:
                        full_aln_left_shift = abs((left_q_start - full_q_start) - (left_ref_start - full_ref_start))
                        full_aln_right_shift = abs((right_ref_end + (span - right_q_end)) - (full_ref_end + (original_length - full_q_end)))
                        #print("Read {0} ins2: left shift: {1} right shift: {2}".format(read_id, full_aln_left_shift, full_aln_right_shift), file=sys.stderr)
                        if full_aln_left_shift < 10 and full_aln_right_shift < 10:
                            #print("CIGAR string: {0}".format(full_aln.cigarstring), file=sys.stderr)
                            candidates = parseCigarTuples(full_aln.cigartuples, insertion_size_estimate, 'ins')
                            if len(candidates) == 1:
                                pos, length = candidates[0]
                                print("{0}\t{1}\t{2}\t{3}\tprecise".format(left_chrId, full_ref_start + pos, full_ref_start + pos + length, length), file=insertion_output)
                                print("Insertion detected: {0}: {1} - {2} (length {3})".format(left_chrId, full_ref_start + pos, full_ref_start + pos + length, length), file=sys.stdout)
                                continue
                    #print("{0}\t{1}\t{2}\t{3}\timprecise".format(left_chrId, left_ref_end, right_ref_start, insertion_size_estimate), file=insertion_output)


            if percent_shift < -0.1:
                deletion_size_estimate = reference_dist - individual_dist + (0.04 * original_length)
                #print("Deletion detected between {0} and {1} on {2} (estimated size: {3} bps)".format(left_ref_end, right_ref_start, left_chrId, deletion_size_estimate), file=sys.stderr)
                if deletion_size_estimate < 10000:
                    if not full_aln is None:
                        full_aln_left_shift = abs((left_q_start - full_q_start) - (left_ref_start - full_ref_start))
                        full_aln_right_shift = abs((right_ref_end + (span - right_q_end)) - (full_ref_end + (original_length - full_q_end)))
                        #print("Read {0} del2: left shift: {1} right shift: {2}".format(read_id, full_aln_left_shift, full_aln_right_shift), file=sys.stderr)
                        if full_aln_left_shift < 10 and full_aln_right_shift < 10:
                            #print("CIGAR string: {0}".format(full_aln.cigarstring), file=sys.stderr)
                            candidates = parseCigarTuples(full_aln.cigartuples, deletion_size_estimate, 'del')
                            if len(candidates) == 1:
                                pos, length = candidates[0]
                                print("{0}\t{1}\t{2}\t{3}\tprecise".format(left_chrId, full_ref_start + pos, full_ref_start + pos + length, length), file=deletion_output)
                                print("Deletion detected: {0}: {1} - {2} (length {3})".format(left_chrId, full_ref_start + pos, full_ref_start + pos + length, length), file=sys.stdout)
                                continue
                    #print("{0}\t{1}\t{2}\t{3}\timprecise".format(left_chrId, left_ref_end, right_ref_start, deletion_size_estimate), file=deletion_output)
        else:
            pass #Todo
    elif left_aln.is_reverse and not right_aln.is_reverse:
        print("Inversion detected between {0} and {1} on {2}".format(left_chrId, left_ref_start, right_ref_end), file=sys.stdout)
        print("{0}\t{1}\t{2}".format(left_chrId, left_ref_start, right_ref_end), file=inversion_output)
    else:
        print("Inversion detected between {0} and {1} on {2}".format(left_chrId, left_ref_start, right_ref_end), file=sys.stdout)
        print("{0}\t{1}\t{2}".format(left_chrId, left_ref_start, right_ref_end), file=inversion_output)
