from __future__ import print_function

import sys
import argparse
import os
from subprocess import call, Popen, PIPE

import pysam
from Bio import SeqIO

from callPacParams import callPacParams
from kmerCounting import find_svs
from callPacIndels import parse_sam_file, form_partitions, clusters_from_partitions, consolidate_clusters


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""""")
    parser.add_argument('fasta', type=argparse.FileType('r'))
    parser.add_argument('genome', default="/scratch/cluster/heller_d/genomes/hg19/hg19.fa", type=str)
    parser.add_argument('temp_dir', type=str, help='temp directory')
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
        seq_nr = 0
        for line in fasta:
            if line.startswith('>'):
                seq_nr += 1
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


def run_alignments(temp_dir, genome):
    """Align left and right read tails with BWA MEM."""
    if not os.path.exists(temp_dir + '/left_aln.sorted.bam'):
        bwa = Popen(['/scratch/ngsvin/bin/bwa.kit/bwa',
                     'mem', '-x', 'pacbio', '-t', '30', genome, temp_dir + '/left.fa'], stdout=PIPE)
        view = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'view', '-b', '-@', '10'], stdin=bwa.stdout, stdout=PIPE)
        sort = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'sort', '-@', '10', '-o', temp_dir + '/left_aln.sorted.bam'], stdin=view.stdout)
        sort.wait()
        if call(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                 'index', temp_dir + '/left_aln.sorted.bam']) != 0:
            print("ERROR: Calling samtools index on left sequences failed", file=sys.stderr)
    else:
        print("WARNING: Alignment for left sequences exists. Skip", file=sys.stderr)

    if not os.path.exists(temp_dir + '/right_aln.sorted.bam'):
        bwa = Popen(['/scratch/ngsvin/bin/bwa.kit/bwa',
                     'mem', '-x', 'pacbio', '-t', '30', genome, temp_dir + '/right.fa'], stdout=PIPE)
        view = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'view', '-b', '-@', '10'], stdin=bwa.stdout, stdout=PIPE)
        sort = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'sort', '-@', '10', '-o', temp_dir + '/right_aln.sorted.bam'], stdin=view.stdout)
        sort.wait()
        if call(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                 'index', temp_dir + '/right_aln.sorted.bam']) != 0:
            print("ERROR: Calling samtools index on right sequences failed", file=sys.stderr)
    else:
        print("WARNING: Alignment for right sequences exists. Skip", file=sys.stderr)

    print("INFO: Alignment finished", file=sys.stderr)


def search_svs(temp_dir, genome, fasta, parameters):
    """Search for SVs using aligned read tails and the read and reference regions in-between."""
    left_sam = pysam.AlignmentFile(temp_dir + '/left_aln.sorted.bam')
    right_sam = pysam.AlignmentFile(temp_dir + '/right_aln.sorted.bam')
    left_contigs = left_sam.references

    reference = SeqIO.index(genome, "fasta")
    reads = SeqIO.index(fasta.name, "fasta")
    print("INFO: Opening output and reference finished", file=sys.stderr)

    found_svs = []

    try:
        for left_contig in left_contigs:
            print("INFO: Searching for SVs in", left_contig, "..", file=sys.stderr)
            left_aln_dict = parse_sam_file(left_sam, left_contig)
            right_aln_dict = parse_sam_file(right_sam, left_contig)

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
                                        found_svs.append( (left_contig, left_ref_start - end, left_ref_start - start, typ) )
                                    if typ == "ins":
                                        print("Insertion detected: {0}: {1} - {2} (length {3})".format(left_contig, left_ref_start - start, left_ref_start - start + (end - start), end - start), file=sys.stdout)
                                        found_svs.append( (left_contig, left_ref_start - start, left_ref_start - start + (end - start), typ) )

                        if percent_shift < parameters.tail_min_deviation:
                            deletion_size_estimate = reference_dist - individual_dist + (0.04 * original_length)
                            print("Deletion detected between {0} and {1} on {2} (estimated size: {3} bps, read {4})".format(right_ref_end, left_ref_start, left_contig, deletion_size_estimate, read_id), file=sys.stderr)
                            if deletion_size_estimate < 10000:
                                read_snippet = str(reads[read].seq[parameters.tail_span - left_q_start : original_length - right_q_end].upper())
                                ref_snippet = str(reference[left_contig].seq[right_ref_end:left_ref_start].upper().reverse_complement())
                                sv_results = find_svs(ref_snippet, read_snippet, parameters, debug = False)
                                for typ, start, end in sv_results:
                                    if typ == "del":
                                        print("Deletion detected: {0}: {1} - {2} (length {3})".format(left_contig, left_ref_start - end, left_ref_start - start, end - start), file=sys.stdout)
                                        found_svs.append( (left_contig, left_ref_start - end, left_ref_start - start, typ) )
                                    if typ == "ins":
                                        print("Insertion detected: {0}: {1} - {2} (length {3})".format(left_contig, left_ref_start - start, left_ref_start - start + (end - start), end - start), file=sys.stdout)
                                        found_svs.append( (left_contig, left_ref_start - start, left_ref_start - start + (end - start), typ) )
                    else:
                        pass #Todo
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
                                        found_svs.append( (left_contig, left_ref_end + start, left_ref_end + end, typ) )
                                    if typ == "ins":
                                        #print("{0}\t{1}\t{2}\t{3}\tprecise".format(left_contig, left_ref_end + (start * 50), left_ref_end + (end * 50), (end - start) * 50), file=insertion_output)
                                        print("Insertion detected: {0}: {1} - {2} (length {3})".format(left_contig, left_ref_end + start, left_ref_end + end, end - start), file=sys.stdout)
                                        found_svs.append( (left_contig, left_ref_end + start, left_ref_end + end, typ) )

                        if percent_shift < parameters.tail_min_deviation:
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
                                        found_svs.append( (left_contig, left_ref_end + start, left_ref_end + end, typ) )
                                    if typ == "ins":
                                        #print("{0}\t{1}\t{2}\t{3}\tprecise".format(left_contig, left_ref_end + (start * 50), left_ref_end + (end * 50), (end - start) * 50), file=insertion_output)
                                        print("Insertion detected: {0}: {1} - {2} (length {3})".format(left_contig, left_ref_end + start, left_ref_end + end, end - start), file=sys.stdout)
                                        found_svs.append( (left_contig, left_ref_end + start, left_ref_end + end, typ) )
                    else:
                        pass #Todo
                elif left_aln.is_reverse and not right_aln.is_reverse:
                    print("Inversion detected on {0} between {1} and {2}".format(left_contig, left_ref_start, right_ref_end), file=sys.stdout)
                else:
                    print("Inversion detected on {0} between {1} and {2}".format(left_contig, left_ref_start, right_ref_end), file=sys.stdout)
    except KeyboardInterrupt:
        print("ERROR: Keyboard interrupt received. Interrupting SV detection..", file=sys.stderr)
    return found_svs


def main():
    options = parse_arguments()
    parameters = callPacParams()

    # Run SV search only if raw SV results do not exist
    if not os.path.exists(options.temp_dir + '/svs_raw.tsv'):
        create_temp_files(options.temp_dir, options.fasta, parameters.tail_span)
        run_alignments(options.temp_dir, options.genome)
        svs_raw = search_svs(options.temp_dir, options.genome, options.fasta, parameters)
        
        raw_file = open(options.temp_dir + '/svs_raw.tsv', 'w')
        for sv in svs_raw:
            print("\t".join(map(str, sv)), file=raw_file)
        raw_file.close()
    else:
        print("WARNING: Result file with raw SVs already exists. Skip", file=sys.stderr)
        raw_file = open(options.temp_dir + '/svs_raw.tsv', 'r')
        svs_raw = []
        for line in raw_file:
            c, s, e, t = line.strip().split("\t")
            svs_raw.append((c, int(s), int(e), t))
        raw_file.close()

    # Cluster raw SVs
    partitions_raw = form_partitions(svs_raw)
    print("Formed {0} partitions".format(len(partitions_raw)), file=sys.stderr)
    clusters_raw = clusters_from_partitions(partitions_raw, sorted([e - s for c, s, e, t in svs_raw])[-1])
    print("Subdivided partition into {0} clusters".format(len(clusters_raw)), file=sys.stderr)
    clusters_consolidated = consolidate_clusters(clusters_raw)
    partitions_consolidated = consolidate_clusters(partitions_raw)

    # Print output
    partition_output = open(options.temp_dir + '/partitions.bed', 'w')
    insertion_output = open(options.temp_dir + '/ins.bed', 'w')
    deletion_output = open(options.temp_dir + '/del.bed', 'w')
    
    for typ, contig, start, end, support, members in partitions_consolidated:
        print("{0}\t{1}\t{2}\t{3}\t{4}".format(contig, start, end, support, members), file=partition_output)
        
    for typ, contig, start, end, support, members in clusters_consolidated:
        if typ == 'ins':
            print("{0}\t{1}\t{2}\t{3}\t{4}".format(contig, start, end, support, members), file=insertion_output)
        if typ == 'del':
            print("{0}\t{1}\t{2}\t{3}\t{4}".format(contig, start, end, support, members), file=deletion_output)


if __name__ == "__main__":
    sys.exit(main())
