from __future__ import print_function

import sys
import argparse
import os
from subprocess import call, Popen, PIPE

import pysam
from Bio import SeqIO
import networkx as nx

from callPacParams import callPacParams
from kmerCounting import find_svs
from callPacIndels import parse_sam_file, find_indels_in_cigar_tuples

class SVEvidence:
    def __init__(self, contig, start, end, type, evidence, source):
        self.contig = contig
        self.start = start
        self.end = end
        self.type = type
        self.evidence = evidence
        self.source = source

    def as_tuple(self):
        return (self.contig, self.start, self.end, self.type, self.evidence, self.source)

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

    # Align full reads with NGM-LR
    if not os.path.exists(temp_dir + '/full_aln.chained.sorted.bam'):
        ngmlr = Popen(['/home/heller_d/bin/miniconda2/bin/ngmlr',
                       '-t', '30', '-r', genome, '-q', os.path.realpath(fasta.name), ], stdout=PIPE)
        view = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'view', '-b', '-@', '10'], stdin=ngmlr.stdout, stdout=PIPE)
        sort = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'sort', '-n', '-@', '10', '-o', temp_dir + '/full_aln.querysorted.bam'],
                     stdin=view.stdout)
        sort.wait()
        if call(['python', '/home/heller_d/bin/bamChain', temp_dir + '/full_aln.querysorted.bam',
                 temp_dir + '/full_aln.chained.bam', '--minmapq', '40']) != 0:
            print("ERROR: Calling bamchain on full sequences failed", file=sys.stderr)
        if call(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools', 'sort', '-@', '10',
                 temp_dir + '/full_aln.chained.bam', '-o',
                 temp_dir + '/full_aln.chained.sorted.bam']) != 0:
            print("ERROR: Calling samtools sort on full sequences failed", file=sys.stderr)
        if call(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                 'index', temp_dir + '/full_aln.chained.sorted.bam']) != 0:
            print("ERROR: Calling samtools index on full sequences failed", file=sys.stderr)
    else:
        print("WARNING: Alignment for full sequences exists. Skip", file=sys.stderr)

    print("INFO: Alignment finished", file=sys.stderr)


def search_svs(temp_dir, genome, fasta, parameters):
    """Search for SVs using aligned read tails and the read and reference regions in-between."""
    left_sam = pysam.AlignmentFile(temp_dir + '/left_aln.sorted.bam')
    right_sam = pysam.AlignmentFile(temp_dir + '/right_aln.sorted.bam')
    full_sam = pysam.AlignmentFile(temp_dir + '/full_aln.chained.sorted.bam', 'rb')
    left_contigs = left_sam.references

    reference = SeqIO.index(genome, "fasta")
    reads = SeqIO.index(fasta.name, "fasta")
    print("INFO: Opening output and reference finished", file=sys.stderr)

    sv_evidences = []

    try:
        for left_contig in left_contigs:
            print("INFO: Searching for SVs in", left_contig, "..", file=sys.stderr)
            left_aln_dict = parse_sam_file(left_sam, left_contig)
            right_aln_dict = parse_sam_file(right_sam, left_contig)
            full_aln_dict = parse_sam_file(full_sam, left_contig)

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
        starts = [member.start for member in cluster]
        ends = [member.end for member in cluster]
        consolidated_clusters.append((cluster[0].type, cluster[0].contig,
                                      int(round(sum(starts) / length)), int(round(sum(ends)) / length),
                                      length, map(lambda x: x.as_tuple(), cluster)))
    return consolidated_clusters


def main():
    options = parse_arguments()
    parameters = callPacParams()
    parameters.set_with_options(options)

    # Run SV search only if raw SV results do not exist
    if not os.path.exists(options.temp_dir + '/sv_evidences.tsv'):
        create_temp_files(options.temp_dir, options.fasta, parameters.tail_span)
        run_alignments(options.temp_dir, options.genome, options.fasta)
        sv_evidences = search_svs(options.temp_dir, options.genome, options.fasta, parameters)
        
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
    
    for typ, contig, start, end, support, members in partitions_consolidated:
        print("{0}\t{1}\t{2}\t{3}\t{4}".format(contig, start, end, support, members), file=partition_output)
        
    for typ, contig, start, end, support, members in clusters_consolidated:
        if typ == 'ins':
            print("{0}\t{1}\t{2}\t{3}\t{4}".format(contig, start, end, support, members), file=insertion_output)
        if typ == 'del':
            print("{0}\t{1}\t{2}\t{3}\t{4}".format(contig, start, end, support, members), file=deletion_output)


if __name__ == "__main__":
    sys.exit(main())
