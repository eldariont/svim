from __future__ import print_function

import sys
import argparse
import os
from subprocess import call, Popen, PIPE
from collections import defaultdict

import pysam
import networkx as nx


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""""")
    parser.add_argument('fasta', type=argparse.FileType('r'))
    parser.add_argument('genome', default="/scratch/cluster/heller_d/genomes/hg19/hg19.fa", type=str)
    parser.add_argument('temp_dir', type=str, help='temp directory')
    return parser.parse_args()


def parse_sam_file(sam, contig):
    """Parses a SAM file and returns a dict of reads (list of alignments for each read) for a given reference contig"""
    alns = sam.fetch(reference=contig)
    aln_dict = defaultdict(list)
    for aln in alns:
        aln_dict[aln.query_name].append(aln)
    return aln_dict


def find_indels_in_cigar_tuples(tuples, min_length=50):
    """Parses CIGAR tuples (op, len) and returns Indels with a length > minLength"""
    pos = 0
    indels = []
    for operation, length in tuples:
        if operation == 0:                     # alignment match
            pos += length
        elif operation == 1:                   # insertion
            if length >= min_length:
                indels.append((pos, length, 'ins'))
        elif operation == 2:                   # deletion
            if length >= min_length:
                indels.append((pos, length, 'del'))
            pos += length
        elif operation == 7 or operation == 8:        # match or mismatch
            pos += length
    return indels


def mean_distance(indel1, indel2):
    """Return distance between means of two indels.
       An indel is a tuple of (contig, start, end, type)."""
    if indel1[0] == indel2[0] and indel1[3] == indel2[3]:
        return abs(((indel1[1] + indel1[2]) / 2) - ((indel2[1] + indel2[2]) / 2))
    else:
        return float("inf")


def gowda_diday_distance(indel1, indel2, largest_indel_size):
    """Return Gowda-Diday distance between two indels.
       An indel is a tuple of (contig, start, end, type)."""
    # different chromosomes
    if indel1[0] != indel2[0]:
        return float("inf")
    # different SV type
    if indel1[3] != indel2[3]:
        return float("inf")
    # non-intersecting
    if indel1[2] <= indel2[1] or indel2[2] <= indel1[1]:
        return float("inf")
    dist_pos = abs(indel1[1] - indel2[1]) / float(largest_indel_size)
    span1 = abs(indel1[2] - indel1[1])
    span2 = abs(indel2[2] - indel2[1])
    span_total = abs(max(indel1[2], indel2[2]) - min(indel1[1], indel2[1]))
    dist_span = abs(span1 - span2) / float(span_total)
    inter = min(indel1[2], indel2[2]) - max(indel1[1], indel2[1])
    dist_content = (span1 + span2 - 2 * inter) / float(span_total)
    return dist_pos + dist_span + dist_content


def form_partitions(indels, max_delta=1000):
    """Form partitions of indels using mean distance.
       An indel is a tuple of (contig, start, end, type)."""
    # Sort indels by their type, contig and mean
    sorted_indels = sorted(indels, key=lambda ind: (ind[3], ind[0], (ind[1] + ind[2]) / 2))
    partitions = []
    current_partition = []
    for indel in sorted_indels:
        if len(current_partition) < 1:
            current_partition.append(indel)
            continue
        if mean_distance(current_partition[0], indel) > max_delta:
            partitions.append(current_partition[:])
            while len(current_partition) > 0 and mean_distance(current_partition[0], indel) > max_delta:
                current_partition.pop(0)
        current_partition.append(indel)
    partitions.append(current_partition[:])
    return partitions


def clusters_from_partitions(partitions, largest_indel_size, max_delta=1):
    """Form clusters in partitions using Gowda-Diday distance and clique finding in a distance graph."""
    clusters_full = []
    # Find clusters in each partition individually.
    for num, partition in enumerate(partitions):
        # print("Process partition", num, "of length", len(partition), file=sys.stderr)
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
    """Consolidate clusters to a list of (typ, contig, mean start, mean end, cluster size, members) tuples."""
    consolidated_clusters = []
    for cluster in clusters:
        length = len(cluster)
        starts = [member[1] for member in cluster]
        ends = [member[2] for member in cluster]
        consolidated_clusters.append((cluster[0][3], cluster[0][0],
                                      int(round(sum(starts) / length)), int(round(sum(ends)) / length),
                                      length, cluster))
    return consolidated_clusters


def main():
    options = parse_arguments()

    if not os.path.exists(options.temp_dir):
        print("ERROR: Given temp directory does not exist", file=sys.stderr)
        sys.exit()

    # Create temporary FASTA file
    if not os.path.exists(options.temp_dir + '/full.fa'):
        full_file = open(options.temp_dir + '/full.fa', 'w')
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

                print(">" + read_name, file=full_file)
                print(sequence, file=full_file)
        full_file.close()
    else:
        print("WARNING: Temp file for full sequences exists. Skip", file=sys.stderr)

    print("INFO: Temporary files written", file=sys.stderr)

    # Align full reads with NGM-LR
    if not os.path.exists(options.temp_dir + '/full_aln.chained.sorted.bam'):
        ngmlr = Popen(['/home/heller_d/bin/miniconda2/bin/ngmlr',
                       '-t', '30', '-r', options.genome, '-q', options.temp_dir + '/full.fa', ], stdout=PIPE)
        view = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'view', '-b', '-@', '10'], stdin=ngmlr.stdout, stdout=PIPE)
        sort = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'sort', '-n', '-@', '10', '-o', options.temp_dir + '/full_aln.querysorted.bam'], stdin=view.stdout)
        sort.wait()
        if call(['python', '/home/heller_d/bin/bamChain', options.temp_dir + '/full_aln.querysorted.bam',
                 options.temp_dir + '/full_aln.chained.bam', '--minmapq', '40']) != 0:
            print("ERROR: Calling bamchain on full sequences failed", file=sys.stderr)
        if call(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools', 'sort', '-@', '10',
                 options.temp_dir + '/full_aln.chained.bam', '-o', options.temp_dir + '/full_aln.chained.sorted.bam']) != 0:
            print("ERROR: Calling samtools sort on full sequences failed", file=sys.stderr)
        if call(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                 'index', options.temp_dir + '/full_aln.chained.sorted.bam']) != 0:
            print("ERROR: Calling samtools index on full sequences failed", file=sys.stderr)
    else:
        print("WARNING: Alignment for full sequences exists. Skip", file=sys.stderr)

    print("INFO: Alignment finished", file=sys.stderr)

    sam = pysam.AlignmentFile(options.temp_dir + '/full_aln.chained.sorted.bam', 'rb')
    contigs = sam.references

    indel_output = open(options.temp_dir + '/indels.bed', 'w')
    partition_output = open(options.temp_dir + '/partitions.bed', 'w')
    insertion_output = open(options.temp_dir + '/ins.bed', 'w')
    deletion_output = open(options.temp_dir + '/del.bed', 'w')

    large_indels = []
    largest_indel_size = -1

    for contig in contigs:
        print("INFO: Searching for SVs in", contig, "..", file=sys.stderr)
        full_aln_dict = parse_sam_file(sam, contig)
        for read in full_aln_dict.keys():
            full_aln = full_aln_dict[read][0]
            indels = find_indels_in_cigar_tuples(full_aln.cigartuples)
            for pos, length, typ in indels:
                large_indels.append(
                    (contig, full_aln.reference_start + pos, full_aln.reference_start + pos + length, typ))
                if length > largest_indel_size:
                    largest_indel_size = length

    for contig, start, end, typ in large_indels:
        print("{0}\t{1}\t{2}\t{3}".format(contig, start, end, typ), file=indel_output)

    partitions = form_partitions(large_indels)
    print("Formed {0} partitions".format(len(partitions)), file=sys.stderr)
    raw_clusters = clusters_from_partitions(partitions, largest_indel_size)
    print("Subdivided partition into {0} clusters".format(len(raw_clusters)), file=sys.stderr)
    clusters = consolidate_clusters(raw_clusters)
    consolidated_partitions = consolidate_clusters(partitions)

    for typ, contig, start, end, support, members in consolidated_partitions:
        print("{0}\t{1}\t{2}\t{3}\t{4}".format(contig, start, end, support, members), file=partition_output)

    for typ, contig, start, end, support, members in clusters:
        if typ == 'ins':
            print("{0}\t{1}\t{2}\t{3}\t{4}".format(contig, start, end, support, members), file=insertion_output)
        if typ == 'del':
            print("{0}\t{1}\t{2}\t{3}\t{4}".format(contig, start, end, support, members), file=deletion_output)


if __name__ == "__main__":
    sys.exit(main())
