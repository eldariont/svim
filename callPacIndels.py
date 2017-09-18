import sys
import argparse
import os
import pysam
from subprocess import call, Popen
import subprocess
from collections import defaultdict
import networkx as nx

def parseArguments(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""""")
    parser.add_argument('fasta', type=argparse.FileType('r'))
    parser.add_argument('genome', default = "/scratch/cluster/heller_d/genomes/hg19/hg19.fa", type=str)
    parser.add_argument('temp_dir', type=str, help='temp directory')
    return parser.parse_args()


def parseAlignmentFile(sam, contig):
    """Parses SAM and BAM files and returns two dictionaries: dict of reads (list of alignments for each read) and dict for mapping between chrID and ChrName"""
    alns = sam.fetch(reference = contig)
    aln_dict = defaultdict(list)
    for aln in alns:
        aln_dict[aln.query_name].append(aln)
        #print(aln.query_name)
    return aln_dict


def findIndelsInCigarTuples(tuples, minLength = 50):
    """Parses CIGAR tuples (op, len) and returns Indels with a length > minLength"""
    pos = 0
    indels = []
    for op, l in tuples:
        if op == 0: #alignment match
            pos += l
        elif op == 1: #insertion
            if l >= minLength:
                indels.append((pos, l, 'ins'))
        elif op == 2: #deletion
            if l >= minLength:
                indels.append((pos, l, 'del'))
            pos += l
        elif op == 7 or op == 8: #match or mismatch
            pos += l
    return indels


def meanDistance(indel1, indel2):
    if indel1[0] == indel2[0] and indel1[3] == indel2[3]:
        return abs( ( ( indel1[1] + indel1[2] ) / 2 ) - ( ( indel2[1] + indel2[2] ) / 2 ) )
    else:
        return float("inf")


def hausdorffDistance(indel1, indel2):
    if indel1[0] == indel2[0] and indel1[3] == indel2[3]:
        return max( abs( indel1[1] - indel2[1] ), abs( indel1[2] - indel2[2] ) )
    else:
        return float("inf")


def gowdaDidayDistance(indel1, indel2, largestIndelSize):
    #different chromosomes
    if indel1[0] != indel2[0]:
        return float("inf")
    #different SV type
    if indel1[3] != indel2[3]:
        return float("inf")
    #non-intersecting
    if indel1[2] <= indel2[1] or indel2[2] <= indel1[1]:
        return float("inf")
    distPos = abs(indel1[1] - indel2[1]) / float(largestIndelSize)
    span1 = abs(indel1[2] - indel1[1])
    span2 = abs(indel2[2] - indel2[1])
    spanTotal = abs(max(indel1[2], indel2[2]) - min(indel1[1], indel2[1]))
    distSpan = abs(span1 - span2) / float(spanTotal)
    inter = min(indel1[2], indel2[2]) - max(indel1[1], indel2[1])
    distContent = (span1 + span2 - 2 * inter) / float(spanTotal)
    return distPos + distSpan + distContent


def formPartitions(largeIndels, maxDelta = 5000):
    #sort indels by their chromosome and mean
    sortedIndels = sorted(largeIndels, key=lambda indel: (indel[3], indel[0], ( indel[1] + indel[2] ) / 2 ))
    partitions = []
    currentPartition = []
    for indel in sortedIndels:
        if len(currentPartition) < 1:
            currentPartition.append(indel)
            continue
        if meanDistance(currentPartition[0], indel) > maxDelta:
            partitions.append(currentPartition[:])
            while len(currentPartition) > 0 and meanDistance(currentPartition[0], indel) > maxDelta:
                currentPartition.pop(0)
        currentPartition.append(indel)
    partitions.append(currentPartition[:])
    return partitions


def clustersFromPartitions(partitions, largestIndelSize, maxDelta = 1):
    clusters = set()
    for partition in partitions:
        for i1 in range(len(partition)):
            currentCluster = []
            for i2 in range(len(partition)):
                if i1 == i2 or gowdaDidayDistance(partition[i1], partition[i2], largestIndelSize) <= maxDelta:
                    currentCluster.append(partition[i2])
            clusters.add(tuple(elem for elem in currentCluster))
    return list(clusters)


def clustersFromPartitions2(partitions, largestIndelSize, maxDelta = 1):
    clusters_full = []
    for num, partition in enumerate(partitions):
        print("Process partition", num, file=sys.stderr)
        connectionGraph = nx.Graph()
        for i1 in range(len(partition)):
            for i2 in range(len(partition)):
                if i1 == i2 or gowdaDidayDistance(partition[i1], partition[i2], largestIndelSize) > maxDelta:
                    pass
                else:
                    connectionGraph.add_edge(i1, i2)
        #print(connectionGraph)
        clusters_indices = nx.find_cliques(connectionGraph)
        for cluster in clusters_indices:
            clusters_full.append([partition[index] for index in cluster])
    return clusters_full


def bronKerboschClustering(R, P, X, g):
    """taken from https://stackoverflow.com/a/13905694"""
    if not any((P, X)):
        yield R
    else:
        u = sorted([(numN(x, g), x) for x in P + X], key=lambda tup: tup[0], reverse=True)[0][1]
        u_n = N(u, g)
    for v in [elem for elem in P if elem not in u_n]:
        R_v = R + [v]
        P_v = [v1 for v1 in P if v1 in N(v, g)]
        X_v = [v1 for v1 in X if v1 in N(v, g)]
        for r in bronKerboschClustering(R_v, P_v, X_v, g):
            yield r
        P.remove(v)
        X.append(v)


def N(v, g):
    return [i for i, n_v in enumerate(g[v]) if n_v]


def numN(v, g):
    return sum(g[v])


def consolidateClusters(clusters):
    consolidatedClusters = []
    for cluster in clusters:
        l = len(cluster)
        starts = [member[1] for member in cluster]
        ends = [member[2] for member in cluster]
        consolidatedClusters.append( ( cluster[0][3], cluster[0][0], round(sum(starts) / l), round(sum(ends) / l), l, cluster) )
    return consolidatedClusters


options = parseArguments(sys.argv)

if not os.path.exists(options.temp_dir):
    print("ERROR: Given temp directory does not exist", file=sys.stderr)
    sys.exit()

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


options.fasta.close()

print("INFO: Temporary files written", file=sys.stderr)

if not os.path.exists(options.temp_dir + '/full_aln.chained.sorted.bam'):
    ngmlr = Popen(['/home/heller_d/bin/miniconda2/bin/ngmlr', '-t', '30', '-r', options.genome, '-q', options.temp_dir + '/full.fa', ], stdout=subprocess.PIPE)
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

print("INFO: Alignment finished", file=sys.stderr)

sam = pysam.AlignmentFile(options.temp_dir + '/full_aln.chained.sorted.bam', 'rb')
contigs = sam.references

insertion_output = open(options.temp_dir + '/ins.bed', 'w')
deletion_output = open(options.temp_dir + '/del.bed', 'w')
partition_output = open(options.temp_dir + '/partitions.bed', 'w')
indel_output = open(options.temp_dir + '/indels.bed', 'w')

largeIndels = []
largestIndelSize = -1

for contig in contigs:
    print("INFO: Searching for SVs in", contig, "..", file=sys.stderr)
    full_aln_dict = parseAlignmentFile(sam, contig)
    for read in full_aln_dict.keys():
        read_id = read.split("_")[1]

        #Read full read alignment
        try:
            full_aln = full_aln_dict[read][0]
            full_ref_start = full_aln.reference_start
            full_ref_end = full_aln.reference_end
            full_q_start = full_aln.query_alignment_start
            full_q_end = full_aln.query_alignment_end
            full_chrId = contig

            indels = findIndelsInCigarTuples(full_aln.cigartuples)
            for pos, l, typ in indels:
                largeIndels.append( (full_chrId, full_ref_start + pos, full_ref_start + pos + l, typ) )
                if l > largestIndelSize:
                    largestIndelSize = l
        except IndexError:
            full_aln = None

for chr, start, end, typ in largeIndels:
    print("{0}\t{1}\t{2}\t{3}".format(chr, start, end, typ), file = indel_output)

partitions = formPartitions(largeIndels)
print("Formed {0} partitions".format(len(partitions)), file=sys.stderr)
rawClusters = clustersFromPartitions2(partitions, largestIndelSize)
print("Subdivided partition into {0} clusters".format(len(rawClusters)), file=sys.stderr)
clusters = consolidateClusters(rawClusters)
consolidatedPartitions = consolidateClusters(partitions)

for typ, chr, start, end, support, members in consolidatedPartitions:
    print("{0}\t{1}\t{2}\t{3}\t{4}".format(chr, start, end, support, members), file = partition_output)

for typ, chr, start, end, support, members in clusters:
    if typ == 'ins':
        print("{0}\t{1}\t{2}\t{3}\t{4}".format(chr, start, end, support, members), file = insertion_output)
    if typ == 'del':
        print("{0}\t{1}\t{2}\t{3}\t{4}".format(chr, start, end, support, members), file = deletion_output)
