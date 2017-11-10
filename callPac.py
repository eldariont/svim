from __future__ import print_function

import sys
import argparse
import os
import re
from subprocess import call, Popen, PIPE
from collections import defaultdict

import pysam
from Bio import SeqIO

from callPacParams import callPacParams
from SVEvidence import EvidenceDeletion, EvidenceInsertion, EvidenceInversion, EvidenceTranslocation, EvidenceDuplicationTandem, EvidenceInsertionFrom
from SVCandidate import CandidateDeletion, CandidateInsertion, CandidateInversion, CandidateDuplicationTandem

from callPacFull import analyze_full_read
from callPacTails import analyze_pair_of_read_tails
from callPacCluster import partition_and_cluster


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
    parser.add_argument('--read_name', type=str, default="all", help='read name filter (default: all)')
    return parser.parse_args()


def parse_sam_file(sam, contig):
    """Parses a SAM file and returns a dict of reads (list of alignments for each read) for a given reference contig"""
    alns = sam.fetch(reference=contig)
    aln_dict = defaultdict(list)
    for aln in alns:
        aln_dict[aln.query_name].append(aln)
    return aln_dict


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
    """Splits a read name into a tuple of strings and numbers. This facilitates the sort order applied by samtools -n
       See https://www.biostars.org/p/102735/"""
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
                sv_evidences.extend(analyze_pair_of_read_tails(left_iterator_object, right_iterator_object, left_bam, right_bam, reads, reference, parameters))
                left_iterator_object = left_it.next()
        except StopIteration:
            break
        except KeyboardInterrupt:
            print('Execution interrupted by user. Stop detection and continue with clustering..')
            break
    return sv_evidences


def analyze_full_reads(temp_dir, genome, fasta, parameters):
    full_bam = pysam.AlignmentFile(temp_dir + '/full_aln.chained.rsorted.bam')
    full_it = bam_iterator(full_bam)

    sv_evidences = []

    while True:
        try:
            full_iterator_object = full_it.next()
            if int(full_iterator_object[0].split("_")[1]) % 1000 == 0:
                print("INFO: Processed read", full_iterator_object[0].split("_")[1], file=sys.stderr)
            sv_evidences.extend(analyze_full_read(full_iterator_object, full_bam, parameters))
        except StopIteration:
            break
        except KeyboardInterrupt:
            print('Execution interrupted by user. Stop detection and continue with clustering..')
            break
    return sv_evidences


def analyze_specific_read(temp_dir, genome, fasta, parameters, read_name):
    left_bam = pysam.AlignmentFile(temp_dir + '/left_aln.rsorted.bam')
    right_bam = pysam.AlignmentFile(temp_dir + '/right_aln.rsorted.bam')
    left_it = bam_iterator(left_bam)
    right_it = bam_iterator(right_bam)

    reads = SeqIO.index(fasta.name, "fasta")
    reference = SeqIO.index(genome, "fasta")
    print("INFO: Indexing reads and reference finished", file=sys.stderr)

    sv_evidences = []

    left_iterator_object = left_it.next()
    while left_iterator_object[0] != read_name:
        try:
            left_iterator_object = left_it.next()
        except StopIteration:
            break
    right_iterator_object = right_it.next()
    while right_iterator_object[0] != read_name:
        try:
            right_iterator_object = right_it.next()
        except StopIteration:
            break

    result = analyze_pair_of_read_tails(left_iterator_object, right_iterator_object, left_bam, right_bam, reads, reference, parameters)


def read_evidences_file(path):
    raw_file = open(path, 'r')
    sv_evidences = []
    for line in raw_file:
        contig, start, end, name, score, read = line.strip().split("\t")
        name_fields = name.split(";")
        if name_fields[0] == "del":
            sv_evidences.append(EvidenceDeletion(contig, int(start), int(end), name_fields[1], read))
        elif name_fields[0] == "ins":
            sv_evidences.append(EvidenceInsertion(contig, int(start), int(end), name_fields[1], read))
        elif name_fields[0] == "inv":
            sv_evidences.append(EvidenceInversion(contig, int(start), int(end), name_fields[1], read))
        elif name_fields[0] == "dup":
            sv_evidences.append(EvidenceDuplicationTandem(contig, int(start), int(end), name_fields[1], name_fields[2],read))
        elif name_fields[0] == "tra":
            contig2, pos2 = name_fields[1].split(":")
            sv_evidences.append(EvidenceTranslocation(contig, int(start), contig2, int(pos2), name_fields[2], read))
        elif name_fields[0] == "ins_dup":
            source_contig = name_fields[1].split(":")[0]
            source_start, source_end = name_fields[1].split(":")[1].split("-")
            sv_evidences.append(EvidenceInsertionFrom(source_contig, int(source_start), int(source_end), contig, int(start), name_fields[2], read))
    raw_file.close()
    return sv_evidences


def main():
    options = parse_arguments()
    parameters = callPacParams()
    parameters.set_with_options(options)

    # Run SV search only if raw SV results do not exist
    if not os.path.exists(options.temp_dir + '/sv_evidences.tsv'):
        create_temp_files(options.temp_dir, options.fasta, parameters.tail_span)
        run_alignments(options.temp_dir, options.genome, options.fasta)
        if options.read_name != "all":
            analyze_specific_read(options.temp_dir, options.genome, options.fasta, parameters, options.read_name)
            return
        else:
            sv_evidences = analyze_read_tails(options.temp_dir, options.genome, options.fasta, parameters)
            sv_evidences.extend(analyze_full_reads(options.temp_dir, options.genome, options.fasta, parameters))

        raw_file = open(options.temp_dir + '/sv_evidences.tsv', 'w')
        for sv_evidence in sv_evidences:
            print(sv_evidence.get_bed_entry(), file=raw_file)
        raw_file.close()
    else:
        print("WARNING: Result file with SV evidences already exists. Skip", file=sys.stderr)
        sv_evidences = read_evidences_file(options.temp_dir + '/sv_evidences.tsv')

    deletion_evidences = [ev for ev in sv_evidences if ev.type == 'del']
    insertion_evidences = [ev for ev in sv_evidences if ev.type == 'ins']
    inversion_evidences = [ev for ev in sv_evidences if ev.type == 'inv']
    tandem_duplication_evidences = [ev for ev in sv_evidences if ev.type == 'dup']
    translocation_evidences = [ev for ev in sv_evidences if ev.type == 'tra']
    insertion_from_evidences = [ev for ev in sv_evidences if ev.type == 'ins_dup']

    print("INFO: Found {0}/{1}/{2}/{3}/{4}/{5} evidences for deletions, insertions, inversions, tandem duplications, translocations, and insertion_from, respectively.".format(
        len(deletion_evidences), len(insertion_evidences), len(inversion_evidences), len(tandem_duplication_evidences), len(translocation_evidences), len(insertion_from_evidences)), file=sys.stderr)
    
    # Cluster raw SVs
    deletion_evidence_clusters = partition_and_cluster(deletion_evidences)
    insertion_evidence_clusters = partition_and_cluster(insertion_evidences)
    inversion_evidence_clusters = partition_and_cluster(inversion_evidences)
    tandem_duplication_evidence_clusters = partition_and_cluster(tandem_duplication_evidences)
    insertion_from_evidence_clusters = partition_and_cluster(insertion_from_evidences)

    # Print output
    deletion_output = open(options.temp_dir + '/del.bed', 'w')
    insertion_output = open(options.temp_dir + '/ins.bed', 'w')
    inversion_output = open(options.temp_dir + '/inv.bed', 'w')
    tandem_duplication_output = open(options.temp_dir + '/dup_tan.bed', 'w')
    translocation_output = open(options.temp_dir + '/trans.bed', 'w')
    insertion_from_output = open(options.temp_dir + '/ins_dup.bed', 'w')

    for cluster in deletion_evidence_clusters:
        print(cluster.get_bed_entry(), file=deletion_output)
    for cluster in insertion_evidence_clusters:
        print(cluster.get_bed_entry(), file=insertion_output)
    for cluster in inversion_evidence_clusters:
        print(cluster.get_bed_entry(), file=inversion_output)
    for cluster in tandem_duplication_evidence_clusters:
        print(cluster.get_bed_entry(), file=tandem_duplication_output)
    for cluster in insertion_from_evidence_clusters:
        print(cluster.get_bed_entry(), file=insertion_from_output)

    for translocation in translocation_evidences:
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(translocation.contig1, translocation.pos1, translocation.pos1+1, ">{0}:{1}".format(translocation.contig2, translocation.pos2), translocation.evidence, translocation.read), file=translocation_output)
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(translocation.contig2, translocation.pos2, translocation.pos2+1, ">{0}:{1}".format(translocation.contig1, translocation.pos1), translocation.evidence, translocation.read), file=translocation_output)


if __name__ == "__main__":
    sys.exit(main())
