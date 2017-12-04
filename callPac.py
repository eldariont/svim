from __future__ import print_function

import sys
import argparse
import os
import re
from subprocess import call, Popen, PIPE
from collections import defaultdict
import pickle
import gzip

import pysam
from Bio import SeqIO

from callPacParams import callPacParams
from SVEvidence import EvidenceDeletion, EvidenceInsertion, EvidenceInversion, EvidenceTranslocation, EvidenceDuplicationTandem, EvidenceInsertionFrom
from SVCandidate import CandidateDeletion, CandidateInsertion, CandidateInversion, CandidateDuplicationTandem, CandidateDuplicationInterspersed

from callPacFull import analyze_full_read_indel, analyze_full_read_segments
from callPacTails import analyze_pair_of_read_tails
from callPacCluster import partition_and_cluster_unilocal, partition_and_cluster_bilocal
from callPacMerge import merge_insertions_from, merge_translocations

def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""callPac is a tool for accurate detection of structural variants (SVs).""")
    subparsers = parser.add_subparsers(help='subcommands', dest='sub')

    parser_bam = subparsers.add_parser('load', help='Load existing .obj file from working directory')
    parser_bam.add_argument('working_dir', type=str, help='working directory')

    parser_bam = subparsers.add_parser('alignment', help='Detect SVs from an existing alignment')
    parser_bam.add_argument('working_dir', type=str, help='working directory')
    parser_bam.add_argument('bam_file', type=argparse.FileType('r'), help='SAM/BAM file with aligned long reads (must be query-sorted)')
    parser_bam.add_argument('--skip_indel', action='store_true', help='disable indel part')
    parser_bam.add_argument('--skip_segment', action='store_true', help='disable segment part')

    parser_bam.add_argument('--tail_min_mapq', type=int, default=30, help='minimum mapping quality')

    parser_fasta = subparsers.add_parser('reads', help='Detect SVs from raw reads')
    parser_fasta.add_argument('working_dir', type=str, help='working directory')
    parser_fasta.add_argument('reads_file', type=str, help='Read file (FASTA, FASTQ, gzipped FASTA and FASTQ)')
    parser_fasta.add_argument('genome', type=str, default="/scratch/cluster/heller_d/genomes/hg19/hg19.fa", help='Reference genome file (FASTA)')

    parser_fasta.add_argument('--skip_kmer', action='store_true', help='disable kmer counting')
    parser_fasta.add_argument('--skip_indel', action='store_true', help='disable indel part')
    parser_fasta.add_argument('--skip_segment', action='store_true', help='disable segment part')
    parser_fasta.add_argument('--read_name', type=str, default="all", help='read name filter (default: all)')

    parser_fasta.add_argument('--cores', type=int, default=4, help='number of CPU cores to utilize for alignment')
    parser_fasta.add_argument('--tail_span', type=int, default=1000, help='length of read tails')
    parser_fasta.add_argument('--tail_min_mapq', type=int, default=30, help='minimum mapping quality')
    parser_fasta.add_argument('--tail_min_deviation', type=float, default=-0.1, help='minimum deviation')
    parser_fasta.add_argument('--tail_max_deviation', type=float, default=0.2, help='maximum deviation')
    parser_fasta.add_argument('--count_win_size', type=int, default=50, help='window size for k-mer counting')
    parser_fasta.add_argument('--count_k', type=int, default=7, help='k for k-mer counting')
    parser_fasta.add_argument('--count_band', type=float, default=0.5, help='band width')
    parser_fasta.add_argument('--stretch_threshold', type=int, default=7, help='z-score threshold')
    parser_fasta.add_argument('--stretch_tolerance', type=int, default=2, help='tolerance for stretch finding')
    parser_fasta.add_argument('--stretch_min_length', type=int, default=3, help='minimum stretch length')
    parser_fasta.add_argument('--path_constant_gap_cost', type=float, default=0, help='constant gap cost for path finding')
    parser_fasta.add_argument('--path_linear_gap_cost', type=float, default=0, help='linear gap cost for path finding')
    parser_fasta.add_argument('--path_convex_gap_cost', type=float, default=3, help='convex gap cost for path finding')
    parser_fasta.add_argument('--path_root_gap_cost', type=float, default=0, help='root gap cost for path finding')
    parser_fasta.add_argument('--path_tolerance', type=int, default=2, help='tolerance for overlapping segments')
    parser_fasta.add_argument('--align_costs_match', type=int, default=3, help='match cost for alignment')
    parser_fasta.add_argument('--align_costs_mismatch', type=int, default=-12, help='mismatch cost for alignment')
    parser_fasta.add_argument('--align_costs_gap', type=int, default=-12, help='gap cost for alignment')
    return parser.parse_args()


def parse_sam_file(sam, contig):
    """Parses a SAM file and returns a dict of reads (list of alignments for each read) for a given reference contig"""
    alns = sam.fetch(reference=contig)
    aln_dict = defaultdict(list)
    for aln in alns:
        aln_dict[aln.query_name].append(aln)
    return aln_dict


def guess_file_type(reads_path):
    if reads_path.endswith(".fa") or reads_path.endswith(".fasta") or reads_path.endswith(".FA"):
        print("INFO: Recognized reads file as FASTA format.", file=sys.stderr)
        return "fasta"
    elif reads_path.endswith(".fq") or reads_path.endswith(".fastq") or reads_path.endswith(".FQ"):
        print("INFO: Recognized reads file as FASTQ format.", file=sys.stderr)
        return "fastq"
    elif reads_path.endswith(".fa.gz") or reads_path.endswith(".fasta.gz") or reads_path.endswith(".fa.gzip") or reads_path.endswith(".fasta.gzip"):
        print("INFO: Recognized reads file as gzipped FASTA format.", file=sys.stderr)
        return "fasta_gzip"
    elif reads_path.endswith(".fq.gz") or reads_path.endswith(".fastq.gz") or reads_path.endswith(".fq.gzip") or reads_path.endswith(".fastq.gzip"):
        print("INFO: Recognized reads file as gzipped FASTQ format.", file=sys.stderr)
        return "fastq_gzip"
    else:
        print("ERROR: Unknown file ending of file {0}. Exiting.".format(reads_path), file=sys.stderr)
        return "unknown"

def create_tail_files(working_dir, reads_path, reads_type, span):
    """Create FASTA files with read tails and full reads if they do not exist."""
    if not os.path.exists(working_dir):
        print("ERROR: Given working directory does not exist", file=sys.stderr)
        sys.exit()

    if not os.path.exists(working_dir + '/left.fa'):
        left_file = open(working_dir + '/left.fa', 'w')
        write_left = True
    else:
        write_left = False
        print("WARNING: FASTA file for left tails exists. Skip", file=sys.stderr)

    if not os.path.exists(working_dir + '/right.fa'):
        right_file = open(working_dir + '/right.fa', 'w')
        write_right = True
    else:
        write_right = False
        print("WARNING: FASTA file for right tails exists. Skip", file=sys.stderr)

    if reads_type == "fasta_gzip" or reads_type == "fastq_gzip":
        full_reads_path = working_dir + '/full.fa'
        if not os.path.exists(working_dir + '/full.fa'):
            full_file = open(working_dir + '/full.fa', 'w')
            write_full = True
        else:
            write_full = False
            print("WARNING: FASTA file for full reads exists. Skip", file=sys.stderr)
    else:
        full_reads_path = reads_path
        write_full = False

    if write_left or write_right or write_full:
        if reads_type == "fasta" or reads_type == "fastq":
            reads_file = open(reads_path, "r")
        elif reads_type == "fasta_gzip" or reads_type == "fastq_gzip":
            reads_file = gzip.open(reads_path, "rb")
        
        if reads_type == "fasta" or reads_type == "fasta_gzip": 
            for line in reads_file:
                if line.startswith('>'):
                    read_name = line.strip()[1:]
                else:
                    sequence = line.strip()
                    if len(sequence) > 2000:
                        if write_left:
                            prefix = sequence[:span]
                            print(">" + read_name, file=left_file)
                            print(prefix, file=left_file)
                        if write_right:
                            suffix = sequence[-span:]
                            print(">" + read_name, file=right_file)
                            print(suffix, file=right_file)
                    if write_full:
                        print(">" + read_name, file=full_file)
                        print(sequence, file=full_file)
            reads_file.close()
        elif reads_type == "fastq" or reads_type == "fastq_gzip": 
            sequence_line_is_next = False
            for line in reads_file:
                if line.startswith('@'):
                    read_name = line.strip()[1:]
                    sequence_line_is_next = True
                elif sequence_line_is_next:
                    sequence_line_is_next = False
                    sequence = line.strip()
                    if len(sequence) > 2000:
                        if write_left:
                            prefix = sequence[:span]
                            print(">" + read_name, file=left_file)
                            print(prefix, file=left_file)
                        if write_right:
                            suffix = sequence[-span:]
                            print(">" + read_name, file=right_file)
                            print(suffix, file=right_file)
                    if write_full:
                        print(">" + read_name, file=full_file)
                        print(sequence, file=full_file)
            reads_file.close()
    
    if write_left:
        left_file.close()
    if write_right:
        right_file.close()
    if write_full:
        full_file.close()

    print("INFO: Read tail and full read files written", file=sys.stderr)
    return full_reads_path


def run_alignments(working_dir, genome, reads_path, cores):
    """Align full reads and read tails with NGM-LR and BWA MEM, respectively."""
    if not os.path.exists(working_dir + '/left_aln.rsorted.bam'):
        bwa = Popen(['/scratch/ngsvin/bin/bwa.kit/bwa',
                     'mem', '-x', 'pacbio', '-t', str(cores), genome, working_dir + '/left.fa'], stdout=PIPE)
        view = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'view', '-b', '-@', str(cores)], stdin=bwa.stdout, stdout=PIPE)
        sort = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'sort', '-@', str(cores), '-n', '-o', working_dir + '/left_aln.rsorted.bam'], stdin=view.stdout)
        sort.wait()
    else:
        print("WARNING: Alignment for left sequences exists. Skip", file=sys.stderr)

    if not os.path.exists(working_dir + '/right_aln.rsorted.bam'):
        bwa = Popen(['/scratch/ngsvin/bin/bwa.kit/bwa',
                     'mem', '-x', 'pacbio', '-t', str(cores), genome, working_dir + '/right.fa'], stdout=PIPE)
        view = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'view', '-b', '-@', str(cores)], stdin=bwa.stdout, stdout=PIPE)
        sort = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'sort', '-@', str(cores), '-n', '-o', working_dir + '/right_aln.rsorted.bam'], stdin=view.stdout)
        sort.wait()
    else:
        print("WARNING: Alignment for right sequences exists. Skip", file=sys.stderr)

    # Align full reads with NGM-LR
    if not os.path.exists(working_dir + '/full_aln.chained.rsorted.bam'):
        ngmlr = Popen(['/home/heller_d/bin/miniconda2/bin/ngmlr',
                       '-t', str(cores), '-r', genome, '-q', os.path.realpath(reads_path), ], stdout=PIPE)
        view = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'view', '-b', '-@', str(cores)], stdin=ngmlr.stdout, stdout=PIPE)
        sort = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'sort', '-n', '-@', str(cores), '-o', working_dir + '/full_aln.rsorted.bam'],
                     stdin=view.stdout)
        sort.wait()
        if call(['python', '/home/heller_d/bin/bamChain', working_dir + '/full_aln.rsorted.bam',
                 working_dir + '/full_aln.chained.bam', '--minmapq', '40']) != 0:
            print("ERROR: Calling bamchain on full sequences failed", file=sys.stderr)
        if call(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools', 'sort', '-n', '-@', str(cores),
                 working_dir + '/full_aln.chained.bam', '-o',
                 working_dir + '/full_aln.chained.rsorted.bam']) != 0:
            print("ERROR: Calling samtools sort on full sequences failed", file=sys.stderr)
    else:
        print("WARNING: Alignment for full sequences exists. Skip", file=sys.stderr)

    print("INFO: Alignment finished", file=sys.stderr)


def natural_representation(qname): 
    """Splits a read name into a tuple of strings and numbers. This facilitates the sort order applied by samtools -n
       See https://www.biostars.org/p/102735/"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', qname)]


def bam_iterator(bam):
    """Returns an iterator for the given SAM/BAM file (must be query-sorted). 
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


def analyze_read_tails(working_dir, genome, reads_path, reads_type, parameters):
    left_bam = pysam.AlignmentFile(working_dir + '/left_aln.rsorted.bam')
    right_bam = pysam.AlignmentFile(working_dir + '/right_aln.rsorted.bam')
    left_it = bam_iterator(left_bam)
    right_it = bam_iterator(right_bam)

    if reads_type == "fasta" or reads_type == "fasta_gzip" or reads_type == "fastq_gzip":
        reads = SeqIO.index_db(reads_path + ".idx", reads_path, "fasta")
    elif reads_type == "fastq":
        reads = SeqIO.index_db(reads_path + ".idx", reads_path, "fastq")
    reference =SeqIO.index_db(genome + ".idx", genome, "fasta")
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


def analyze_indel(bam_path, parameters):
    full_bam = pysam.AlignmentFile(bam_path)
    full_it = bam_iterator(full_bam)

    sv_evidences = []
    read_nr = 0

    while True:
        try:
            full_iterator_object = full_it.next()
            read_nr += 1
            if read_nr % 10000 == 0:
                print("INFO: Processed read", read_nr, file=sys.stderr)
            sv_evidences.extend(analyze_full_read_indel(full_iterator_object, full_bam, parameters))
        except StopIteration:
            break
        except KeyboardInterrupt:
            print('Execution interrupted by user. Stop detection and continue with clustering..')
            break
    return sv_evidences


def analyze_segments(bam_path, parameters):
    full_bam = pysam.AlignmentFile(bam_path)
    full_it = bam_iterator(full_bam)

    sv_evidences = []
    read_nr = 0

    while True:
        try:
            full_iterator_object = full_it.next()
            read_nr += 1
            if read_nr % 10000 == 0:
                print("INFO: Processed read", read_nr, file=sys.stderr)
            sv_evidences.extend(analyze_full_read_segments(full_iterator_object, full_bam, parameters))
        except StopIteration:
            break
        except KeyboardInterrupt:
            print('Execution interrupted by user. Stop detection and continue with clustering..')
            break
    return sv_evidences


def analyze_specific_read(working_dir, genome, reads_path, reads_type, parameters, read_name):
    left_bam = pysam.AlignmentFile(working_dir + '/left_aln.rsorted.bam')
    right_bam = pysam.AlignmentFile(working_dir + '/right_aln.rsorted.bam')
    left_it = bam_iterator(left_bam)
    right_it = bam_iterator(right_bam)

    if reads_type == "fasta" or reads_type == "fasta_gzip":
        reads = SeqIO.index_db(reads_path + ".idx", reads_path, "fasta")
    elif reads_type == "fastq" or reads_type == "fastq_gzip":
        reads = SeqIO.index_db(reads_path + ".idx", reads_path, "fastq")
    reference = SeqIO.index_db(genome + ".idx", genome, "fasta")
    print("INFO: Indexing reads and reference finished", file=sys.stderr)

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

    analyze_pair_of_read_tails(left_iterator_object, right_iterator_object, left_bam, right_bam, reads, reference, parameters)


def post_processing(sv_evidences, working_dir):
    deletion_evidences = [ev for ev in sv_evidences if ev.type == 'del']
    insertion_evidences = [ev for ev in sv_evidences if ev.type == 'ins']
    inversion_evidences = [ev for ev in sv_evidences if ev.type == 'inv']
    tandem_duplication_evidences = [ev for ev in sv_evidences if ev.type == 'dup']
    translocation_evidences = [ev for ev in sv_evidences if ev.type == 'tra']
    insertion_from_evidences = [ev for ev in sv_evidences if ev.type == 'ins_dup']

    print("INFO: Found {0}/{1}/{2}/{3}/{4}/{5} evidences for deletions, insertions, inversions, tandem duplications, translocations, and insertion_from, respectively.".format(
        len(deletion_evidences), len(insertion_evidences), len(inversion_evidences), len(tandem_duplication_evidences), len(translocation_evidences), len(insertion_from_evidences)), file=sys.stderr)
    
    # Cluster SV evidences
    print("INFO: Cluster deletion evidences..", file=sys.stderr)
    deletion_evidence_clusters = partition_and_cluster_unilocal(deletion_evidences)
    print("INFO: Cluster insertion evidences..", file=sys.stderr)
    insertion_evidence_clusters = partition_and_cluster_unilocal(insertion_evidences)
    print("INFO: Cluster inversion evidences..", file=sys.stderr)
    inversion_evidence_clusters = partition_and_cluster_unilocal(inversion_evidences)
    print("INFO: Cluster tandem duplication evidences..", file=sys.stderr)
    tandem_duplication_evidence_clusters = partition_and_cluster_bilocal(tandem_duplication_evidences)
    print("INFO: Cluster insertion evidences with source..", file=sys.stderr)
    insertion_from_evidence_clusters = partition_and_cluster_bilocal(insertion_from_evidences)

    # Print SV evidence clusters
    if not os.path.exists(working_dir + '/evidences'):
        os.mkdir(working_dir + '/evidences')
    deletion_evidence_output = open(working_dir + '/evidences/del.bed', 'w')
    insertion_evidence_output = open(working_dir + '/evidences/ins.bed', 'w')
    inversion_evidence_output = open(working_dir + '/evidences/inv.bed', 'w')
    tandem_duplication_evidence_source_output = open(working_dir + '/evidences/dup_tan_source.bed', 'w')
    tandem_duplication_evidence_dest_output = open(working_dir + '/evidences/dup_tan_dest.bed', 'w')
    translocation_evidence_output = open(working_dir + '/evidences/trans.bed', 'w')
    insertion_from_evidence_output = open(working_dir + '/evidences/ins_dup.bed', 'w')

    for cluster in deletion_evidence_clusters:
        print(cluster.get_bed_entry(), file=deletion_evidence_output)
    for cluster in insertion_evidence_clusters:
        print(cluster.get_bed_entry(), file=insertion_evidence_output)
    for cluster in inversion_evidence_clusters:
        print(cluster.get_bed_entry(), file=inversion_evidence_output)
    for cluster in tandem_duplication_evidence_clusters:
        bed_entries = cluster.get_bed_entries()
        print(bed_entries[0], file=tandem_duplication_evidence_source_output)
        print(bed_entries[1], file=tandem_duplication_evidence_dest_output)
    for cluster in insertion_from_evidence_clusters:
        bed_entries = cluster.get_bed_entries()
        print(bed_entries[0], file=insertion_from_evidence_output)
        print(bed_entries[1], file=insertion_from_evidence_output)

    
    # Merge evidences to candidates
    insertion_candidates = []
    int_duplication_candidates = []        

    #Merge insertions with source
    new_insertion_candidates,new_int_duplication_candidates = merge_insertions_from(insertion_from_evidence_clusters, deletion_evidence_clusters)
    insertion_candidates.extend(new_insertion_candidates)
    int_duplication_candidates.extend(new_int_duplication_candidates)

    #Merge translocation breakpoints
    new_insertion_candidates,new_int_duplication_candidates = merge_translocations(translocation_evidences, deletion_evidence_clusters, insertion_evidence_clusters)
    insertion_candidates.extend(new_insertion_candidates)
    int_duplication_candidates.extend(new_int_duplication_candidates)

    for translocation in translocation_evidences:
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(translocation.contig1, translocation.pos1, translocation.pos1+1, ">{0}:{1}".format(translocation.contig2, translocation.pos2), translocation.evidence, translocation.read), file=translocation_evidence_output)
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(translocation.contig2, translocation.pos2, translocation.pos2+1, ">{0}:{1}".format(translocation.contig1, translocation.pos1), translocation.evidence, translocation.read), file=translocation_evidence_output)

    if not os.path.exists(working_dir + '/candidates'):
        os.mkdir(working_dir + '/candidates')
    deletion_candidate_output = open(working_dir + '/candidates/del.bed', 'w')
    insertion_candidate_source_output = open(working_dir + '/candidates/cand_insertions_source.bed', 'w')
    insertion_candidate_dest_output = open(working_dir + '/candidates/cand_insertions_dest.bed', 'w')
    #inversion_candidate_output = open(working_dir + '/candidates/inv.bed', 'w')
    #tandem_duplication_candidate_output = open(working_dir + '/candidates/dup_tan.bed', 'w')
    interspersed_duplication_candidate_source_output = open(working_dir + '/candidates/cand_int_duplications_source.bed', 'w')
    interspersed_duplication_candidate_dest_output = open(working_dir + '/candidates/cand_int_duplications_dest.bed', 'w')
    #insertion_from_candidate_output = open(working_dir + '/candidates/ins_dup.bed', 'w')

    for candidate in insertion_candidates:
        bed_entries = candidate.get_bed_entries()
        print(bed_entries[0], file=insertion_candidate_source_output)
        print(bed_entries[1], file=insertion_candidate_dest_output)
    for candidate in int_duplication_candidates:
        bed_entries = candidate.get_bed_entries()
        print(bed_entries[0], file=interspersed_duplication_candidate_source_output)
        print(bed_entries[1], file=interspersed_duplication_candidate_dest_output)
    for cluster in deletion_evidence_clusters:
        print(cluster.get_bed_entry(), file=deletion_candidate_output)


def main():
    #Fetch command-line options and set parameters accordingly
    options = parse_arguments()
    parameters = callPacParams()
    parameters.set_with_options(options)

    #Search for SV evidences
    if options.sub == 'load':
        print("INFO: Load sv_evidences.obj with SV evidences..", file=sys.stderr)
        evidences_file = open(options.working_dir + '/sv_evidences.obj', 'r')
        sv_evidences = pickle.load(evidences_file)
        evidences_file.close()
    elif options.sub == 'reads':
        reads_type = guess_file_type(options.reads_file)
        if reads_type == "unknown":
            return
        full_reads_path = create_tail_files(options.working_dir, options.reads_file, reads_type, parameters.tail_span)
        run_alignments(options.working_dir, options.genome, full_reads_path, options.cores)
        if options.read_name != "all":
            analyze_specific_read(options.working_dir, options.genome, full_reads_path, reads_type, parameters, options.read_name)
            return
        sv_evidences = []
        if not options.skip_kmer:
            sv_evidences.extend(analyze_read_tails(options.working_dir, options.genome, full_reads_path, reads_type, parameters))
        if not options.skip_indel:
            sv_evidences.extend(analyze_indel(options.working_dir + '/full_aln.chained.rsorted.bam', parameters))
        if not options.skip_segment:
            sv_evidences.extend(analyze_segments(options.working_dir + '/full_aln.chained.rsorted.bam', parameters))
        evidences_file = open(options.working_dir + '/sv_evidences.obj', 'w')
        pickle.dump(sv_evidences, evidences_file) 
        evidences_file.close()
    elif options.sub == 'alignment':
        sv_evidences = []
        if not options.skip_indel:
            sv_evidences.extend(analyze_indel(options.bam_file.name, parameters))
        if not options.skip_segment:
            sv_evidences.extend(analyze_segments(options.bam_file.name, parameters))
        evidences_file = open(options.working_dir + '/sv_evidences.obj', 'w')
        pickle.dump(sv_evidences, evidences_file)
        evidences_file.close()

    #Post-process SV evidences
    post_processing(sv_evidences, options.working_dir)    

if __name__ == "__main__":
    sys.exit(main())
