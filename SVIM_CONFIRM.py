from __future__ import print_function

__version__ = '0.2'
__author__ = 'David Heller'

import sys
import argparse
import os
import pickle
import gzip
import logging
import ConfigParser

from subprocess import Popen, PIPE, call
from collections import defaultdict
from time import strftime, localtime
from math import pow, sqrt

import pysam
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Bio import SeqIO
from matplotlib.backends.backend_pdf import PdfPages

from SVIM_clustering import form_partitions, partition_and_cluster_candidates
from SVCandidate import CandidateInversion
from SVIM_merging import merge_insertions_from, merge_translocations_at_deletions, merge_translocations_at_insertions
from SVIM_readtails import confirm_del, confirm_ins, confirm_inv


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""SVIM (pronounced SWIM) is a structural variant caller for long reads. 
It combines full alignment analysis, split-read mapping, and read-tail mapping to 
distinguish five classes of structural variants. SVIM discriminates between similar 
SV classes such as interspersed duplications and cut&paste insertions and is unique 
in its capability of extracting both the genomic origin and destination of insertions 
and duplications.

SVIM consists of two programs SWIM-COLLECT and SWIM-CONFIRM. You are running SWIM-CONFIRM 
which confirms SV evidences using read-tail mapping and classifies them into distinct 
SV types. To confirm SV evidences, you need to supply both the --reads and --genome 
arguments.""")
    parser.add_argument('--version', '-v', action='version', version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('working_dir', type=str, help='working directory')
    parser.add_argument('--config', type=str, default="{0}/default_config.cfg".format(os.path.dirname(os.path.realpath(__file__))), help='configuration file, default: {0}/default_config.cfg'.format(os.path.dirname(os.path.realpath(__file__))))
    parser.add_argument('--obj_file', '-i', type=argparse.FileType('r'), help='Path of .obj file to load (default: working_dir/sv_evidences.obj')
    parser.add_argument('--reads', type=str, help='Read file (FASTA, FASTQ, gzipped FASTA and FASTQ)')
    parser.add_argument('--genome', type=str, help='Reference genome file (FASTA)')
    parser.add_argument('--debug_confirm', action='store_true', help='print dot plots when confirming SV evidence clusters')
    return parser.parse_args()


def read_parameters(options):
    config = ConfigParser.RawConfigParser()
    config.read(options.config)

    parameters = dict()
    parameters["cores"] = config.getint("alignment", "cores")

    parameters["min_mapq"] = config.getint("detection", "min_mapq")
    parameters["max_sv_size"] = config.getint("detection", "max_sv_size")

    parameters["min_length"] = config.getint("split read", "min_length")
    parameters["max_segment_gap_tolerance"] = config.getint("split read", "max_segment_gap_tolerance")
    parameters["max_deletion_size"] = config.getint("split read", "max_deletion_size")
    parameters["segment_overlap_tolerance"] = config.getint("split read", "segment_overlap_tolerance")

    parameters["partition_max_distance"] = config.getint("clustering", "partition_max_distance")
    parameters["cluster_max_distance"] = config.getfloat("clustering", "cluster_max_distance")

    parameters["del_ins_dup_max_distance"] = config.getfloat("merging", "del_ins_dup_max_distance")
    parameters["trans_destination_partition_max_distance"] = config.getint("merging", "trans_destination_partition_max_distance")
    parameters["trans_partition_max_distance"] = config.getint("merging", "trans_partition_max_distance")
    parameters["trans_sv_max_distance"] = config.getint("merging", "trans_sv_max_distance")

    parameters["tail_span"] = config.getint("confirmation", "tail_span")
    parameters["tail_min_deviation"] = config.getfloat("confirmation", "tail_min_deviation")
    parameters["tail_max_deviation"] = config.getfloat("confirmation", "tail_max_deviation")
    parameters["count_win_size"] = config.getint("confirmation", "count_win_size")
    parameters["count_k"] = config.getint("confirmation", "count_k")
    parameters["count_band"] = config.getfloat("confirmation", "count_band")
    parameters["stretch_threshold"] = config.getint("confirmation", "stretch_threshold")
    parameters["stretch_tolerance"] = config.getint("confirmation", "stretch_tolerance")
    parameters["stretch_min_length"] = config.getint("confirmation", "stretch_min_length")
    parameters["path_constant_gap_cost"] = config.getint("confirmation", "path_constant_gap_cost")
    parameters["path_linear_gap_cost"] = config.getint("confirmation", "path_linear_gap_cost")
    parameters["path_convex_gap_cost"] = config.getint("confirmation", "path_convex_gap_cost")
    parameters["path_root_gap_cost"] = config.getint("confirmation", "path_root_gap_cost")
    parameters["path_tolerance"] = config.getint("confirmation", "path_tolerance")
    parameters["align_costs"] = (config.getint("confirmation", "align_costs_match"), config.getint("confirmation", "align_costs_mismatch"), config.getint("confirmation", "align_costs_gap"))

    try:
        parameters["skip_indel"] =  options.skip_indel
    except AttributeError:
        parameters["skip_indel"] =  False
    try:
        parameters["skip_segment"] =  options.skip_segment
    except AttributeError:
        parameters["skip_segment"] =  False
    try:
        parameters["skip_confirm"] =  options.skip_confirm
    except AttributeError:
        parameters["skip_confirm"] =  False
    try:
        parameters["debug_confirm"] =  options.debug_confirm
    except AttributeError:
        parameters["debug_confirm"] =  False

    return parameters


def guess_file_type(reads_path):
    if reads_path.endswith(".fa") or reads_path.endswith(".fasta") or reads_path.endswith(".FA"):
        logging.info("Recognized reads file as FASTA format.")
        return "fasta"
    elif reads_path.endswith(".fq") or reads_path.endswith(".fastq") or reads_path.endswith(".FQ"):
        logging.info("Recognized reads file as FASTQ format.")
        return "fastq"
    elif reads_path.endswith(".fa.gz") or reads_path.endswith(".fasta.gz") or reads_path.endswith(".fa.gzip") or reads_path.endswith(".fasta.gzip"):
        logging.info("Recognized reads file as gzipped FASTA format.")
        return "fasta_gzip"
    elif reads_path.endswith(".fq.gz") or reads_path.endswith(".fastq.gz") or reads_path.endswith(".fq.gzip") or reads_path.endswith(".fastq.gzip"):
        logging.info("Recognized reads file as gzipped FASTQ format.")
        return "fastq_gzip"
    elif reads_path.endswith(".fa.fn"):
        logging.info("Recognized reads file as FASTA file list format.")
        return "list" 
    else:
        logging.error("Unknown file ending of file {0}. Exiting.".format(reads_path))
        return "unknown"


def create_tail_files(working_dir, reads_path, reads_type, span):
    """Create FASTA files with read tails and full reads if they do not exist."""
    if not os.path.exists(working_dir):
        logging.error("Given working directory does not exist")
        sys.exit()

    reads_file_prefix = os.path.splitext(os.path.basename(reads_path))[0]

    if not os.path.exists("{0}/{1}_left.fa".format(working_dir, reads_file_prefix)):
        left_file = open("{0}/{1}_left.fa".format(working_dir, reads_file_prefix), 'w')
        write_left = True
    else:
        write_left = False
        logging.warning("FASTA file for left tails exists. Skip")

    if not os.path.exists("{0}/{1}_right.fa".format(working_dir, reads_file_prefix)):
        right_file = open("{0}/{1}_right.fa".format(working_dir, reads_file_prefix), 'w')
        write_right = True
    else:
        write_right = False
        logging.warning("FASTA file for right tails exists. Skip")

    if reads_type == "fasta_gzip" or reads_type == "fastq_gzip":
        full_reads_path = "{0}/{1}.fa".format(working_dir, reads_file_prefix)
        if not os.path.exists(full_reads_path):
            full_file = open(full_reads_path, 'w')
            write_full = True
        else:
            write_full = False
            logging.warning("FASTA file for full reads exists. Skip")
    else:
        full_reads_path = reads_path
        write_full = False

    if write_left or write_right or write_full:
        if reads_type == "fasta" or reads_type == "fastq":
            reads_file = open(reads_path, "r")
        elif reads_type == "fasta_gzip" or reads_type == "fastq_gzip":
            reads_file = gzip.open(reads_path, "rb")
        
        if reads_type == "fasta" or reads_type == "fasta_gzip":
            sequence = ""
            for line in reads_file:
                if line.startswith('>'):
                    if sequence != "":
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
                    read_name = line.strip()[1:]
                    sequence = ""
                else:
                    sequence += line.strip()
            if sequence != "":
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

    logging.info("Read tail and full read files written")
    return full_reads_path


def run_tail_alignments(working_dir, genome, reads_path, cores):
    """Align read tails with NGM-LR."""
    reads_file_prefix = os.path.splitext(os.path.basename(reads_path))[0]
    left_fa = "{0}/{1}_left.fa".format(working_dir, reads_file_prefix)
    left_aln = "{0}/{1}_left_aln.coordsorted.bam".format(working_dir, reads_file_prefix)
    left_aln_index = "{0}/{1}_left_aln.coordsorted.bam.bai".format(working_dir, reads_file_prefix)
    right_fa = "{0}/{1}_right.fa".format(working_dir, reads_file_prefix)
    right_aln = "{0}/{1}_right_aln.coordsorted.bam".format(working_dir, reads_file_prefix)
    right_aln_index = "{0}/{1}_right_aln.coordsorted.bam.bai".format(working_dir, reads_file_prefix)

    if not os.path.exists(left_aln):
        bwa = Popen(['ngmlr',
                     '-t', str(cores), '-r', genome, '-q', left_fa], stdout=PIPE)
        view = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'view', '-b', '-@', str(cores)], stdin=bwa.stdout, stdout=PIPE)
        sort = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'sort', '-@', str(cores), '-o', left_aln], stdin=view.stdout)
        sort.wait()
    else:
        logging.warning("Alignment for left sequences exists. Skip")

    if not os.path.exists(left_aln_index):
        call(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'index', left_aln])
    else:
        logging.warning("Alignment index for left sequences exists. Skip")

    if not os.path.exists(right_aln):
        bwa = Popen(['ngmlr',
                     '-t', str(cores), '-r', genome, '-q', right_fa], stdout=PIPE)
        view = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'view', '-b', '-@', str(cores)], stdin=bwa.stdout, stdout=PIPE)
        sort = Popen(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'sort', '-@', str(cores), '-o', right_aln], stdin=view.stdout)
        sort.wait()
    else:
        logging.warning("Alignment for right sequences exists. Skip")

    if not os.path.exists(right_aln_index):
        call(['/scratch/ngsvin/bin/samtools/samtools-1.3.1/samtools',
                      'index', right_aln])
    else:
        logging.warning("Alignment index for right sequences exists. Skip")

    logging.info("Tail alignments finished")


def calculate_score_inversion(direction_counts, inversion_length, successful_confirmations, total_confirmations, parameters):
    left_evidences = direction_counts[0] + direction_counts[1]
    right_evidences = direction_counts[2] + direction_counts[3]
    valid_suppl_evidences = min(left_evidences, right_evidences) + direction_counts[4]
    if inversion_length > parameters["max_sv_size"]:
        return 0
    else:
        if total_confirmations > 0:
            confirmation_rate = successful_confirmations / float(total_confirmations)
            if confirmation_rate > 0.5:
                return valid_suppl_evidences + confirmation_rate * 20
            elif confirmation_rate < 0.3:
                return 0
            else:
                return valid_suppl_evidences
        else:
            return valid_suppl_evidences


def cluster_sv_candidates(insertion_candidates, int_duplication_candidates, parameters):
    """Takes a list of SVCandidates and splits them up by type. The SVCandidates of each type are clustered and returned as a tuple of
    (deletion_evidence_clusters, insertion_evidence_clusters, inversion_evidence_clusters, tandem_duplication_evidence_clusters, insertion_from_evidence_clusters, completed_translocation_evidences)."""

    logging.info("Found {0}/{1} candidates for insertions and interspersed duplications.".format(
        len(insertion_candidates), len(int_duplication_candidates)))

    logging.info("Cluster insertion candidates:")
    final_insertion_candidates = partition_and_cluster_candidates(insertion_candidates, parameters)
    logging.info("Cluster interspersed duplication candidates:")
    final_int_duplication_candidates = partition_and_cluster_candidates(int_duplication_candidates, parameters)

    return (final_insertion_candidates, final_int_duplication_candidates)


def write_candidates(working_dir, candidates):
    insertion_candidates, int_duplication_candidates, inversion_candidates = candidates

    if not os.path.exists(working_dir + '/candidates'):
        os.mkdir(working_dir + '/candidates')
    #deletion_candidate_output = open(working_dir + '/candidates/candidates_deletions.bed', 'w')
    insertion_candidate_source_output = open(working_dir + '/candidates/candidates_insertions_source.bed', 'w')
    insertion_candidate_dest_output = open(working_dir + '/candidates/candidates_insertions_dest.bed', 'w')
    inversion_candidate_output = open(working_dir + '/candidates/candidates_inversions.bed', 'w')
    # tandem_duplication_candidate_output = open(working_dir + '/candidates/dup_tan.bed', 'w')
    interspersed_duplication_candidate_source_output = open(working_dir + '/candidates/candidates_int_duplications_source.bed', 'w')
    interspersed_duplication_candidate_dest_output = open(working_dir + '/candidates/candidates_int_duplications_dest.bed', 'w')

    for candidate in insertion_candidates:
        bed_entries = candidate.get_bed_entries()
        print(bed_entries[0], file=insertion_candidate_source_output)
        print(bed_entries[1], file=insertion_candidate_dest_output)
    for candidate in int_duplication_candidates:
        bed_entries = candidate.get_bed_entries()
        print(bed_entries[0], file=interspersed_duplication_candidate_source_output)
        print(bed_entries[1], file=interspersed_duplication_candidate_dest_output)
    for candidate in inversion_candidates:
        print(candidate.get_bed_entry(), file=inversion_candidate_output)

    insertion_candidate_source_output.close()
    insertion_candidate_dest_output.close()
    inversion_candidate_output.close()
    interspersed_duplication_candidate_source_output.close()
    interspersed_duplication_candidate_dest_output.close()


def write_final_vcf(working_dir, genome, insertion_candidates, int_duplication_candidates, inversion_candidates, deletion_evidence_clusters, tandem_duplication_evidence_clusters):
    vcf_output = open(working_dir + '/final_results.vcf', 'w')

    # Write header lines
    print("##fileformat=VCFv4.3", file=vcf_output)
    print("##source=SVIMV{0}".format(__version__), file=vcf_output)
    print("##reference={0}".format(genome), file=vcf_output)
    print("##ALT=<ID=DEL,Description=\"Deletion\">", file=vcf_output)
    print("##ALT=<ID=INV,Description=\"Inversion\">", file=vcf_output)
    print("##ALT=<ID=DUP,Description=\"Duplication\">", file=vcf_output)
    print("##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">", file=vcf_output)
    print("##ALT=<ID=DUP:INT,Description=\"Interspersed Duplication\">", file=vcf_output)
    print("##ALT=<ID=INS,Description=\"Insertion\">", file=vcf_output)
    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">", file=vcf_output)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">", file=vcf_output)
    print("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">", file=vcf_output)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=vcf_output)

    vcf_entries = []
    for cluster in deletion_evidence_clusters:
        vcf_entries.append((cluster.get_source(), cluster.get_vcf_entry()))
    for cluster in tandem_duplication_evidence_clusters:
        vcf_entries.append((cluster.get_source(), cluster.get_vcf_entry()))
    for candidate in insertion_candidates:
        vcf_entries.append((candidate.get_destination(), candidate.get_vcf_entry()))
    for candidate in int_duplication_candidates:
        vcf_entries.append((candidate.get_destination(), candidate.get_vcf_entry()))
    for candidate in inversion_candidates:
        vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry()))

    # Sort and write entries to VCF
    for source, entry in sorted(vcf_entries, key=lambda pair: pair[0]):
        print(entry, file=vcf_output)

    vcf_output.close()


def main():
    # Fetch command-line options and configuration file values and set parameters accordingly
    options = parse_arguments()
    parameters = read_parameters(options)

    # Set up logging
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-7.7s]  %(message)s")
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    fileHandler = logging.FileHandler("{0}/SWIM-CONFIRM_{1}.log".format(options.working_dir, strftime("%y%m%d_%H%M%S", localtime())), mode="w")
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    logging.info("****************** Start SVIM-CONFIRM, version {0} ******************".format(__version__))
    logging.info("CMD: python {0}".format(" ".join(sys.argv)))
    logging.info("WORKING DIR: {0}".format(os.path.abspath(options.working_dir)))

    if options.obj_file:
        logging.info("INPUT: {0}".format(os.path.abspath(options.obj_file.name)))
        evidences_file = options.obj_file
    else:
        logging.info("INPUT: {0}".format(os.path.abspath(options.working_dir + '/sv_evidences.obj')))
        evidences_file = open(options.working_dir + '/sv_evidences.obj', 'r')
    evidence_clusters = pickle.load(evidences_file)
    deletion_evidence_clusters, insertion_evidence_clusters, inversion_evidence_clusters, tandem_duplication_evidence_clusters, insertion_from_evidence_clusters, completed_translocations = evidence_clusters
    evidences_file.close()

    if options.reads and options.genome:
        reads_type = guess_file_type(options.reads)
        full_reads_path = create_tail_files(options.working_dir, options.reads, reads_type, parameters["tail_span"])
        run_tail_alignments(options.working_dir, options.genome, full_reads_path, parameters["cores"])
        parameters["skip_confirm"] = False
    else:
        parameters["skip_confirm"] = True
        logging.info("Skipping confirmation with read tails because reads or genome are missing.")

    ####################
    # Confirm clusters #
    ####################
    if not parameters["skip_confirm"]:
        del_confirmation_threshold = np.percentile([cluster.score for cluster in deletion_evidence_clusters], 25, interpolation='higher')
        ins_confirmation_threshold = np.percentile([cluster.score for cluster in insertion_evidence_clusters], 50, interpolation='higher')
        inv_confirmation_threshold = np.percentile([cluster.score for cluster in inversion_evidence_clusters], 25, interpolation='higher')

        logging.info("Confirming deleted, inserted and inverted regions with a score smaller than {0}, {1} and {2}, respectively.".format(del_confirmation_threshold, ins_confirmation_threshold, inv_confirmation_threshold))

        reads_file_prefix = os.path.splitext(os.path.basename(full_reads_path))[0]
        left_aln = "{0}/{1}_left_aln.coordsorted.bam".format(options.working_dir, reads_file_prefix)
        right_aln = "{0}/{1}_right_aln.coordsorted.bam".format(options.working_dir, reads_file_prefix)
        left_bam = pysam.AlignmentFile(left_aln)
        right_bam = pysam.AlignmentFile(right_aln)

        reads = SeqIO.index_db(full_reads_path + ".idx", full_reads_path, "fasta")
        reference =SeqIO.index_db(options.genome + ".idx", options.genome, "fasta")
        logging.info("Indexing reads and reference finished")

        if parameters["debug_confirm"]:
            del_pdf = PdfPages('dotplots_deletions.pdf')
            ins_pdf = PdfPages('dotplots_insertions.pdf')

        for del_cluster in deletion_evidence_clusters:
            if del_cluster.score <= del_confirmation_threshold:
                if parameters["debug_confirm"]:
                    fig = plt.figure()
                    fig.suptitle('Deleted region cluster (score {0}) {1}:{2}-{3}'.format(del_cluster.score, *del_cluster.get_source()), fontsize=10)
                successful_confirmations, total_confirmations = confirm_del(left_bam, right_bam, del_cluster, reads, reference, parameters)
                if total_confirmations > 0:
                    confirmation_rate = successful_confirmations / float(total_confirmations)
                    if confirmation_rate > 0.5:
                        del_cluster.score += int(confirmation_rate * 20)
                    elif confirmation_rate < 0.3 and del_cluster.end - del_cluster.start > parameters["count_win_size"] * 3:
                        del_cluster.score = 0
                if parameters["debug_confirm"]:
                    del_pdf.savefig(fig)
                    plt.close(fig)
        if parameters["debug_confirm"]:
            del_pdf.close()

        for ins_cluster in insertion_evidence_clusters:
            if ins_cluster.score <= ins_confirmation_threshold:
                if parameters["debug_confirm"]:
                    fig = plt.figure()
                    fig.suptitle('Inserted region cluster (score {0}) {1}:{2}-{3}'.format(ins_cluster.score, *ins_cluster.get_source()), fontsize=10)
                successful_confirmations, total_confirmations = confirm_ins(left_bam, right_bam, ins_cluster, reads, reference, parameters)
                if total_confirmations > 0:
                    confirmation_rate = successful_confirmations / float(total_confirmations)
                    if confirmation_rate > 0.5:
                        ins_cluster.score += int(confirmation_rate * 20)
                    elif confirmation_rate < 0.3 and ins_cluster.end - ins_cluster.start > parameters["count_win_size"] * 3:
                        ins_cluster.score = 0
                if parameters["debug_confirm"]:
                    ins_pdf.savefig(fig)
                    plt.close(fig)
        if parameters["debug_confirm"]:
            ins_pdf.close()

    inversion_candidates = []
    for inv_cluster in inversion_evidence_clusters:
        directions = [ev.direction for ev in inv_cluster.members]
        direction_counts = [0, 0, 0, 0, 0]
        for direction in directions:
            if direction == "left_fwd": direction_counts[0] += 1
            if direction == "left_rev": direction_counts[1] += 1
            if direction == "right_fwd": direction_counts[2] += 1
            if direction == "right_rev": direction_counts[3] += 1
            if direction == "all": direction_counts[4] += 1
        contig, start, end = inv_cluster.get_source()

        if not parameters["skip_confirm"] and inv_cluster.score <= inv_confirmation_threshold:
            successful_confirmations, total_confirmations = confirm_inv(left_bam, right_bam, inv_cluster, reads, reference, parameters)
            score = calculate_score_inversion(direction_counts, end - start, successful_confirmations, total_confirmations, parameters)
        else:
            score = calculate_score_inversion(direction_counts, end - start, 0, 0, parameters)
        inversion_candidates.append(CandidateInversion(contig, start, end, inv_cluster.members, score))

    ##################
    # Write clusters #
    ##################
    if not parameters["skip_confirm"]:
        logging.info("Write confirmed evidence clusters..")
        deletion_evidence_output = open(options.working_dir + '/evidences/del_confirmed.bed', 'w')
        insertion_evidence_output = open(options.working_dir + '/evidences/ins_confirmed.bed', 'w')
        inversion_evidence_output = open(options.working_dir + '/evidences/inv_confirmed.bed', 'w')

        for cluster in deletion_evidence_clusters:
            print(cluster.get_bed_entry(), file=deletion_evidence_output)
        for cluster in insertion_evidence_clusters:
            print(cluster.get_bed_entry(), file=insertion_evidence_output)
        for cluster in inversion_evidence_clusters:
            print(cluster.get_bed_entry(), file=inversion_evidence_output)

        deletion_evidence_output.close()
        insertion_evidence_output.close()
        inversion_evidence_output.close()

    ###################################
    # Merge translocation breakpoints #
    ###################################

    # Cluster translocations by contig and pos1
    logging.info("Cluster translocations..")
    translocation_partitions = form_partitions(completed_translocations, parameters["trans_partition_max_distance"])

    logging.info("Compile translocation dict..")
    translocation_partitions_dict = defaultdict(list)
    for partition in translocation_partitions:
        translocation_partitions_dict[partition[0].contig1].append(partition)

    logging.info("Compute translocation means and std deviations..")
    translocation_partition_means_dict = {}
    translocation_partition_stds_dict = {}
    for contig in translocation_partitions_dict.keys():
        translocation_partition_means_dict[contig] = [int(round(sum([ev.pos1 for ev in partition]) / float(len(partition)))) for partition in translocation_partitions_dict[contig]]
        translocation_partition_stds_dict[contig] = [int(round(sqrt(sum([pow(abs(ev.pos1 - translocation_partition_means_dict[contig][index]), 2) for ev in partition]) / float(len(partition))))) for index, partition in enumerate(translocation_partitions_dict[contig])]

    insertion_candidates = []
    int_duplication_candidates = []

    logging.info("Merge translocations at deletions..")
    new_insertion_candidates = merge_translocations_at_deletions(translocation_partitions_dict, translocation_partition_means_dict, translocation_partition_stds_dict, deletion_evidence_clusters, parameters)
    insertion_candidates.extend(new_insertion_candidates)

    logging.info("Merge translocations at insertions..")
    insertion_from_evidence_clusters.extend(merge_translocations_at_insertions(translocation_partitions_dict, translocation_partition_means_dict, translocation_partition_stds_dict, insertion_evidence_clusters, deletion_evidence_clusters, parameters))
    # insertion_candidates.extend(new_insertion_candidates)
    # int_duplication_candidates.extend(new_int_duplication_candidates)

    # Merge insertions with source
    logging.info("Classify insertion/duplication evidence clusters..")
    new_insertion_candidates, new_int_duplication_candidates = merge_insertions_from(insertion_from_evidence_clusters, deletion_evidence_clusters, parameters)
    insertion_candidates.extend(new_insertion_candidates)
    int_duplication_candidates.extend(new_int_duplication_candidates)

    # Cluster candidates
    logging.info("Cluster SV candidates..")
    final_insertion_candidates, final_int_duplication_candidates = cluster_sv_candidates(insertion_candidates, int_duplication_candidates, parameters)

    #Write candidates
    logging.info("Write SV candidates..")
    write_candidates(options.working_dir, (final_insertion_candidates, final_int_duplication_candidates, inversion_candidates))
    write_final_vcf(options.working_dir, options.genome, final_insertion_candidates, final_int_duplication_candidates, inversion_candidates, deletion_evidence_clusters, tandem_duplication_evidence_clusters)

if __name__ == "__main__":
    sys.exit(main())
