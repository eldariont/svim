from __future__ import print_function

__version__ = '0.1'
__author__ = 'David Heller'

import sys
import argparse
import os
import re
import pickle
import gzip
import logging
import configparser

from subprocess import Popen, PIPE
from collections import defaultdict
from time import strftime, localtime

import pysam

from SVIM_fullread import analyze_full_read_indel
from SVIM_splitread import analyze_full_read_segments
from SVIM_postprocessing import cluster_sv_evidences, write_evidence_clusters_bed, write_evidence_clusters_vcf, plot_histograms


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""SVIM (pronounced SWIM) is a structural variant caller for long reads. 
It combines full alignment analysis, split-read mapping, and read-tail mapping to 
distinguish five classes of structural variants. SVIM discriminates between similar 
SV classes such as interspersed duplications and cut&paste insertions and is unique 
in its capability of extracting both the genomic origin and destination of insertions 
and duplications.

SVIM consists of three programs SVIM-COLLECT, SVIM-CONFIRM, and SVIM-COMBINE. You are running SVIM-COLLECT 
which analyzes read alignments to collect evidences for SVs.

SVIM-COLLECT performs three steps to detect SVs: 
1) Alignment, 
2) SV detection,
3) Clustering""")
    subparsers = parser.add_subparsers(help='modes', dest='sub')
    parser.add_argument('--version', '-v', action='version', version='%(prog)s {version}'.format(version=__version__))

    parser_fasta = subparsers.add_parser('reads', help='Detect SVs from raw reads. Perform steps 1-3.')
    parser_fasta.add_argument('working_dir', type=str, help='working directory')
    parser_fasta.add_argument('reads', type=str, help='Read file (FASTA, FASTQ, gzipped FASTA and FASTQ)')
    parser_fasta.add_argument('genome', type=str, help='Reference genome file (FASTA)')
    parser_fasta.add_argument('--config', type=str, default="{0}/default_config.cfg".format(os.path.dirname(os.path.realpath(__file__))), help='configuration file, default: {0}/default_config.cfg'.format(os.path.dirname(os.path.realpath(__file__))))
    parser_fasta.add_argument('--skip_indel', action='store_true', help='disable indel part')
    parser_fasta.add_argument('--skip_segment', action='store_true', help='disable segment part')
    parser_fasta.add_argument('--cores', type=int, default=1, help='CPU cores to use for alignment')

    parser_bam = subparsers.add_parser('alignment', help='Detect SVs from an existing alignment. Perform steps 2-3.')
    parser_bam.add_argument('working_dir', type=os.path.abspath, help='working directory')
    parser_bam.add_argument('bam_file', type=argparse.FileType('r'), help='SAM/BAM file with aligned long reads (must be query-sorted)')
    parser_bam.add_argument('--config', type=str, default="{0}/default_config.cfg".format(os.path.dirname(os.path.realpath(__file__))), help='configuration file, default: {0}/default_config.cfg'.format(os.path.dirname(os.path.realpath(__file__))))
    parser_bam.add_argument('--skip_indel', action='store_true', help='disable indel part')
    parser_bam.add_argument('--skip_segment', action='store_true', help='disable segment part')

    return parser.parse_args()


def read_parameters(options):
    config = configparser.RawConfigParser(inline_comment_prefixes=';')
    config.read(options.config)

    parameters = dict()
    parameters["min_mapq"] = config.getint("detection", "min_mapq")
    parameters["max_sv_size"] = config.getint("detection", "max_sv_size")
    parameters["min_sv_size"] = config.getint("detection", "min_sv_size")

    parameters["max_segment_gap_tolerance"] = config.getint("split read", "max_segment_gap_tolerance")
    parameters["max_deletion_size"] = config.getint("split read", "max_deletion_size")
    parameters["segment_overlap_tolerance"] = config.getint("split read", "segment_overlap_tolerance")

    parameters["distance_metric"] = config.get("clustering", "distance_metric")
    parameters["distance_normalizer"] = config.getint("clustering", "distance_normalizer")
    parameters["partition_max_distance"] = config.getint("clustering", "partition_max_distance")
    parameters["cluster_max_distance"] = config.getfloat("clustering", "cluster_max_distance")

    parameters["del_ins_dup_max_distance"] = config.getfloat("merging", "del_ins_dup_max_distance")
    parameters["trans_destination_partition_max_distance"] = config.getint("merging", "trans_destination_partition_max_distance")
    parameters["trans_partition_max_distance"] = config.getint("merging", "trans_partition_max_distance")
    parameters["trans_sv_max_distance"] = config.getint("merging", "trans_sv_max_distance")

    parameters["max_confirmation_number"] = config.getint("confirmation", "max_confirmation_number")
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


def parse_sam_file(sam, contig):
    """Parses a SAM file and returns a dict of reads (list of alignments for each read) for a given reference contig"""
    alns = sam.fetch(reference=contig)
    aln_dict = defaultdict(list)
    for aln in alns:
        aln_dict[aln.query_name].append(aln)
    return aln_dict


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


def create_full_file(working_dir, reads_path, reads_type):
    """Create FASTA file with full reads if it does not exist."""
    if not os.path.exists(working_dir):
        logging.error("Given working directory does not exist")
        sys.exit()

    reads_file_prefix = os.path.splitext(os.path.basename(reads_path))[0]

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

    if write_full:
        if reads_type == "fasta" or reads_type == "fastq":
            reads_file = open(reads_path, "r")
        elif reads_type == "fasta_gzip" or reads_type == "fastq_gzip":
            reads_file = gzip.open(reads_path, "rb")

        if reads_type == "fasta" or reads_type == "fasta_gzip":
            sequence = ""
            for line in reads_file:
                if line.startswith('>'):
                    if sequence != "":
                        if write_full:
                            print(">" + read_name, file=full_file)
                            print(sequence, file=full_file)
                    read_name = line.strip()[1:]
                    sequence = ""
                else:
                    sequence += line.strip()
            if sequence != "":
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
                    if write_full:
                        print(">" + read_name, file=full_file)
                        print(sequence, file=full_file)
            reads_file.close()

    if write_full:
        full_file.close()

    logging.info("Full read files written")
    return full_reads_path


def run_full_alignment(working_dir, genome, reads_path, cores):
    """Align full reads with NGM-LR."""
    reads_file_prefix = os.path.splitext(os.path.basename(reads_path))[0]
    full_aln = "{0}/{1}_aln.querysorted.bam".format(working_dir, reads_file_prefix)

    if not os.path.exists(full_aln):
        ngmlr = Popen(['ngmlr',
                       '-t', str(cores), '-r', genome, '-q', os.path.realpath(reads_path)], stdout=PIPE)
        view = Popen(['samtools',
                      'view', '-b', '-@', str(cores)], stdin=ngmlr.stdout, stdout=PIPE)
        sort = Popen(['samtools',
                      'sort', '-n', '-@', str(cores), '-o', full_aln],
                     stdin=view.stdout)
        sort.wait()
        logging.info("Alignment finished")
    else:
        logging.warning("Alignment for full sequences exists. Skip")


def natural_representation(qname): 
    """Splits a read name into a tuple of strings and numbers. This facilitates the sort order applied by samtools -n
       See https://www.biostars.org/p/102735/"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', qname)]


def bam_iterator(bam):
    """Returns an iterator for the given SAM/BAM file (must be query-sorted). 
    In each call, the alignments of a single read are yielded as a 4-tuple: (read_name, list of primary pysam.AlignedSegment, list of supplementary pysam.AlignedSegment, list of secondary pysam.AlignedSegment)."""
    alignments = bam.fetch(until_eof=True)
    current_aln = next(alignments)
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
            next_aln = next(alignments)
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


def analyze_alignment(bam_path, parameters):
    full_bam = pysam.AlignmentFile(bam_path)
    full_it = bam_iterator(full_bam)

    sv_evidences = []
    read_nr = 0

    while True:
        try:
            full_iterator_object = next(full_it)
            read_nr += 1
            if read_nr % 10000 == 0:
                logging.info("Processed read {0}".format(read_nr))
            if not parameters["skip_indel"]:
                sv_evidences.extend(analyze_full_read_indel(full_iterator_object, full_bam, parameters))
            if not parameters["skip_segment"]:
                sv_evidences.extend(analyze_full_read_segments(full_iterator_object, full_bam, parameters))
        except StopIteration:
            break
        except KeyboardInterrupt:
            logging.warning('Execution interrupted by user. Stop detection and continue with next step..')
            break
    return sv_evidences


def read_file_list(path):
    file_list = open(path, "r")
    for line in file_list:
        yield line.strip()
    file_list.close()


def main():
    # Fetch command-line options and configuration file values and set parameters accordingly
    options = parse_arguments()
    parameters = read_parameters(options)

    # Set up logging
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-7.7s]  %(message)s")
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    fileHandler = logging.FileHandler("{0}/SVIM-COLLECT_{1}.log".format(options.working_dir, strftime("%y%m%d_%H%M%S", localtime())), mode="w")
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    logging.info("****************** Start SVIM-COLLECT, version {0} ******************".format(__version__))
    logging.info("CMD: python {0}".format(" ".join(sys.argv)))
    logging.info("WORKING DIR: {0}".format(os.path.abspath(options.working_dir)))

    # Search for SV evidences
    if options.sub == 'reads':
        logging.info("MODE: reads")
        logging.info("INPUT: {0}".format(os.path.abspath(options.reads)))
        logging.info("GENOME: {0}".format(os.path.abspath(options.genome)))
        reads_type = guess_file_type(options.reads)
        if reads_type == "unknown":
            return
        elif reads_type == "list":
            # List of read files
            sv_evidences = []
            for file_path in read_file_list(options.reads):
                reads_type = guess_file_type(file_path)
                full_reads_path = create_full_file(options.working_dir, file_path, reads_type)
                run_full_alignment(options.working_dir, options.genome, full_reads_path, options.cores)
                reads_file_prefix = os.path.splitext(os.path.basename(full_reads_path))[0]
                full_aln = "{0}/{1}_aln.querysorted.bam".format(options.working_dir, reads_file_prefix)
                sv_evidences.extend(analyze_alignment(full_aln, parameters))
        else:
            # Single read file
            full_reads_path = create_full_file(options.working_dir, options.reads, reads_type)
            run_full_alignment(options.working_dir, options.genome, full_reads_path, options.cores)
            reads_file_prefix = os.path.splitext(os.path.basename(full_reads_path))[0]
            full_aln = "{0}/{1}_aln.querysorted.bam".format(options.working_dir, reads_file_prefix)
            sv_evidences = analyze_alignment(full_aln, parameters)
    elif options.sub == 'alignment':
        logging.info("MODE: alignment")
        logging.info("INPUT: {0}".format(os.path.abspath(options.bam_file.name)))
        sv_evidences = analyze_alignment(options.bam_file.name, parameters)

    # Cluster SV evidences
    logging.info("Cluster SV evidences..")
    evidence_clusters = cluster_sv_evidences(sv_evidences, parameters)

    # Write SV evidence clusters
    logging.info("Write evidence clusters..")
    write_evidence_clusters_bed(options.working_dir, evidence_clusters)
    write_evidence_clusters_vcf(options.working_dir, evidence_clusters)

    # Create result plots
    plot_histograms(options.working_dir, evidence_clusters)

    # Dump obj file
    evidences_file = open(options.working_dir + '/sv_evidences.obj', 'wb')
    logging.info("Storing collected evidence clusters into sv_evidences.obj..")
    pickle.dump(evidence_clusters, evidences_file)
    evidences_file.close()


if __name__ == "__main__":
    sys.exit(main())