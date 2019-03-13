import sys
import os
import logging
import argparse


def parse_arguments(program_version, arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""SVIM (pronounced SWIM) is a structural variant caller for long reads. 
It discriminates five different variant classes: deletions, tandem and interspersed duplications, 
inversions and novel element insertions. SVIM is unique in its capability of extracting both the genomic origin and 
destination of duplications.

SVIM consists of three major steps:
- COLLECT detects signatures for SVs in long read alignments
- CLUSTER merges signatures that come from the same SV
- COMBINE combines clusters from different genomic regions and classifies them into distinct SV types

SVIM can process two types of input. Firstly, it can detect SVs from raw reads by aligning them to a given reference genome first ("SVIM.py reads [options] working_dir reads genome").
Alternatively, it can detect SVs from existing reads alignments in SAM/BAM format ("SVIM.py alignment [options] working_dir bam_file").
""")
    subparsers = parser.add_subparsers(help='modes', dest='sub')
    parser.add_argument('--version', '-v', action='version', version='%(prog)s {version}'.format(version=program_version))

    parser_fasta = subparsers.add_parser('reads', help='Detect SVs from raw reads. Align reads to given reference genome first.')
    parser_fasta.add_argument('working_dir', type=str, help='working directory')
    parser_fasta.add_argument('reads', type=str, help='Read file (FASTA, FASTQ, gzipped FASTA and FASTQ)')
    parser_fasta.add_argument('genome', type=str, help='Reference genome file (FASTA)')
    group_fasta_collect = parser_fasta.add_argument_group('COLLECT')
    group_fasta_collect.add_argument('--min_mapq', type=int, default=20, help='Minimum mapping quality of reads to consider (default: %(default)s)')
    group_fasta_collect.add_argument('--min_sv_size', type=int, default=40, help='Minimum SV size to detect (default: %(default)s)')
    group_fasta_collect.add_argument('--max_sv_size', type=int, default=100000, help='Maximum SV size to detect (default: %(default)s)')
    group_fasta_collect.add_argument('--skip_indel', action='store_true', help='disable signature collection from within read alignments (default: %(default)s)')
    group_fasta_collect.add_argument('--skip_segment', action='store_true', help='disable signature collection from between read alignments (default: %(default)s)')
    group_fasta_collect.add_argument('--cores', type=int, default=1, help='CPU cores to use for alignment with ngmlr (default: %(default)s)')
    group_fasta_collect.add_argument('--aligner', type=str, default="ngmlr", choices=["ngmlr", "minimap2"], help='tool for read alignment: ngmlr or minimap2 (default: %(default)s)')
    group_fasta_collect.add_argument('--nanopore', action='store_true', help='use Nanopore settings for read alignment (default: %(default)s)')
    group_fasta_collect.add_argument('--segment_gap_tolerance', type=int, default=10, help='Maximum tolerated gap between adjacent alignment segments (default: %(default)s)')
    group_fasta_collect.add_argument('--segment_overlap_tolerance', type=int, default=5, help='Maximum tolerated overlap between adjacent alignment segments (default: %(default)s)')
    group_fasta_cluster = parser_fasta.add_argument_group('CLUSTER')
    group_fasta_cluster.add_argument('--partition_max_distance', type=int, default=5000, help='Maximum distance in bp between SVs in a partition (default: %(default)s)')
    group_fasta_cluster.add_argument('--distance_normalizer', type=int, default=900, help='Distance normalizer used for span-position distance (default: %(default)s)')
    group_fasta_cluster.add_argument('--cluster_max_distance', type=float, default=0.3, help='Maximum span-position distance between SVs in a cluster (default: %(default)s)')
    group_fasta_combine = parser_fasta.add_argument_group('COMBINE')
    group_fasta_combine.add_argument('--del_ins_dup_max_distance', type=float, default=1.0, help='Maximum span-position distance between the origin of an insertion and a deletion to be flagged as a potential cut&paste insertion (default: %(default)s)')
    group_fasta_combine.add_argument('--trans_destination_partition_max_distance', type=int, default=1000, help='Maximum distance in bp between translocation breakpoint destinations in a partition (default: %(default)s)')
    group_fasta_combine.add_argument('--trans_partition_max_distance', type=int, default=200, help='Maximum distance in bp between translocation breakpoints in a partition (default: %(default)s)')
    group_fasta_combine.add_argument('--trans_sv_max_distance', type=int, default=500, help='Maximum distance in bp between a translocation breakpoint and an SV signature to be combined (default: %(default)s)')
    group_fasta_combine.add_argument('--sample', type=str, default="Sample", help='Sample ID to include in output vcf (default: %(default)s)')

    parser_bam = subparsers.add_parser('alignment', help='Detect SVs from an existing alignment')
    parser_bam.add_argument('working_dir', type=os.path.abspath, help='working directory')
    parser_bam.add_argument('bam_file', type=str, help='SAM/BAM file with aligned long reads (sorted either by coordinate or queryname)')
    group_bam_collect = parser_bam.add_argument_group('COLLECT')
    group_bam_collect.add_argument('--min_mapq', type=int, default=20, help='Minimum mapping quality of reads to consider (default: %(default)s)')
    group_bam_collect.add_argument('--min_sv_size', type=int, default=40, help='Minimum SV size to detect (default: %(default)s)')
    group_bam_collect.add_argument('--max_sv_size', type=int, default=100000, help='Maximum SV size to detect (default: %(default)s)')
    group_bam_collect.add_argument('--skip_indel', action='store_true', help='disable signature collection from within read alignments (default: %(default)s)')
    group_bam_collect.add_argument('--skip_segment', action='store_true', help='disable signature collection from between read alignments (default: %(default)s)')
    group_bam_collect.add_argument('--segment_gap_tolerance', type=int, default=10, help='Maximum tolerated gap between adjacent alignment segments (default: %(default)s)')
    group_bam_collect.add_argument('--segment_overlap_tolerance', type=int, default=5, help='Maximum tolerated overlap between adjacent alignment segments (default: %(default)s)')
    group_bam_cluster = parser_bam.add_argument_group('CLUSTER')
    group_bam_cluster.add_argument('--partition_max_distance', type=int, default=5000, help='Maximum distance in bp between SVs in a partition (default: %(default)s)')
    group_bam_cluster.add_argument('--distance_normalizer', type=int, default=900, help='Distance normalizer used for span-position distance (default: %(default)s)')
    group_bam_cluster.add_argument('--cluster_max_distance', type=float, default=0.3, help='Maximum span-position distance between SVs in a cluster (default: %(default)s)')
    group_bam_combine = parser_bam.add_argument_group('COMBINE')
    group_bam_combine.add_argument('--del_ins_dup_max_distance', type=float, default=1.0, help='Maximum span-position distance between the origin of an insertion and a deletion to be flagged as a potential cut&paste insertion (default: %(default)s)')
    group_bam_combine.add_argument('--trans_destination_partition_max_distance', type=int, default=1000, help='Maximum distance in bp between translocation breakpoint destinations in a partition (default: %(default)s)')
    group_bam_combine.add_argument('--trans_partition_max_distance', type=int, default=200, help='Maximum distance in bp between translocation breakpoints in a partition (default: %(default)s)')
    group_bam_combine.add_argument('--trans_sv_max_distance', type=int, default=500, help='Maximum distance in bp between a translocation breakpoint and an SV signature to be combined (default: %(default)s)')
    group_bam_combine.add_argument('--sample', type=str, default="Sample", help='Sample ID to include in output vcf (default: %(default)s)')

    return parser.parse_args(arguments)


def guess_file_type(reads_path):
    if reads_path.endswith(".fa") or reads_path.endswith(".fasta") or reads_path.endswith(".FA"):
        logging.info("Recognized reads file as FASTA format.")
        return "fasta"
    elif reads_path.endswith(".fq") or reads_path.endswith(".fastq") or reads_path.endswith(".FQ"):
        logging.info("Recognized reads file as FASTQ format.")
        return "fastq"
    elif reads_path.endswith(".fa.gz") or reads_path.endswith(".fasta.gz") or reads_path.endswith(".FA.gz") or reads_path.endswith(".fa.gzip") or reads_path.endswith(".fasta.gzip") or reads_path.endswith(".FA.gzip"):
        logging.info("Recognized reads file as gzipped FASTA format.")
        return "fasta_gzip"
    elif reads_path.endswith(".fq.gz") or reads_path.endswith(".fastq.gz") or reads_path.endswith(".FQ.gz") or reads_path.endswith(".fq.gzip") or reads_path.endswith(".fastq.gzip") or reads_path.endswith(".FQ.gzip"):
        logging.info("Recognized reads file as gzipped FASTQ format.")
        return "fastq_gzip"
    elif reads_path.endswith(".fa.fn") or reads_path.endswith(".fasta.fn") or reads_path.endswith(".FA.fn") or reads_path.endswith(".fq.fn") or reads_path.endswith(".fastq.fn") or reads_path.endswith(".FQ.fn"):
        logging.info("Recognized reads file as file list format.")
        return "list"
    else:
        logging.error("Unknown file ending of file {0}. See github.com/eldariont/svim/wiki/ for supported file endings. Exiting.".format(reads_path))
        return "unknown"


def read_file_list(path):
    file_list = open(path, "r")
    for line in file_list:
        yield line.strip()
    file_list.close()