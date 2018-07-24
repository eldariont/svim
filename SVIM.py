__version__ = '0.3'
__author__ = 'David Heller'

import sys
import argparse
import os
import re
import pickle
import gzip
import logging
import configparser

from time import strftime, localtime

from SVIM_COLLECT import guess_file_type, read_file_list, create_full_file, run_full_alignment, analyze_alignment
from SVIM_CLUSTER import cluster_sv_signatures, write_signature_clusters_bed, write_signature_clusters_vcf, plot_histograms
from SVIM_COMBINE import combine_clusters


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""SVIM (pronounced SWIM) is a structural variant caller for long reads. 
It discriminates six different variant classes: deletions, cut&paste insertions, tandem and interspersed duplications, 
inversions and novel element insertions. SVIM is unique in its capability of extracting both the genomic origin and 
destination of insertions and duplications.

SVIM consists of three major steps:
- COLLECT detects signatures for SVs in long read alignments
- CLUSTER merges signatures that come from the same SV
- COMBINE combines clusters from different genomic regions and classifies them into distinct SV types

SVIM can process two types of input. Firstly, it can detect SVs from raw reads by aligning them to a given reference genome first ("SVIM.py reads [options] working_dir reads genome").
Alternatively, it can detect SVs from existing reads alignments in SAM/BAM format ("SVIM.py alignment [options] working_dir bam_file").
""")
    subparsers = parser.add_subparsers(help='modes', dest='sub')
    parser.add_argument('--version', '-v', action='version', version='%(prog)s {version}'.format(version=__version__))

    parser_fasta = subparsers.add_parser('reads', help='Detect SVs from raw reads')
    parser_fasta.add_argument('working_dir', type=str, help='working directory')
    parser_fasta.add_argument('reads', type=str, help='Read file (FASTA, FASTQ, gzipped FASTA and FASTQ)')
    parser_fasta.add_argument('genome', type=str, help='Reference genome file (FASTA)')
    group_fasta_collect = parser_fasta.add_argument_group('COLLECT')
    group_fasta_collect.add_argument('--min_mapq', type=int, default=20, help='Minimum mapping quality of reads to consider')
    group_fasta_collect.add_argument('--min_sv_size', type=int, default=40, help='Minimum SV size to detect')
    group_fasta_collect.add_argument('--max_sv_size', type=int, default=100000, help='Maximum SV size to detect')
    group_fasta_collect.add_argument('--skip_indel', action='store_true', help='disable indel part')
    group_fasta_collect.add_argument('--skip_segment', action='store_true', help='disable segment part')
    group_fasta_collect.add_argument('--cores', type=int, default=1, help='CPU cores to use for alignment')
    group_fasta_collect.add_argument('--segment_gap_tolerance', type=int, default=10, help='Maximum tolerated gap between adjacent alignment segments')
    group_fasta_collect.add_argument('--segment_overlap_tolerance', type=int, default=5, help='Maximum tolerated overlap between adjacent alignment segments')
    group_fasta_cluster = parser_fasta.add_argument_group('CLUSTER')
    group_fasta_cluster.add_argument('--partition_max_distance', type=int, default=5000, help='Maximum distance in bp between SVs in a partition')
    group_fasta_cluster.add_argument('--distance_normalizer', type=int, default=900, help='Distance normalizer used for span-position distance')
    group_fasta_cluster.add_argument('--cluster_max_distance', type=float, default=0.7, help='Maximum span-position distance between SVs in a cluster')
    group_fasta_combine = parser_fasta.add_argument_group('COMBINE')
    group_fasta_combine.add_argument('--del_ins_dup_max_distance', type=float, default=1.0, help='')
    group_fasta_combine.add_argument('--trans_destination_partition_max_distance', type=int, default=1000, help='')
    group_fasta_combine.add_argument('--trans_partition_max_distance', type=int, default=200, help='')
    group_fasta_combine.add_argument('--trans_sv_max_distance', type=int, default=500, help='')

    parser_bam = subparsers.add_parser('alignment', help='Detect SVs from an existing alignment')
    parser_bam.add_argument('working_dir', type=os.path.abspath, help='working directory')
    parser_bam.add_argument('bam_file', type=argparse.FileType('r'), help='SAM/BAM file with aligned long reads (must be query-sorted)')
    group_bam_collect = parser_bam.add_argument_group('COLLECT')
    group_bam_collect.add_argument('--min_mapq', type=int, default=20, help='Minimum mapping quality of reads to consider')
    group_bam_collect.add_argument('--min_sv_size', type=int, default=40, help='Minimum SV size to detect')
    group_bam_collect.add_argument('--max_sv_size', type=int, default=100000, help='Maximum SV size to detect')
    group_bam_collect.add_argument('--skip_indel', action='store_true', help='disable indel part')
    group_bam_collect.add_argument('--skip_segment', action='store_true', help='disable segment part')
    group_bam_collect.add_argument('--segment_gap_tolerance', type=int, default=10, help='Maximum tolerated gap between adjacent alignment segments')
    group_bam_collect.add_argument('--segment_overlap_tolerance', type=int, default=5, help='Maximum tolerated overlap between adjacent alignment segments')
    group_bam_cluster = parser_bam.add_argument_group('CLUSTER')
    group_bam_cluster.add_argument('--partition_max_distance', type=int, default=5000, help='Maximum distance in bp between SVs in a partition')
    group_bam_cluster.add_argument('--distance_normalizer', type=int, default=900, help='Distance normalizer used for span-position distance')
    group_bam_cluster.add_argument('--cluster_max_distance', type=float, default=0.7, help='Maximum span-position distance between SVs in a cluster')
    group_bam_combine = parser_bam.add_argument_group('COMBINE')
    group_bam_combine.add_argument('--del_ins_dup_max_distance', type=float, default=1.0, help='')
    group_bam_combine.add_argument('--trans_destination_partition_max_distance', type=int, default=1000, help='')
    group_bam_combine.add_argument('--trans_partition_max_distance', type=int, default=200, help='')
    group_bam_combine.add_argument('--trans_sv_max_distance', type=int, default=500, help='')

    return parser.parse_args()


def main():
    # Fetch command-line options
    options = parse_arguments()
    options.distance_metric = "sl" 

    if not options.sub:
        print("Please choose one of the two modes ('reads' or 'alignment'). See --help for more information.")
        return

    # Set up logging
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-7.7s]  %(message)s")
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    # Create working dir if it does not exist
    if not os.path.exists(options.working_dir):
        os.makedirs(options.working_dir)

    # Create log file
    fileHandler = logging.FileHandler("{0}/SVIM_{1}.log".format(options.working_dir, strftime("%y%m%d_%H%M%S", localtime())), mode="w")
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    logging.info("****************** Start SVIM, version {0} ******************".format(__version__))
    logging.info("CMD: python3 {0}".format(" ".join(sys.argv)))
    logging.info("WORKING DIR: {0}".format(os.path.abspath(options.working_dir)))
    logging.info("****************** STEP 1: COLLECT ******************")

    # Search for SV signatures
    if options.sub == 'reads':
        logging.info("MODE: reads")
        logging.info("INPUT: {0}".format(os.path.abspath(options.reads)))
        logging.info("GENOME: {0}".format(os.path.abspath(options.genome)))
        reads_type = guess_file_type(options.reads)
        if reads_type == "unknown":
            return
        elif reads_type == "list":
            # List of read files
            sv_signatures = []
            for file_path in read_file_list(options.reads):
                reads_type = guess_file_type(file_path)
                full_reads_path = create_full_file(options.working_dir, file_path, reads_type)
                run_full_alignment(options.working_dir, options.genome, full_reads_path, options.cores)
                reads_file_prefix = os.path.splitext(os.path.basename(full_reads_path))[0]
                full_aln = "{0}/{1}_aln.querysorted.bam".format(options.working_dir, reads_file_prefix)
                sv_signatures.extend(analyze_alignment(full_aln, options))
        else:
            # Single read file
            full_reads_path = create_full_file(options.working_dir, options.reads, reads_type)
            run_full_alignment(options.working_dir, options.genome, full_reads_path, options.cores)
            reads_file_prefix = os.path.splitext(os.path.basename(full_reads_path))[0]
            full_aln = "{0}/{1}_aln.querysorted.bam".format(options.working_dir, reads_file_prefix)
            sv_signatures = analyze_alignment(full_aln, options)
    elif options.sub == 'alignment':
        logging.info("MODE: alignment")
        logging.info("INPUT: {0}".format(os.path.abspath(options.bam_file.name)))
        sv_signatures = analyze_alignment(options.bam_file.name, options)

    deletion_signatures = [ev for ev in sv_signatures if ev.type == 'del']
    insertion_signatures = [ev for ev in sv_signatures if ev.type == 'ins']
    inversion_signatures = [ev for ev in sv_signatures if ev.type == 'inv']
    tandem_duplication_signatures = [ev for ev in sv_signatures if ev.type == 'dup']
    translocation_signatures = [ev for ev in sv_signatures if ev.type == 'tra']
    insertion_from_signatures = [ev for ev in sv_signatures if ev.type == 'ins_dup']

    logging.info("Found {0} signatures for deleted regions.".format(len(deletion_signatures)))
    logging.info("Found {0} signatures for inserted regions.".format(len(insertion_signatures)))
    logging.info("Found {0} signatures for inverted regions.".format(len(inversion_signatures)))
    logging.info("Found {0} signatures for tandem duplicated regions.".format(len(tandem_duplication_signatures)))
    logging.info("Found {0} signatures for translocation breakpoints.".format(len(translocation_signatures)))
    logging.info("Found {0} signatures for inserted regions with detected region of origin.".format(len(insertion_from_signatures)))
    
    # Cluster SV signatures
    logging.info("****************** STEP 2: CLUSTER ******************")
    signature_clusters = cluster_sv_signatures(sv_signatures, options)

    # Write SV signature clusters
    logging.info("Finished clustering. Writing signature clusters..")
    write_signature_clusters_bed(options.working_dir, signature_clusters)
    write_signature_clusters_vcf(options.working_dir, signature_clusters, __version__)

    # Create result plots
    plot_histograms(options.working_dir, signature_clusters)

    # Dump obj file
    # signatures_file = open(options.working_dir + '/sv_signatures.obj', 'wb')
    # logging.info("Storing collected signature clusters into sv_signatures.obj..")
    # pickle.dump(signature_clusters, signatures_file)
    # signatures_file.close()

    logging.info("****************** STEP 3: COMBINE ******************")
    combine_clusters(signature_clusters, options.working_dir, options, __version__)

if __name__ == "__main__":
    sys.exit(main())