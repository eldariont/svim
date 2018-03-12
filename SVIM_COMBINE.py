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

from collections import defaultdict
from time import strftime, localtime
from math import pow, sqrt

from SVIM_clustering import form_partitions, partition_and_cluster_candidates
from SVCandidate import CandidateInversion
from SVIM_merging import merge_insertions_from, merge_translocations_at_deletions, merge_translocations_at_insertions


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""SVIM (pronounced SWIM) is a structural variant caller for long reads. 
It combines full alignment analysis, split-read mapping, and read-tail mapping to 
distinguish five classes of structural variants. SVIM discriminates between similar 
SV classes such as interspersed duplications and cut&paste insertions and is unique 
in its capability of extracting both the genomic origin and destination of insertions 
and duplications.

SVIM consists of three programs SVIM-COLLECT, SVIM-CONFIRM, and SVIM-COMBINE. You are 
running SVIM-COMBINE which combines different SV evidences and classifies them into 
distinct SV types.""")
    parser.add_argument('--version', '-v', action='version', version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('working_dir', type=str, help='working directory')
    parser.add_argument('--config', type=str, default="{0}/default_config.cfg".format(os.path.dirname(os.path.realpath(__file__))), help='configuration file, default: {0}/default_config.cfg'.format(os.path.dirname(os.path.realpath(__file__))))
    parser.add_argument('--obj_file', '-i', type=argparse.FileType('r'), help='Path of .obj file to load (default: working_dir/sv_confirmed_evidences.obj')
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


def write_final_vcf(working_dir, insertion_candidates, int_duplication_candidates, inversion_candidates, deletion_evidence_clusters, tandem_duplication_evidence_clusters):
    vcf_output = open(working_dir + '/final_results.vcf', 'w')

    # Write header lines
    print("##fileformat=VCFv4.3", file=vcf_output)
    print("##source=SVIMV{0}".format(__version__), file=vcf_output)
    #print("##reference={0}".format(genome), file=vcf_output)
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

    fileHandler = logging.FileHandler("{0}/SVIM-COMBINE_{1}.log".format(options.working_dir, strftime("%y%m%d_%H%M%S", localtime())), mode="w")
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    logging.info("****************** Start SVIM-COMBINE, version {0} ******************".format(__version__))
    logging.info("CMD: python {0}".format(" ".join(sys.argv)))
    logging.info("WORKING DIR: {0}".format(os.path.abspath(options.working_dir)))

    if options.obj_file:
        logging.info("INPUT: {0}".format(os.path.abspath(options.obj_file.name)))
        evidences_file = options.obj_file
    else:
        logging.info("INPUT: {0}".format(os.path.abspath(options.working_dir + '/sv_evidences.obj')))
        evidences_file = open(options.working_dir + '/sv_evidences.obj', 'r')

    logging.info("Loading object file created by SVIM-COLLECT or SVIM_COMBINE.")
    evidence_clusters = pickle.load(evidences_file)
    deletion_evidence_clusters, insertion_evidence_clusters, inversion_evidence_clusters, tandem_duplication_evidence_clusters, insertion_from_evidence_clusters, completed_translocations = evidence_clusters
    evidences_file.close()

    ###############################
    # Create inversion candidates #
    ###############################
    inversion_candidates = []
    for inv_cluster in inversion_evidence_clusters:
        if inv_cluster.score > 0:
            inversion_candidates.append(CandidateInversion(inv_cluster.contig, inv_cluster.start, inv_cluster.end, inv_cluster.members, inv_cluster.score))

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
    write_final_vcf(options.working_dir, final_insertion_candidates, final_int_duplication_candidates, inversion_candidates, deletion_evidence_clusters, tandem_duplication_evidence_clusters)

if __name__ == "__main__":
    sys.exit(main())
