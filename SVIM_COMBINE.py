from __future__ import print_function

__version__ = '0.2'
__author__ = 'David Heller'

import sys
import argparse
import os
import pickle
import gzip
import logging
import configparser

from collections import defaultdict
from time import strftime, localtime
from math import pow, sqrt

from SVIM_clustering import form_partitions, partition_and_cluster_candidates
from SVCandidate import CandidateInversion, CandidateDuplicationTandem, CandidateDeletion, CandidateNovelInsertion
from SVIM_merging import merge_insertions_from, merge_translocations_at_deletions, merge_translocations_at_insertions
from SVIM_COLLECT import read_parameters

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
    parser.add_argument('--obj_file', '-i', type=argparse.FileType('rb'), help='Path of .obj file to load (default: working_dir/sv_evidences.obj')
    return parser.parse_args()


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
    insertion_candidates, int_duplication_candidates, inversion_candidates, tan_duplication_candidates, deletion_candidates, novel_insertion_candidates = candidates

    if not os.path.exists(working_dir + '/candidates'):
        os.mkdir(working_dir + '/candidates')
    deletion_candidate_output = open(working_dir + '/candidates/candidates_deletions.bed', 'w')
    insertion_candidate_source_output = open(working_dir + '/candidates/candidates_insertions_source.bed', 'w')
    insertion_candidate_dest_output = open(working_dir + '/candidates/candidates_insertions_dest.bed', 'w')
    inversion_candidate_output = open(working_dir + '/candidates/candidates_inversions.bed', 'w')
    tandem_duplication_candidate_source_output = open(working_dir + '/candidates/candidates_tan_duplications_source.bed', 'w')
    tandem_duplication_candidate_dest_output = open(working_dir + '/candidates/candidates_tan_duplications_dest.bed', 'w')
    interspersed_duplication_candidate_source_output = open(working_dir + '/candidates/candidates_int_duplications_source.bed', 'w')
    interspersed_duplication_candidate_dest_output = open(working_dir + '/candidates/candidates_int_duplications_dest.bed', 'w')
    novel_insertion_candidate_output = open(working_dir + '/candidates/candidates_novel_insertions.bed', 'w')

    for candidate in deletion_candidates:
        print(candidate.get_bed_entry(), file=deletion_candidate_output)
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
    for candidate in tan_duplication_candidates:
        bed_entries = candidate.get_bed_entries()
        print(bed_entries[0], file=tandem_duplication_candidate_source_output)
        print(bed_entries[1], file=tandem_duplication_candidate_dest_output)
    for candidate in novel_insertion_candidates:
        print(candidate.get_bed_entry(), file=novel_insertion_candidate_output)

    deletion_candidate_output.close()
    insertion_candidate_source_output.close()
    insertion_candidate_dest_output.close()
    inversion_candidate_output.close()
    interspersed_duplication_candidate_source_output.close()
    interspersed_duplication_candidate_dest_output.close()
    tandem_duplication_candidate_source_output.close()
    tandem_duplication_candidate_dest_output.close()


def write_final_vcf(working_dir, insertion_candidates, int_duplication_candidates, inversion_candidates, tandem_duplication_candidates, deletion_candidates, novel_insertion_candidates):
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
    print("##ALT=<ID=INS:NOVEL,Description=\"Novel Insertion\">", file=vcf_output)
    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">", file=vcf_output)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">", file=vcf_output)
    print("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">", file=vcf_output)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=vcf_output)

    vcf_entries = []
    for candidate in deletion_candidates:
        vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry()))
    for candidate in inversion_candidates:
        vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry()))
    for candidate in tandem_duplication_candidates:
        vcf_entries.append((candidate.get_destination(), candidate.get_vcf_entry()))
    for candidate in insertion_candidates:
        vcf_entries.append((candidate.get_destination(), candidate.get_vcf_entry()))
    for candidate in int_duplication_candidates:
        vcf_entries.append((candidate.get_destination(), candidate.get_vcf_entry()))
    for candidate in novel_insertion_candidates:
        vcf_entries.append((candidate.get_destination(), candidate.get_vcf_entry()))

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
        evidences_file = open(options.working_dir + '/sv_evidences.obj', 'rb')

    logging.info("Loading object file created by SVIM-COLLECT or SVIM_CONFIRM.")
    evidence_clusters = pickle.load(evidences_file)
    deletion_evidence_clusters, insertion_evidence_clusters, inversion_evidence_clusters, tandem_duplication_evidence_clusters, insertion_from_evidence_clusters, completed_translocations = evidence_clusters
    evidences_file.close()

    ###############################
    # Create inversion candidates #
    ###############################
    inversion_candidates = []
    for inv_cluster in inversion_evidence_clusters:
        inversion_candidates.append(CandidateInversion(inv_cluster.contig, inv_cluster.start, inv_cluster.end, inv_cluster.members, inv_cluster.score, inv_cluster.std_span, inv_cluster.std_pos))

    ########################################
    # Create tandem duplication candidates #
    ########################################
    tan_dup_candidates = []
    for tan_dup_cluster in tandem_duplication_evidence_clusters:
        source_contig, source_start, source_end = tan_dup_cluster.get_source()
        dest_contig, dest_start, dest_end = tan_dup_cluster.get_destination()
        num_copies = int(round((dest_end - dest_start) / (source_end - source_start)))
        tan_dup_candidates.append(CandidateDuplicationTandem(tan_dup_cluster.source_contig, tan_dup_cluster.source_start, tan_dup_cluster.source_end, num_copies, tan_dup_cluster.members, tan_dup_cluster.score, tan_dup_cluster.std_span, tan_dup_cluster.std_pos))

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
        translocation_partition_means_dict[contig] = [int(round(sum([ev.pos1 for ev in partition]) / len(partition))) for partition in translocation_partitions_dict[contig]]
        translocation_partition_stds_dict[contig] = [int(round(sqrt(sum([pow(abs(ev.pos1 - translocation_partition_means_dict[contig][index]), 2) for ev in partition]) / len(partition)))) for index, partition in enumerate(translocation_partitions_dict[contig])]

    insertion_candidates = []
    int_duplication_candidates = []

    logging.info("Merge translocations at deletions..")
    new_insertion_candidates, deleted_regions_to_remove_1 = merge_translocations_at_deletions(translocation_partitions_dict, translocation_partition_means_dict, translocation_partition_stds_dict, deletion_evidence_clusters, parameters)
    insertion_candidates.extend(new_insertion_candidates)

    logging.info("Merge translocations at insertions..")
    new_insertion_from_clusters, inserted_regions_to_remove_1 = merge_translocations_at_insertions(translocation_partitions_dict, translocation_partition_means_dict, translocation_partition_stds_dict, insertion_evidence_clusters, parameters)
    insertion_from_evidence_clusters.extend(new_insertion_from_clusters)

    ###################################
    # Classify insertions with source #
    ###################################

    logging.info("Classify insertion/duplication evidence clusters..")
    new_insertion_candidates, new_int_duplication_candidates, deleted_regions_to_remove_2 = merge_insertions_from(insertion_from_evidence_clusters, deletion_evidence_clusters, parameters)
    insertion_candidates.extend(new_insertion_candidates)
    int_duplication_candidates.extend(new_int_duplication_candidates)

    ##################################
    # Remove deleted region clusters #
    ##################################
    all_deleted_regions_to_remove = sorted(list(set(deleted_regions_to_remove_1 + deleted_regions_to_remove_2)), reverse=True)
    for del_index in all_deleted_regions_to_remove:
        del(deletion_evidence_clusters[del_index])

    ###################################
    # Remove inserted region clusters #
    ###################################

    #find all inserted regions overlapping insertion, interspersed duplication or tandem duplication candidates
    insertion_iterator = iter(sorted(insertion_candidates, key=lambda cand: cand.get_destination()))
    int_duplication_iterator = iter(sorted(int_duplication_candidates, key=lambda cand: cand.get_destination()))
    tan_duplication_iterator = iter(sorted(tan_dup_candidates, key=lambda cand: cand.get_destination()))
    insertions_end = False
    int_duplications_end = False
    tan_duplications_end = False
    inserted_regions_to_remove_2 = []

    try:
        current_insertion = next(insertion_iterator)
    except StopIteration:
        insertions_end = True

    try:
        current_int_duplication = next(int_duplication_iterator)
    except StopIteration:
        int_duplications_end = True

    try:
        current_tan_duplication = next(tan_duplication_iterator)
    except StopIteration:
        tan_duplications_end = True

    for inserted_region_index, inserted_region in enumerate(insertion_evidence_clusters):
        contig1, start1, end1 = inserted_region.get_source()
        if not insertions_end:
            contig2, start2, end2 = current_insertion.get_destination()
            while contig2 < contig1 or (contig2 == contig1 and end2 < start1):
                try:
                    current_insertion = next(insertion_iterator)
                    contig2, start2, end2 = current_insertion.get_destination()
                except StopIteration:
                    insertions_end = True
                    break

        #if overlapping insertion
        if not insertions_end and contig2 == contig1 and start2 < end1:
            inserted_regions_to_remove_2.append(inserted_region_index)
        else:
            if not int_duplications_end:
                contig2, start2, end2 = current_int_duplication.get_destination()
                while contig2 < contig1 or (contig2 == contig1 and end2 < start1):
                    try:
                        current_int_duplication = next(int_duplication_iterator)
                        contig2, start2, end2 = current_int_duplication.get_destination()
                    except StopIteration:
                        int_duplications_end = True
                        break
            #if overlapping interspersed duplication
            if not int_duplications_end and contig2 == contig1 and start2 < end1:
                inserted_regions_to_remove_2.append(inserted_region_index)
            else:
                if not tan_duplications_end:
                    contig2, start2, end2 = current_tan_duplication.get_destination()
                    while contig2 < contig1 or (contig2 == contig1 and end2 < start1):
                        try:
                            current_tan_duplication = next(tan_duplication_iterator)
                            contig2, start2, end2 = current_tan_duplication.get_destination()
                        except StopIteration:
                            tan_duplications_end = True
                            break
                #if overlapping tandem duplication
                if not tan_duplications_end and contig2 == contig1 and start2 < end1:
                    inserted_regions_to_remove_2.append(inserted_region_index)

    # remove found inserted regions
    all_inserted_regions_to_remove = sorted(list(set(inserted_regions_to_remove_1 + inserted_regions_to_remove_2)), reverse=True)
    for ins_index in all_inserted_regions_to_remove:
        del(insertion_evidence_clusters[ins_index])

    ##############################
    # Create deletion candidates #
    ##############################
    deletion_candidates = []
    for del_cluster in deletion_evidence_clusters:
        if del_cluster.score > 0:
            deletion_candidates.append(CandidateDeletion(del_cluster.contig, del_cluster.start, del_cluster.end, del_cluster.members, del_cluster.score, del_cluster.std_span, del_cluster.std_pos))

    #####################################
    # Create novel insertion candidates #
    #####################################
    novel_insertion_candidates = []
    for ins_cluster in insertion_evidence_clusters:
        if ins_cluster.score > 0:
            novel_insertion_candidates.append(CandidateNovelInsertion(ins_cluster.contig, ins_cluster.start, ins_cluster.end, ins_cluster.members, ins_cluster.score, ins_cluster.std_span, ins_cluster.std_pos))

    ######################
    # Cluster candidates #
    ######################
    logging.info("Cluster SV candidates..")
    final_insertion_candidates, final_int_duplication_candidates = cluster_sv_candidates(insertion_candidates, int_duplication_candidates, parameters)

    ####################
    # Write candidates #
    ####################
    logging.info("Write SV candidates..")
    write_candidates(options.working_dir, (final_insertion_candidates, final_int_duplication_candidates, inversion_candidates, tan_dup_candidates, deletion_candidates, novel_insertion_candidates))
    write_final_vcf(options.working_dir, final_insertion_candidates, final_int_duplication_candidates, inversion_candidates, tan_dup_candidates, deletion_candidates, novel_insertion_candidates)

if __name__ == "__main__":
    sys.exit(main())
