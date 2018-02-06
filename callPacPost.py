from __future__ import print_function

__version__ = '0.1'
__author__ = 'David Heller'

import sys
import os
import logging

from collections import defaultdict
from math import pow, sqrt

from callPacCluster import partition_and_cluster_unilocal, partition_and_cluster_bilocal, partition_and_cluster_candidates, form_partitions
from SVEvidence import EvidenceTranslocation
from callPacMerge import merge_insertions_from, merge_translocations_at_deletions, merge_translocations_at_insertions, filter_inversions


def complete_translocations(translocation_evidences):
    """Generate a complete list of translocation by adding all reversed translocations"""

    reversed_translocations = []
    for evidence in translocation_evidences:
        reversed_translocations.append(EvidenceTranslocation(evidence.contig2, evidence.pos2, evidence.contig1, evidence.pos1, evidence.evidence, evidence.read))
    return translocation_evidences + reversed_translocations


def cluster_sv_evidences(sv_evidences, parameters):
    """Takes a list of SVEvidences and splits them up by type. The SVEvidences of each type are clustered and returned as a tuple of
    (deletion_evidence_clusters, insertion_evidence_clusters, inversion_evidence_clusters, tandem_duplication_evidence_clusters, insertion_from_evidence_clusters, completed_translocation_evidences)."""

    deletion_evidences = [ev for ev in sv_evidences if ev.type == 'del']
    insertion_evidences = [ev for ev in sv_evidences if ev.type == 'ins']
    inversion_evidences = [ev for ev in sv_evidences if ev.type == 'inv']
    tandem_duplication_evidences = [ev for ev in sv_evidences if ev.type == 'dup']
    translocation_evidences = [ev for ev in sv_evidences if ev.type == 'tra']
    insertion_from_evidences = [ev for ev in sv_evidences if ev.type == 'ins_dup']

    logging.info("Found {0}/{1}/{2}/{3}/{4}/{5} evidences for deletions, insertions, inversions, tandem duplications, translocations, and insertion_from, respectively.".format(
        len(deletion_evidences), len(insertion_evidences), len(inversion_evidences), len(tandem_duplication_evidences), len(translocation_evidences), len(insertion_from_evidences)))
    
    # Cluster SV evidences
    logging.info("Cluster deletion evidences:")
    deletion_evidence_clusters = partition_and_cluster_unilocal(deletion_evidences, parameters)
    logging.info("Cluster insertion evidences:")
    insertion_evidence_clusters = partition_and_cluster_unilocal(insertion_evidences, parameters)
    logging.info("Cluster inversion evidences:")
    inversion_evidence_clusters = partition_and_cluster_unilocal(inversion_evidences, parameters)
    logging.info("Cluster tandem duplication evidences:")
    tandem_duplication_evidence_clusters = partition_and_cluster_bilocal(tandem_duplication_evidences, parameters)
    logging.info("Cluster insertion evidences with source:")
    insertion_from_evidence_clusters = partition_and_cluster_bilocal(insertion_from_evidences, parameters)

    return (deletion_evidence_clusters, insertion_evidence_clusters, inversion_evidence_clusters, tandem_duplication_evidence_clusters, insertion_from_evidence_clusters, complete_translocations(translocation_evidences))


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


def write_evidence_clusters_bed(working_dir, clusters):
    """Write evidence clusters into working directory in BED format."""
    deletion_evidence_clusters, insertion_evidence_clusters, inversion_evidence_clusters, tandem_duplication_evidence_clusters, insertion_from_evidence_clusters, completed_translocations = clusters

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
    for translocation in completed_translocations:
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(translocation.contig1, translocation.pos1, translocation.pos1+1, ">{0}:{1}".format(translocation.contig2, translocation.pos2), translocation.evidence, translocation.read), file=translocation_evidence_output)


def write_evidence_clusters_vcf(working_dir, clusters, genome):
    """Write evidence clusters into working directory in BED format."""
    deletion_evidence_clusters, insertion_evidence_clusters, inversion_evidence_clusters, tandem_duplication_evidence_clusters, insertion_from_evidence_clusters, completed_translocations = clusters

    if not os.path.exists(working_dir + '/evidences'):
        os.mkdir(working_dir + '/evidences')
    vcf_output = open(working_dir + '/evidences/all.vcf', 'w')

    # Write header lines
    print("##fileformat=VCFv4.3", file=vcf_output)
    print("##source=callPacV{0}".format(__version__), file=vcf_output)
    print("##reference={0}".format(genome), file=vcf_output)
    print("##ALT=<ID=DEL,Description=\"Deletion\">", file=vcf_output)
    print("##ALT=<ID=INV,Description=\"Inversion\">", file=vcf_output)
    print("##ALT=<ID=DUP,Description=\"Duplication\">", file=vcf_output)
    print("##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">", file=vcf_output)
    print("##ALT=<ID=INS,Description=\"Insertion\">", file=vcf_output)
    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">", file=vcf_output)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">", file=vcf_output)
    print("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">", file=vcf_output)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=vcf_output)

    vcf_entries = []
    for cluster in deletion_evidence_clusters:
        vcf_entries.append((cluster.get_source(), cluster.get_vcf_entry()))
    for cluster in insertion_evidence_clusters:
        vcf_entries.append((cluster.get_source(), cluster.get_vcf_entry()))
    for cluster in inversion_evidence_clusters:
        vcf_entries.append((cluster.get_source(), cluster.get_vcf_entry()))
    for cluster in tandem_duplication_evidence_clusters:
        vcf_entries.append((cluster.get_source(), cluster.get_vcf_entry()))

    # Sort and write entries to VCF
    for source, entry in sorted(vcf_entries, key=lambda pair: pair[0]):
        print(entry, file=vcf_output)

    vcf_output.close()

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


def write_final_vcf(working_dir, genome, insertion_candidates, int_duplication_candidates, inversion_candidates, deletion_evidence_clusters, tandem_duplication_evidence_clusters):
    vcf_output = open(working_dir + '/final_results.vcf', 'w')

    # Write header lines
    print("##fileformat=VCFv4.3", file=vcf_output)
    print("##source=callPacV{0}".format(__version__), file=vcf_output)
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

def post_processing(sv_evidences, working_dir, genome, parameters):
    #####################
    # Cluster evidences #
    #####################
    logging.info("Cluster SV evidences..")
    evidence_clusters = cluster_sv_evidences(sv_evidences, parameters)

    ##################
    # Write clusters #
    ##################
    logging.info("Write evidence clusters..")
    write_evidence_clusters_bed(working_dir, evidence_clusters)
    write_evidence_clusters_vcf(working_dir, evidence_clusters, genome)

    deletion_evidence_clusters, insertion_evidence_clusters, inversion_evidence_clusters, tandem_duplication_evidence_clusters, insertion_from_evidence_clusters, completed_translocations = evidence_clusters

    insertion_candidates = []
    int_duplication_candidates = []

    ###################################
    # Merge translocation breakpoints #
    ###################################

    # Cluster translocations by contig and pos1
    logging.info("Cluster translocations..")
    translocation_partitions = form_partitions(completed_translocations, parameters.trans_partition_max_distance)

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

    logging.info("Filter inversions..")
    inversion_candidates = filter_inversions(inversion_evidence_clusters, parameters)

    # Cluster candidates
    logging.info("Cluster SV candidates..")
    final_insertion_candidates, final_int_duplication_candidates = cluster_sv_candidates(insertion_candidates, int_duplication_candidates, parameters)

    #Write candidates
    logging.info("Write SV candidates..")
    write_candidates(working_dir, (final_insertion_candidates, final_int_duplication_candidates, inversion_candidates))
    write_final_vcf(working_dir, genome, final_insertion_candidates, final_int_duplication_candidates, inversion_candidates, deletion_evidence_clusters, tandem_duplication_evidence_clusters)