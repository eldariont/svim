import os
import logging

from collections import defaultdict
from math import pow, sqrt

from SVIM_clustering import form_partitions, partition_and_cluster_candidates
from SVCandidate import CandidateInversion, CandidateDuplicationTandem, CandidateDeletion, CandidateNovelInsertion
from SVIM_merging import merge_insertions_from, merge_translocations_at_deletions, merge_translocations_at_insertions


def cluster_sv_candidates(insertion_candidates, int_duplication_candidates, options):
    """Takes a list of SVCandidates and splits them up by type. The SVCandidates of each type are clustered and returned as a tuple of
    (deletion_signature_clusters, insertion_signature_clusters, inversion_signature_clusters, tandem_duplication_signature_clusters, insertion_from_signature_clusters, completed_translocation_signatures)."""

    final_insertion_candidates = partition_and_cluster_candidates(insertion_candidates, options, "insertion candidates")
    final_int_duplication_candidates = partition_and_cluster_candidates(int_duplication_candidates, options, "interspersed duplication candidates")

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


def write_final_vcf(working_dir, insertion_candidates, int_duplication_candidates, inversion_candidates, tandem_duplication_candidates, deletion_candidates, novel_insertion_candidates, version):
    vcf_output = open(working_dir + '/final_results.vcf', 'w')

    # Write header lines
    print("##fileformat=VCFv4.3", file=vcf_output)
    print("##source=SVIMV{0}".format(version), file=vcf_output)
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


def combine_clusters(signature_clusters, working_dir, options, version):
    deletion_signature_clusters, insertion_signature_clusters, inversion_signature_clusters, tandem_duplication_signature_clusters, insertion_from_signature_clusters, completed_translocations = signature_clusters

    ###############################
    # Create inversion candidates #
    ###############################
    inversion_candidates = []
    for inv_cluster in inversion_signature_clusters:
        inversion_candidates.append(CandidateInversion(inv_cluster.contig, inv_cluster.start, inv_cluster.end, inv_cluster.members, inv_cluster.score, inv_cluster.std_span, inv_cluster.std_pos))

    ########################################
    # Create tandem duplication candidates #
    ########################################
    tan_dup_candidates = []
    for tan_dup_cluster in tandem_duplication_signature_clusters:
        source_contig, source_start, source_end = tan_dup_cluster.get_source()
        dest_contig, dest_start, dest_end = tan_dup_cluster.get_destination()
        num_copies = int(round((dest_end - dest_start) / (source_end - source_start)))
        tan_dup_candidates.append(CandidateDuplicationTandem(tan_dup_cluster.source_contig, tan_dup_cluster.source_start, tan_dup_cluster.source_end, num_copies, tan_dup_cluster.members, tan_dup_cluster.score, tan_dup_cluster.std_span, tan_dup_cluster.std_pos))

    ###################################
    # Merge translocation breakpoints #
    ###################################

    # Cluster translocations by contig and pos1
    logging.info("Cluster translocation breakpoints..")
    translocations_fwdfwd = [tra for tra in completed_translocations if tra.direction1 == "fwd" and tra.direction2 == "fwd"]
    translocations_revrev = [tra for tra in completed_translocations if tra.direction1 == "rev" and tra.direction2 == "rev"]
    translocation_partitions_fwdfwd = form_partitions(translocations_fwdfwd, options.trans_partition_max_distance)
    translocation_partitions_revrev = form_partitions(translocations_revrev, options.trans_partition_max_distance)

    translocation_partitions_fwdfwd_dict = defaultdict(list)
    translocation_partitions_revrev_dict = defaultdict(list)
    for partition in translocation_partitions_fwdfwd:
        translocation_partitions_fwdfwd_dict[partition[0].contig1].append(partition)
    for partition in translocation_partitions_revrev:
        translocation_partitions_revrev_dict[partition[0].contig1].append(partition)

    translocation_partition_means_fwdfwd_dict = {}
    translocation_partition_stds_fwdfwd_dict = {}
    for contig in translocation_partitions_fwdfwd_dict.keys():
        translocation_partition_means_fwdfwd_dict[contig] = [int(round(sum([ev.pos1 for ev in partition]) / len(partition))) for partition in translocation_partitions_fwdfwd_dict[contig]]
        translocation_partition_stds_fwdfwd_dict[contig] = [int(round(sqrt(sum([pow(abs(ev.pos1 - translocation_partition_means_fwdfwd_dict[contig][index]), 2) for ev in partition]) / len(partition)))) for index, partition in enumerate(translocation_partitions_fwdfwd_dict[contig])]
    translocation_partition_means_revrev_dict = {}
    translocation_partition_stds_revrev_dict = {}
    for contig in translocation_partitions_revrev_dict.keys():
        translocation_partition_means_revrev_dict[contig] = [int(round(sum([ev.pos1 for ev in partition]) / len(partition))) for partition in translocation_partitions_revrev_dict[contig]]
        translocation_partition_stds_revrev_dict[contig] = [int(round(sqrt(sum([pow(abs(ev.pos1 - translocation_partition_means_revrev_dict[contig][index]), 2) for ev in partition]) / len(partition)))) for index, partition in enumerate(translocation_partitions_revrev_dict[contig])]


    insertion_candidates = []
    int_duplication_candidates = []

    logging.info("Combine deleted regions with translocations breakpoints..")
    new_insertion_candidates, deleted_regions_to_remove_1 = merge_translocations_at_deletions(translocation_partitions_fwdfwd_dict, translocation_partition_means_fwdfwd_dict, translocation_partition_stds_fwdfwd_dict, translocation_partitions_revrev_dict, translocation_partition_means_revrev_dict, translocation_partition_stds_revrev_dict, deletion_signature_clusters, options)
    insertion_candidates.extend(new_insertion_candidates)

    logging.info("Combine deleted regions with translocation breakpoints..")
    new_insertion_from_clusters, inserted_regions_to_remove_1 = merge_translocations_at_insertions(translocation_partitions_fwdfwd_dict, translocation_partition_means_fwdfwd_dict, translocation_partition_stds_fwdfwd_dict, translocation_partitions_revrev_dict, translocation_partition_means_revrev_dict, translocation_partition_stds_revrev_dict, insertion_signature_clusters, options)
    insertion_from_signature_clusters.extend(new_insertion_from_clusters)

    ###################################
    # Classify insertions with source #
    ###################################

    logging.info("Classify inserted regions with detected region of origin..")
    new_insertion_candidates, new_int_duplication_candidates, deleted_regions_to_remove_2 = merge_insertions_from(insertion_from_signature_clusters, deletion_signature_clusters, options)
    insertion_candidates.extend(new_insertion_candidates)
    int_duplication_candidates.extend(new_int_duplication_candidates)

    ##################################
    # Remove deleted region clusters #
    ##################################
    all_deleted_regions_to_remove = sorted(list(set(deleted_regions_to_remove_1 + deleted_regions_to_remove_2)), reverse=True)
    for del_index in all_deleted_regions_to_remove:
        del(deletion_signature_clusters[del_index])

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

    for inserted_region_index, inserted_region in enumerate(insertion_signature_clusters):
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
        del(insertion_signature_clusters[ins_index])

    ##############################
    # Create deletion candidates #
    ##############################
    deletion_candidates = []
    for del_cluster in deletion_signature_clusters:
        if del_cluster.score > 0:
            deletion_candidates.append(CandidateDeletion(del_cluster.contig, del_cluster.start, del_cluster.end, del_cluster.members, del_cluster.score, del_cluster.std_span, del_cluster.std_pos))

    #####################################
    # Create novel insertion candidates #
    #####################################
    novel_insertion_candidates = []
    for ins_cluster in insertion_signature_clusters:
        if ins_cluster.score > 0:
            novel_insertion_candidates.append(CandidateNovelInsertion(ins_cluster.contig, ins_cluster.start, ins_cluster.end, ins_cluster.members, ins_cluster.score, ins_cluster.std_span, ins_cluster.std_pos))

    ######################
    # Cluster candidates #
    ######################
    logging.info("Cluster SV candidates one more time..")
    final_insertion_candidates, final_int_duplication_candidates = cluster_sv_candidates(insertion_candidates, int_duplication_candidates, options)

    ####################
    # Write candidates #
    ####################
    logging.info("Write SV candidates..")
    logging.info("Final deletion candidates: {0}".format(len(deletion_candidates)))
    logging.info("Final inversion candidates: {0}".format(len(inversion_candidates)))
    logging.info("Final insertion candidates: {0}".format(len(final_insertion_candidates)))
    logging.info("Final interspersed duplication candidates: {0}".format(len(final_int_duplication_candidates)))
    logging.info("Final tandem duplication candidates: {0}".format(len(tan_dup_candidates)))
    logging.info("Final novel insertion candidates: {0}".format(len(novel_insertion_candidates)))
    write_candidates(working_dir, (final_insertion_candidates, final_int_duplication_candidates, inversion_candidates, tan_dup_candidates, deletion_candidates, novel_insertion_candidates))
    write_final_vcf(working_dir, final_insertion_candidates, final_int_duplication_candidates, inversion_candidates, tan_dup_candidates, deletion_candidates, novel_insertion_candidates, version)
