from __future__ import print_function

import sys
from bisect import bisect_left
from collections import defaultdict
from math import pow, sqrt

from SVEvidence import EvidenceTranslocation, EvidenceInsertionFrom, EvidenceClusterBiLocal
from SVCandidate import CandidateInsertion, CandidateDuplicationInterspersed
from callPacCluster import form_partitions

def merge_insertions_from(insertion_from_evidence_clusters, deletion_evidence_clusters):
    insertion_candidates = []
    int_duplication_candidates = []  
    for ins_cluster in insertion_from_evidence_clusters:    
        distances = [(ind, del_cluster.gowda_diday_distance(ins_cluster, max(ins_cluster.get_source_length(), del_cluster.get_length()))) for ind, del_cluster in enumerate(deletion_evidence_clusters)]
        closest_deletion = sorted(distances, key=lambda obj: obj[1])[0]
        source_contig, source_start, source_end = ins_cluster.get_source()
        dest_contig, dest_start, dest_end = ins_cluster.get_destination()
        if closest_deletion[1] <= 1.0:
            #Insertion
            all_members = ins_cluster.members + deletion_evidence_clusters[closest_deletion[0]].members
            insertion_candidates.append(CandidateInsertion(source_contig, source_start, source_end, dest_contig, dest_start, dest_end, all_members, ins_cluster.score))
            del(deletion_evidence_clusters[closest_deletion[0]])
        else:
            #Duplication
            int_duplication_candidates.append(CandidateDuplicationInterspersed(source_contig, source_start, source_end, dest_contig, dest_start, dest_end, ins_cluster.members, ins_cluster.score))
    return (insertion_candidates, int_duplication_candidates)


def get_closest_index(input_list, input_number):
    """
    Assumes input_list is sorted. Returns index of closest value to input_number.

    If two numbers are equally close, return the index of the smallest number.
    """
    if len(input_list) < 1:
        return None
    pos = bisect_left(input_list, input_number)
    if pos == 0:
        return 0
    if pos == len(input_list):
        return len(input_list) - 1
    before = input_list[pos - 1]
    after = input_list[pos]
    if after - input_number < input_number - before:
       return pos
    else:
       return pos - 1


def distance_positions(position1, position2):
    return float("inf") if position1[0] != position2[0] else abs(position1[1] - position2[1])


def cluster_positions_simple(positions, max_delta=1000):
    """Form partitions of (contig, position) pairs using simple distance."""
    sorted_positions = sorted(positions, key=lambda position: (position[0], position[1]))
    partitions = []
    current_partition = []
    for position in sorted_positions:
        if len(current_partition) < 1:
            current_partition.append(position)
            continue
        if distance_positions(current_partition[0], position) > max_delta:
            partitions.append(current_partition[:])
            while len(current_partition) > 0 and distance_positions(current_partition[0], position) > max_delta:
                current_partition.pop(0)
        current_partition.append(position)
    if len(current_partition) > 0:
        partitions.append(current_partition[:])
    return partitions


def calculate_score_deletion(main_score, translocation_distances, translocation_stds, translocation_deviation):
    """Calculate the score of a merged insertion detected from a deletion.
       Parameters: - main_score - score of the underlying main deletion
                   - translocation_distances - distance of the translocation clusters flanking the main deletion (left, right or both sides)
                   - translocation_stds - standard deviation of the translocation clusters flanking the main deletion (left, right or both sides)
                   - translocation_deviation - distance between left and right translocation destinations"""
    if len(translocation_distances) == 2:
        final_score = int(2 * main_score - 0.5 * sum(translocation_distances) - 0.3 * sum(translocation_stds) - 0.5 * translocation_deviation)
    if len(translocation_distances) == 1:
        final_score = int(main_score - 0.5 * sum(translocation_distances) - 0.3 * sum(translocation_stds) - 0.5 * translocation_deviation)
    # print("Score {0} from parameters: {1}, {2}, {3}, {4}".format(max(0, final_score), main_score, translocation_distances, translocation_stds, translocation_deviation))
    return max(0, final_score)


def calculate_score_insertion(main_score, translocation_distance, translocation_std, destination_stds):
    """Calculate the score of a merged insertion or duplication detected from an insertion.
       Parameters: - main_score - score of the underlying main insertion
                   - translocation_distance - distance of the translocation clusters flanking the main insertion (left)
                   - translocation_std - standard deviation of the translocation clusters flanking the main deletion (left)
                   - destination_stds - standard deviations of the left and right translocation destinations"""
    final_score = int(2 * main_score - 0.5 * translocation_distance - 0.2 * translocation_std - 0.1 * sum(destination_stds))
    # print("Score {0} from parameters: {1}, {2}, {3}, {4}".format(max(0, final_score), main_score, translocation_distance, translocation_std, destination_stds))
    return max(0, final_score)


def merge_translocations_at_deletions(translocation_evidences, deletion_evidence_clusters, parameters):
    insertion_candidates = []
    
    # Cluster translocations by contig and pos1
    translocation_partitions = form_partitions(translocation_evidences)
    translocation_partitions_dict = defaultdict(list)
    for partition in translocation_partitions:
        translocation_partitions_dict[partition[0].contig1].append(partition)
    
    # to_delete = []
    for deletion_index, del_cluster in enumerate(deletion_evidence_clusters):
        del_contig, del_start, del_end = del_cluster.get_source()
        translocation_partition_means = [int(round(sum([ev.pos1 for ev in partition]) / float(len(partition)))) for partition in translocation_partitions_dict[del_contig]]
        translocation_partition_std = [int(round(sqrt(sum([pow(abs(ev.pos1 - translocation_partition_means[index]), 2) for ev in partition]) / float(len(partition))))) for index, partition in enumerate(translocation_partitions_dict[del_contig])]
        if len(translocation_partition_means) < 1:
            continue
        closest_to_start = get_closest_index(translocation_partition_means, del_start)
        closest_to_end = get_closest_index(translocation_partition_means, del_end)
        # if translocations found close to start and end of deletion
        if abs(translocation_partition_means[closest_to_start] - del_start) <= parameters.max_translocation_distance and abs(translocation_partition_means[closest_to_end] - del_end) <= parameters.max_translocation_distance:
            destinations_from_start = cluster_positions_simple([(evidence.contig2, evidence.pos2) for evidence in translocation_partitions_dict[del_contig][closest_to_start]])
            destinations_from_end = cluster_positions_simple([(evidence.contig2, evidence.pos2) for evidence in translocation_partitions_dict[del_contig][closest_to_end]])
            # if translocations close to start and end point to only one destination each
            if len(destinations_from_start) == 1 and len(destinations_from_end) == 1:
                destination_from_start = (destinations_from_start[0][0][0], int(round(sum([pos for contig, pos in destinations_from_start[0]]) / float(len(destinations_from_start[0])))))
                destination_from_end = (destinations_from_end[0][0][0], int(round(sum([pos for contig, pos in destinations_from_end[0]]) / float(len(destinations_from_end[0])))))
                # if translocations close to start and end point to the same destination
                if destination_from_start[0] == destination_from_end[0] and abs(destination_from_start[1] - destination_from_end[1]) < 10:
                    mean_destination = int(round((destination_from_start[1] + destination_from_end[1]) / 2))
                    all_members = del_cluster.members + translocation_partitions_dict[del_contig][closest_to_start] + translocation_partitions_dict[del_contig][closest_to_end]
                    score = calculate_score_deletion(del_cluster.score, [abs(translocation_partition_means[closest_to_start] - del_start), abs(translocation_partition_means[closest_to_end] - del_end)], [translocation_partition_std[closest_to_start], translocation_partition_std[closest_to_end]], abs(destination_from_start[1] - destination_from_end[1]))
                    insertion_candidates.append(CandidateInsertion(del_contig, del_start, del_end, destination_from_start[0], mean_destination, mean_destination + (del_end - del_start), all_members, score))
                    # print("Found Insertion thanks to translocations at {0}:{1}-{2}".format(del_contig, del_start, del_end), file=sys.stderr)
                    # to_delete.append(deletion_index)
        # if translocations found close to start of deletion
        elif abs(translocation_partition_means[closest_to_start] - del_start) <= parameters.max_translocation_distance and del_cluster.score > 6:
            destinations_from_start = cluster_positions_simple([(evidence.contig2, evidence.pos2) for evidence in translocation_partitions_dict[del_contig][closest_to_start]])
            # if translocations close to start point to the same destination
            if len(destinations_from_start) == 1:
                destination_from_start = (destinations_from_start[0][0][0], int(round(sum([pos for contig, pos in destinations_from_start[0]]) / float(len(destinations_from_start[0])))))
                all_members = del_cluster.members + translocation_partitions_dict[del_contig][closest_to_start]
                score = calculate_score_deletion(del_cluster.score, [abs(translocation_partition_means[closest_to_start] - del_start)], [translocation_partition_std[closest_to_start]], 0)
                insertion_candidates.append(CandidateInsertion(del_contig, del_start, del_end, destination_from_start[0], destination_from_start[1], destination_from_start[1] + (del_end - del_start), all_members, score))
                # print("Found Insertion thanks to translocations at {0}:{1}-{2}".format(del_contig, del_start, del_end), file=sys.stderr)
                # to_delete.append(deletion_index)
        # if translocations found close to end of deletion
        elif abs(translocation_partition_means[closest_to_end] - del_end) <= parameters.max_translocation_distance and del_cluster.score > 6:
            destinations_from_end = cluster_positions_simple([(evidence.contig2, evidence.pos2) for evidence in translocation_partitions_dict[del_contig][closest_to_end]])
            # if translocations close to end point to the same destination
            if len(destinations_from_end) == 1:
                destination_from_end = (destinations_from_end[0][0][0], int(round(sum([pos for contig, pos in destinations_from_end[0]]) / float(len(destinations_from_end[0])))))
                all_members = del_cluster.members + translocation_partitions_dict[del_contig][closest_to_end]
                score = calculate_score_deletion(del_cluster.score, [abs(translocation_partition_means[closest_to_end] - del_end)], [translocation_partition_std[closest_to_end]], 0)
                insertion_candidates.append(CandidateInsertion(del_contig, del_start, del_end, destination_from_end[0], destination_from_end[1], destination_from_end[1] + (del_end - del_start), all_members, score))
                # print("Found Insertion thanks to translocations at {0}:{1}-{2}".format(del_contig, del_start, del_end), file=sys.stderr)
                # to_delete.append(deletion_index)

    # for deletion_index in to_delete[::-1]:
    #     del(deletion_evidence_clusters[deletion_index])

    return insertion_candidates


def merge_translocations_at_insertions(translocation_evidences, insertion_evidence_clusters, deletion_evidence_clusters, parameters):
    # Cluster translocations by contig and pos1
    translocation_partitions = form_partitions(translocation_evidences)
    translocation_partitions_dict = defaultdict(list)
    for partition in translocation_partitions:
        translocation_partitions_dict[partition[0].contig1].append(partition)
    
    # to_delete = []
    insertion_from_evidence_clusters = []
    for ins_cluster in insertion_evidence_clusters:
        ins_contig, ins_start, ins_end = ins_cluster.get_source()
        translocation_partition_means = [int(round(sum([ev.pos1 for ev in partition]) / float(len(partition)))) for partition in translocation_partitions_dict[ins_contig]]
        translocation_partition_std = [int(round(sqrt(sum([pow(abs(ev.pos1 - translocation_partition_means[index]), 2) for ev in partition]) / float(len(partition))))) for index, partition in enumerate(translocation_partitions_dict[ins_contig])]
        if len(translocation_partition_means) < 1:
            continue
        closest_to_start = get_closest_index(translocation_partition_means, ins_start)
        # if translocations found close to start of insertion
        if abs(translocation_partition_means[closest_to_start] - ins_start) <= parameters.max_translocation_distance:
            destinations_from_start = cluster_positions_simple([(evidence.contig2, evidence.pos2) for evidence in translocation_partitions_dict[ins_contig][closest_to_start]])
            # if translocations close to start point to two destinations
            if len(destinations_from_start) == 2:
                destination1 = (destinations_from_start[0][0][0], int(round(sum([pos for contig, pos in destinations_from_start[0]]) / float(len(destinations_from_start[0])))))
                destination1_std = int(round(sqrt(sum([pow(abs(pos - destination1[1]), 2) for contig, pos in destinations_from_start[0]]) / float(len(destinations_from_start[0])))))
                destination2 = (destinations_from_start[1][0][0], int(round(sum([pos for contig, pos in destinations_from_start[1]]) / float(len(destinations_from_start[1])))))
                destination2_std = int(round(sqrt(sum([pow(abs(pos - destination2[1]), 2) for contig, pos in destinations_from_start[1]]) / float(len(destinations_from_start[1])))))
                # if the two destinations have the right distance
                distance = abs(destination1[1] - destination2[1])
                if destination1[0] == destination2[0] and 0.95 <= ((ins_end - ins_start) / float(distance)) <= 1.1:
                    members = ins_cluster.members + translocation_partitions_dict[ins_contig][closest_to_start]
                    score = calculate_score_insertion(ins_cluster.score, abs(translocation_partition_means[closest_to_start] - ins_start), translocation_partition_std[closest_to_start], [destination1_std, destination2_std])
                    insertion_from_evidence_clusters.append(EvidenceClusterBiLocal(destination1[0], min(destination1[1], destination2[1]), max(destination1[1], destination2[1]), ins_contig, ins_start, ins_start + distance, score, len(members), members, "ins_dup"))
                    # to_delete.append(insertion_index)

    insertion_candidates = []
    int_duplication_candidates = []
    # print("INFO: Number of insertions/duplications detected from inserted regions and translocations:", len(insertion_from_evidence_clusters))
    for ins_dup_cluster in insertion_from_evidence_clusters:    
        distances = [(ind, del_cluster.gowda_diday_distance(ins_dup_cluster, max(ins_dup_cluster.get_source_length(), del_cluster.get_length()))) for ind, del_cluster in enumerate(deletion_evidence_clusters)]
        closest_deletion = sorted(distances, key=lambda obj: obj[1])[0]
        source_contig, source_start, source_end = ins_dup_cluster.get_source()
        dest_contig, dest_start, dest_end = ins_dup_cluster.get_destination()
        if closest_deletion[1] <= 1.0:
            #Insertion
            all_members = ins_dup_cluster.members + deletion_evidence_clusters[closest_deletion[0]].members
            # print("Insertion {0}:{1}-{2} -> {3}:{4}-{5}".format(source_contig, source_start, source_end, dest_contig, dest_start, dest_end))
            insertion_candidates.append(CandidateInsertion(source_contig, source_start, source_end, dest_contig, dest_start, dest_end, all_members, ins_dup_cluster.score))
        else:
            #Duplication
            # print("Duplication {0}:{1}-{2} -> {3}:{4}-{5}".format(source_contig, source_start, source_end, dest_contig, dest_start, dest_end))
            int_duplication_candidates.append(CandidateDuplicationInterspersed(source_contig, source_start, source_end, dest_contig, dest_start, dest_end, ins_dup_cluster.members, ins_dup_cluster.score))
    # for insertion_index in to_delete[::-1]:
    #     del(insertion_evidence_clusters[insertion_index])
    return (insertion_candidates, int_duplication_candidates)
