from __future__ import print_function

import sys
from bisect import bisect_left
from collections import defaultdict
from math import pow, sqrt

from svim.SVSignature import SignatureTranslocation, SignatureInsertionFrom, SignatureClusterBiLocal
from svim.SVCandidate import CandidateDuplicationInterspersed
from svim.SVIM_clustering import span_position_distance

def flag_cutpaste_candidates(insertion_from_signature_clusters, deletion_signature_clusters, options):
    """Flag duplication signature clusters if they overlap a deletion"""
    int_duplication_candidates = []
    for ins_cluster in insertion_from_signature_clusters:
        # Compute distances of every deletion cluster to the current insertion/duplication
        distances = [(del_index, span_position_distance((del_cluster.get_source()[1], del_cluster.get_source()[2], options.distance_normalizer), \
                                                        (ins_cluster.get_source()[1], ins_cluster.get_source()[2], options.distance_normalizer))) \
                                                        for del_index, del_cluster in enumerate(deletion_signature_clusters)]
        closest_deletion_index, closest_deletion = sorted(distances, key=lambda obj: obj[1])[0]
        source_contig, source_start, source_end = ins_cluster.get_source()
        dest_contig, dest_start, dest_end = ins_cluster.get_destination()
        # If close deletion cluster found
        if closest_deletion <= options.del_ins_dup_max_distance:
            #Potential cut&paste insertion
            int_duplication_candidates.append(CandidateDuplicationInterspersed(source_contig, source_start, source_end, dest_contig, dest_start, dest_end, ins_cluster.members, ins_cluster.score, ins_cluster.std_span, ins_cluster.std_pos, cutpaste=True))
        else:
            #Interspersed duplication
            int_duplication_candidates.append(CandidateDuplicationInterspersed(source_contig, source_start, source_end, dest_contig, dest_start, dest_end, ins_cluster.members, ins_cluster.score, ins_cluster.std_span, ins_cluster.std_pos, cutpaste=False))
    return int_duplication_candidates


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


def cluster_positions_simple(positions, max_delta):
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


def calculate_score_insertion(main_score, translocation_distances, translocation_stds, destination_stds):
    """Calculate the score of a merged insertion or duplication detected from an insertion.
       Parameters: - main_score - score of the underlying main insertion
                   - translocation_distances - mean distance of the translocation clusters flanking the main insertion (left)
                   - translocation_stds - standard deviation of the translocation clusters flanking the main insertion (left)
                   - destination_stds - standard deviations of the left and right translocation destinations"""
    #scale translocation distance to [0, 1] range
    td0 = max(0, 100 - translocation_distances[0]) / 100
    td1 = max(0, 100 - translocation_distances[1]) / 100

    #scale translocation std to [0, 1] range
    ts0 = max(0, 100 - translocation_stds[0]) / 100
    ts1 = max(0, 100 - translocation_stds[1]) / 100

    #scale destination stds to [0, 1] range
    ds0 = max(0, 100 - destination_stds[0]) / 100
    ds1 = max(0, 100 - destination_stds[1]) / 100

    #calculate final score as product of components
    product = td0 * td1 * ts0 * ts1 * ds0 * ds1
    final_score = pow(product, 1/6) * main_score
    return final_score


def merge_translocations_at_insertions(translocation_partitions_fwdfwd_dict, translocation_partition_means_fwdfwd_dict, translocation_partition_stds_fwdfwd_dict, translocation_partitions_revrev_dict, translocation_partition_means_revrev_dict, translocation_partition_stds_revrev_dict, insertion_signature_clusters, options):
    inserted_regions_to_remove = []
    insertion_from_signature_clusters = []
    for insertion_index, ins_cluster in enumerate(insertion_signature_clusters):
        ins_contig, ins_start, ins_end = ins_cluster.get_source()
        try:
            closest_to_start_fwdfwd_index = get_closest_index(translocation_partition_means_fwdfwd_dict[ins_contig], ins_start)
            closest_to_start_fwdfwd_mean = translocation_partition_means_fwdfwd_dict[ins_contig][closest_to_start_fwdfwd_index]
            closest_to_start_revrev_index = get_closest_index(translocation_partition_means_revrev_dict[ins_contig], ins_start)
            closest_to_start_revrev_mean = translocation_partition_means_revrev_dict[ins_contig][closest_to_start_revrev_index]
        except KeyError:
            continue
        # if translocations found close to start of insertion
        if abs(closest_to_start_fwdfwd_mean - ins_start) <= options.trans_sv_max_distance and abs(closest_to_start_revrev_mean - ins_start) <= options.trans_sv_max_distance:
            destinations_from_start_fwdfwd = cluster_positions_simple([(signature.contig2, signature.pos2) for signature in translocation_partitions_fwdfwd_dict[ins_contig][closest_to_start_fwdfwd_index]], options.trans_destination_partition_max_distance)
            destinations_from_start_revrev = cluster_positions_simple([(signature.contig2, signature.pos2) for signature in translocation_partitions_revrev_dict[ins_contig][closest_to_start_revrev_index]], options.trans_destination_partition_max_distance)
            # if translocations point to only one destination each
            if len(destinations_from_start_fwdfwd) == 1 and len(destinations_from_start_revrev) == 1:
                destination_from_start_fwdfwd = (destinations_from_start_fwdfwd[0][0][0], int(round(sum([pos for contig, pos in destinations_from_start_fwdfwd[0]]) / len(destinations_from_start_fwdfwd[0]))))
                destination_from_start_fwdfwd_std = int(round(sqrt(sum([pow(abs(pos - destination_from_start_fwdfwd[1]), 2) for contig, pos in destinations_from_start_fwdfwd[0]]) / len(destinations_from_start_fwdfwd[0]))))
                destination_from_start_revrev = (destinations_from_start_revrev[0][0][0], int(round(sum([pos for contig, pos in destinations_from_start_revrev[0]]) / len(destinations_from_start_revrev[0]))))
                destination_from_start_revrev_std = int(round(sqrt(sum([pow(abs(pos - destination_from_start_revrev[1]), 2) for contig, pos in destinations_from_start_revrev[0]]) / len(destinations_from_start_revrev[0]))))
                # if the two destinations have the right distance
                distance = abs(destination_from_start_revrev[1] - destination_from_start_fwdfwd[1])
                if destination_from_start_revrev[0] == destination_from_start_fwdfwd[0] and 0.95 <= ((ins_end - ins_start + 1) / (distance + 1)) <= 1.1:
                    members = ins_cluster.members + translocation_partitions_fwdfwd_dict[ins_contig][closest_to_start_fwdfwd_index] + translocation_partitions_revrev_dict[ins_contig][closest_to_start_revrev_index]
                    score = calculate_score_insertion(ins_cluster.score, 
                                                      [abs(closest_to_start_fwdfwd_mean - ins_start), abs(closest_to_start_revrev_mean - ins_start)], 
                                                      [translocation_partition_stds_fwdfwd_dict[ins_contig][closest_to_start_fwdfwd_index], translocation_partition_stds_revrev_dict[ins_contig][closest_to_start_revrev_index]], 
                                                      [destination_from_start_fwdfwd_std, destination_from_start_revrev_std])
                    insertion_from_signature_clusters.append(SignatureClusterBiLocal(destination_from_start_revrev[0], min(destination_from_start_revrev[1], destination_from_start_fwdfwd[1]), max(destination_from_start_revrev[1], destination_from_start_fwdfwd[1]), ins_contig, ins_start, ins_start + distance, score, len(members), members, "dup_int", ins_cluster.std_span, ins_cluster.std_pos))
                    inserted_regions_to_remove.append(insertion_index)

    return insertion_from_signature_clusters, inserted_regions_to_remove
