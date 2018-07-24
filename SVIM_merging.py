from __future__ import print_function

import sys
from bisect import bisect_left
from collections import defaultdict
from math import pow, sqrt

from SVSignature import SignatureTranslocation, SignatureInsertionFrom, SignatureClusterBiLocal
from SVCandidate import CandidateInsertion, CandidateDuplicationInterspersed

def merge_insertions_from(insertion_from_signature_clusters, deletion_signature_clusters, parameters):
    """Classify insertion/duplication signature clusters as insertion or duplication"""
    insertion_candidates = []
    int_duplication_candidates = []
    deleted_regions_to_remove = []
    for ins_cluster in insertion_from_signature_clusters:
        # Compute distances of every deletion cluster to the current insertion/duplication
        if parameters["distance_metric"] == "gd":
            distances = [(del_index, del_cluster.gowda_diday_distance(ins_cluster, max(ins_cluster.get_source_length(), del_cluster.get_length()))) for del_index, del_cluster in enumerate(deletion_signature_clusters)]
        elif parameters["distance_metric"] == "sl":
            distances = [(del_index, del_cluster.span_loc_distance(ins_cluster, parameters["distance_normalizer"])) for del_index, del_cluster in enumerate(deletion_signature_clusters)]
        closest_deletion_index, closest_deletion = sorted(distances, key=lambda obj: obj[1])[0]
        source_contig, source_start, source_end = ins_cluster.get_source()
        dest_contig, dest_start, dest_end = ins_cluster.get_destination()
        # If close deletion cluster found
        if closest_deletion <= parameters["del_ins_dup_max_distance"]:
            #Insertion
            all_members = ins_cluster.members + deletion_signature_clusters[closest_deletion_index].members
            insertion_candidates.append(CandidateInsertion(source_contig, source_start, source_end, dest_contig, dest_start, dest_end, all_members, ins_cluster.score, ins_cluster.std_span, ins_cluster.std_pos))
            deleted_regions_to_remove.append(closest_deletion_index)
        else:
            #Duplication
            int_duplication_candidates.append(CandidateDuplicationInterspersed(source_contig, source_start, source_end, dest_contig, dest_start, dest_end, ins_cluster.members, ins_cluster.score, ins_cluster.std_span, ins_cluster.std_pos))
    return (insertion_candidates, int_duplication_candidates, deleted_regions_to_remove)


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


def calculate_score_deletion(main_score, translocation_distances, translocation_stds, translocation_deviation):
    """Calculate the score of a merged insertion detected from a deletion.
       Parameters: - main_score - score of the underlying main deletion
                   - translocation_distances - mean distance of the translocation clusters flanking the main deletion (left, right or both sides)
                   - translocation_stds - standard deviation of the translocation clusters flanking the main deletion (left, right or both sides)
                   - translocation_deviation - distance between left and right translocation destinations"""
    if len(translocation_distances) == 2:
        #scale translocation distances to [0, 1] range
        td0 = max(0, 100 - translocation_distances[0]) / 100
        td1 = max(0, 100 - translocation_distances[1]) / 100

        #scale translocation stds to [0, 1] range
        ts0 = max(0, 100 - translocation_stds[0]) / 100
        ts1 = max(0, 100 - translocation_stds[1]) / 100

        #scale translocation deviation to [0, 1] range
        tv = max(0, 100 - translocation_deviation) / 100

        #calculate final score as product of components
        product = (main_score / 100) * td0 * td1 * ts0 * ts1 * tv
        final_score = pow(product, 1/6) * 100
    if len(translocation_distances) == 1:
        #scale translocation distance to [0, 1] range
        td0 = max(0, 100 - translocation_distances[0]) / 100

        #scale translocation std to [0, 1] range
        ts0 = max(0, 100 - translocation_stds[0]) / 100

        #scale translocation deviation to [0, 1] range
        tv = max(0, 100 - translocation_deviation) / 100

        #calculate final score as product of components (half score)
        product = (main_score / 100) * td0 * td0 * ts0 * ts0 * tv
        final_score = pow(product, 1/6) * 50
    # print("Score {0} from parameters: {1}, {2}, {3}, {4}".format(max(0, final_score), main_score, translocation_distances, translocation_stds, translocation_deviation))
    return final_score


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
    product = (main_score / 100) *  td0 * td1 * ts0 * ts1 * ds0 * ds1
    final_score = pow(product, 1/7) * 100
    # print("Score {0} from parameters: {1}, {2}, {3}, {4}".format(max(0, final_score), main_score, translocation_distance, translocation_std, destination_stds))
    return final_score


def merge_translocations_at_deletions(translocation_partitions_fwdfwd_dict, translocation_partition_means_fwdfwd_dict, translocation_partition_stds_fwdfwd_dict, translocation_partitions_revrev_dict, translocation_partition_means_revrev_dict, translocation_partition_stds_revrev_dict, deletion_signature_clusters, parameters):
    """Analyze all deletion signature clusters and look for flanking translocations that might indicate insertions."""
    insertion_candidates = []
    deleted_regions_to_remove = []
    for deletion_index, del_cluster in enumerate(deletion_signature_clusters):
        del_contig, del_start, del_end = del_cluster.get_source()

        try:
            closest_to_start_index = get_closest_index(translocation_partition_means_revrev_dict[del_contig], del_start)
            closest_to_start_mean = translocation_partition_means_revrev_dict[del_contig][closest_to_start_index]
            closest_to_end_index = get_closest_index(translocation_partition_means_fwdfwd_dict[del_contig], del_end)
            closest_to_end_mean = translocation_partition_means_fwdfwd_dict[del_contig][closest_to_end_index]
        except KeyError:
            continue
        # if translocations found close to start and end of deletion
        if abs(closest_to_start_mean - del_start) < abs(closest_to_start_mean - del_end) and \
           abs(closest_to_end_mean - del_end) < abs(closest_to_end_mean - del_start) and \
           abs(closest_to_start_mean - del_start) <= parameters["trans_sv_max_distance"] and \
           abs(closest_to_end_mean - del_end) <= parameters["trans_sv_max_distance"]:
            destinations_from_start = cluster_positions_simple([(signature.contig2, signature.pos2) for signature in translocation_partitions_revrev_dict[del_contig][closest_to_start_index]], parameters["trans_destination_partition_max_distance"])
            destinations_from_end = cluster_positions_simple([(signature.contig2, signature.pos2) for signature in translocation_partitions_fwdfwd_dict[del_contig][closest_to_end_index]], parameters["trans_destination_partition_max_distance"])
            # if translocations close to start and end point to only one destination each
            if len(destinations_from_start) == 1 and len(destinations_from_end) == 1:
                destination_from_start = (destinations_from_start[0][0][0], int(round(sum([pos for contig, pos in destinations_from_start[0]]) / len(destinations_from_start[0]))))
                destination_from_end = (destinations_from_end[0][0][0], int(round(sum([pos for contig, pos in destinations_from_end[0]]) / len(destinations_from_end[0]))))
                # if translocations close to start and end point to the same destination
                if destination_from_start[0] == destination_from_end[0] and abs(destination_from_start[1] - destination_from_end[1]) < parameters["trans_sv_max_distance"]:
                    mean_destination = int(round((destination_from_start[1] + destination_from_end[1]) / 2))
                    all_members = del_cluster.members + translocation_partitions_revrev_dict[del_contig][closest_to_start_index] + translocation_partitions_fwdfwd_dict[del_contig][closest_to_end_index]
                    score = calculate_score_deletion(del_cluster.score, 
                                                     [abs(closest_to_start_mean - del_start), abs(closest_to_end_mean - del_end)], 
                                                     [translocation_partition_stds_revrev_dict[del_contig][closest_to_start_index], 
                                                     translocation_partition_stds_fwdfwd_dict[del_contig][closest_to_end_index]], 
                                                     abs(destination_from_start[1] - destination_from_end[1]))
                    insertion_candidates.append(CandidateInsertion(del_contig, closest_to_start_mean, closest_to_end_mean, destination_from_start[0], mean_destination, mean_destination + (del_end - del_start), all_members, score, del_cluster.std_span, del_cluster.std_pos))
                    # print("Found Insertion thanks to translocations at {0}:{1}-{2}".format(del_contig, del_start, del_end), file=sys.stderr)
                    deleted_regions_to_remove.append(deletion_index)
        # if translocations found close to start of deletion
        elif abs(closest_to_start_mean - del_start) < abs(closest_to_start_mean - del_end) and \
             abs(closest_to_start_mean - del_start) <= abs(closest_to_end_mean - del_end) and \
             abs(closest_to_start_mean - del_start) <= parameters["trans_sv_max_distance"]:
            destinations_from_start = cluster_positions_simple([(signature.contig2, signature.pos2) for signature in translocation_partitions_revrev_dict[del_contig][closest_to_start_index]], parameters["trans_destination_partition_max_distance"])
            # if translocations close to start point to the same destination
            if len(destinations_from_start) == 1:
                destination_from_start = (destinations_from_start[0][0][0], int(round(sum([pos for contig, pos in destinations_from_start[0]]) / len(destinations_from_start[0]))))
                all_members = del_cluster.members + translocation_partitions_revrev_dict[del_contig][closest_to_start_index]
                score = calculate_score_deletion(del_cluster.score, [abs(closest_to_start_mean - del_start)], [translocation_partition_stds_revrev_dict[del_contig][closest_to_start_index]], 0)
                insertion_candidates.append(CandidateInsertion(del_contig, closest_to_start_mean, del_end, destination_from_start[0], destination_from_start[1], destination_from_start[1] + (del_end - del_start), all_members, score, del_cluster.std_span, del_cluster.std_pos))
                # print("Found Insertion thanks to translocations at {0}:{1}-{2}".format(del_contig, del_start, del_end), file=sys.stderr)
                deleted_regions_to_remove.append(deletion_index)
        # if translocations found close to end of deletion
        elif abs(closest_to_end_mean - del_end) < abs(closest_to_end_mean - del_start) and \
             abs(closest_to_end_mean - del_end) <= parameters["trans_sv_max_distance"]:
            destinations_from_end = cluster_positions_simple([(signature.contig2, signature.pos2) for signature in translocation_partitions_fwdfwd_dict[del_contig][closest_to_end_index]], parameters["trans_destination_partition_max_distance"])
            # if translocations close to end point to the same destination
            if len(destinations_from_end) == 1:
                destination_from_end = (destinations_from_end[0][0][0], int(round(sum([pos for contig, pos in destinations_from_end[0]]) / len(destinations_from_end[0]))))
                all_members = del_cluster.members + translocation_partitions_fwdfwd_dict[del_contig][closest_to_end_index]
                score = calculate_score_deletion(del_cluster.score, [abs(closest_to_end_mean - del_end)], [translocation_partition_stds_fwdfwd_dict[del_contig][closest_to_end_index]], 0)
                insertion_candidates.append(CandidateInsertion(del_contig, del_start, closest_to_end_mean, destination_from_end[0], destination_from_end[1], destination_from_end[1] + (del_end - del_start), all_members, score, del_cluster.std_span, del_cluster.std_pos))
                # print("Found Insertion thanks to translocations at {0}:{1}-{2}".format(del_contig, del_start, del_end), file=sys.stderr)
                deleted_regions_to_remove.append(deletion_index)

    return insertion_candidates, deleted_regions_to_remove


def merge_translocations_at_insertions(translocation_partitions_fwdfwd_dict, translocation_partition_means_fwdfwd_dict, translocation_partition_stds_fwdfwd_dict, translocation_partitions_revrev_dict, translocation_partition_means_revrev_dict, translocation_partition_stds_revrev_dict, insertion_signature_clusters, parameters):
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
        if abs(closest_to_start_fwdfwd_mean - ins_start) <= parameters["trans_sv_max_distance"] and abs(closest_to_start_revrev_mean - ins_start) <= parameters["trans_sv_max_distance"]:
            destinations_from_start_fwdfwd = cluster_positions_simple([(signature.contig2, signature.pos2) for signature in translocation_partitions_fwdfwd_dict[ins_contig][closest_to_start_fwdfwd_index]], parameters["trans_destination_partition_max_distance"])
            destinations_from_start_revrev = cluster_positions_simple([(signature.contig2, signature.pos2) for signature in translocation_partitions_revrev_dict[ins_contig][closest_to_start_revrev_index]], parameters["trans_destination_partition_max_distance"])
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
                    insertion_from_signature_clusters.append(SignatureClusterBiLocal(destination_from_start_revrev[0], min(destination_from_start_revrev[1], destination_from_start_fwdfwd[1]), max(destination_from_start_revrev[1], destination_from_start_fwdfwd[1]), ins_contig, ins_start, ins_start + distance, score, len(members), members, "ins_dup", ins_cluster.std_span, ins_cluster.std_pos))
                    inserted_regions_to_remove.append(insertion_index)

    return insertion_from_signature_clusters, inserted_regions_to_remove
