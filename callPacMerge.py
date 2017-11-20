from __future__ import print_function

import sys
from bisect import bisect_left
from collections import defaultdict

from SVEvidence import EvidenceTranslocation
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


def merge_translocations(translocation_evidences, deletion_evidence_clusters, insertion_evidence_clusters):
    insertion_candidates = []
    int_duplication_candidates = [] 
    
    # Generate a complete list of translocation by adding all reversed translocations
    reversed_translocations = []
    for evidence in translocation_evidences:
        reversed_translocations.append(EvidenceTranslocation(evidence.contig2, evidence.pos2, evidence.contig1, evidence.pos1, evidence.evidence, evidence.read))
    all_translocations = translocation_evidences + reversed_translocations
    
    # Cluster translocations by contig and pos1
    translocation_partitions = form_partitions(all_translocations)
    translocation_partitions_dict = defaultdict(list)
    for partition in translocation_partitions:
        translocation_partitions_dict[partition[0].contig1].append(partition)
    
    to_delete = []
    for deletion_index, del_cluster in enumerate(deletion_evidence_clusters):
        del_contig, del_start, del_end = del_cluster.get_source()
        translocation_partition_means = [int(round(sum([ev.pos1 for ev in partition]) / float(len(partition)))) for partition in translocation_partitions_dict[del_contig]]
        if len(translocation_partition_means) < 1:
            continue
        closest_to_start = get_closest_index(translocation_partition_means, del_start)
        closest_to_end = get_closest_index(translocation_partition_means, del_end)
        # if translocations found close to start and end of deletion
        if abs(translocation_partition_means[closest_to_start] - del_start) < 10 and abs(translocation_partition_means[closest_to_end] - del_end) < 10:
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
                    insertion_candidates.append(CandidateInsertion(del_contig, del_start, del_end, destination_from_start[0], mean_destination, mean_destination + (del_end - del_start), all_members, del_cluster.score))
                    #print("Found Insertion thanks to translocations at {0}:{1}-{2}".format(del_contig, del_start, del_end), file=sys.stderr)
                    to_delete.append(deletion_index)
        # if translocations found close to start of deletion
        elif abs(translocation_partition_means[closest_to_start] - del_start) < 10 and del_cluster.score > 6:
            destinations_from_start = cluster_positions_simple([(evidence.contig2, evidence.pos2) for evidence in translocation_partitions_dict[del_contig][closest_to_start]])
            # if translocations close to start point to the same destination
            if len(destinations_from_start) == 1:
                destination_from_start = (destinations_from_start[0][0][0], int(round(sum([pos for contig, pos in destinations_from_start[0]]) / float(len(destinations_from_start[0])))))
                all_members = del_cluster.members + translocation_partitions_dict[del_contig][closest_to_start]
                insertion_candidates.append(CandidateInsertion(del_contig, del_start, del_end, destination_from_start[0], destination_from_start[1], destination_from_start[1] + (del_end - del_start), all_members, del_cluster.score))
                #print("Found Insertion thanks to translocations at {0}:{1}-{2}".format(del_contig, del_start, del_end), file=sys.stderr)
                to_delete.append(deletion_index)
        # if translocations found close to end of deletion
        elif abs(translocation_partition_means[closest_to_end] - del_end) < 10 and del_cluster.score > 6:
            destinations_from_end = cluster_positions_simple([(evidence.contig2, evidence.pos2) for evidence in translocation_partitions_dict[del_contig][closest_to_end]])
            # if translocations close to end point to the same destination
            if len(destinations_from_end) == 1:
                destination_from_end = (destinations_from_end[0][0][0], int(round(sum([pos for contig, pos in destinations_from_end[0]]) / float(len(destinations_from_end[0])))))
                all_members = del_cluster.members + translocation_partitions_dict[del_contig][closest_to_end]
                insertion_candidates.append(CandidateInsertion(del_contig, del_start, del_end, destination_from_end[0], destination_from_end[1], destination_from_end[1] + (del_end - del_start), all_members, del_cluster.score))
                #print("Found Insertion thanks to translocations at {0}:{1}-{2}".format(del_contig, del_start, del_end), file=sys.stderr)
                to_delete.append(deletion_index)

    for deletion_index in to_delete[::-1]:
        del(deletion_evidence_clusters[deletion_index])
    
    return (insertion_candidates, int_duplication_candidates)