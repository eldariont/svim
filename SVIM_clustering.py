from __future__ import print_function

import sys
import logging

import networkx as nx
from random import sample

from SVEvidence import EvidenceClusterUniLocal, EvidenceClusterBiLocal
from SVCandidate import CandidateInsertion, CandidateDuplicationInterspersed


def form_partitions(sv_evidences, max_delta):
    """Form partitions of evidences using mean distance."""
    sorted_evidences = sorted(sv_evidences, key=lambda evi: evi.get_key())
    partitions = []
    current_partition = []
    for evidence in sorted_evidences:
        if len(current_partition) > 0 and current_partition[0].mean_distance_to(evidence) > max_delta:
            partitions.append(current_partition[:])
            current_partition = []
        current_partition.append(evidence)
    if len(current_partition) > 0:
        partitions.append(current_partition[:])
    return partitions


def clusters_from_partitions(partitions, parameters):
    """Form clusters in partitions using Gowda-Diday distance and clique finding in a distance graph."""
    clusters_full = []
    # Find clusters in each partition individually.
    for num, partition in enumerate(partitions):
        if len(partition) > 100:
            partition_sample = sample(partition, 100)
        else:
            partition_sample = partition
        largest_evidence = sorted(partition_sample, key=lambda evi: (evi.get_source()[2] - evi.get_source()[1]))[-1]
        largest_indel_size = largest_evidence.get_source()[2] - largest_evidence.get_source()[1]
        connection_graph = nx.Graph()
        connection_graph.add_nodes_from(range(len(partition_sample)))
        for i1 in range(len(partition_sample)):
            for i2 in range(len(partition_sample)):
                if i1 != i2:
                    if parameters["distance_metric"] == "gd":
                        if partition_sample[i1].gowda_diday_distance(partition_sample[i2], largest_indel_size) <= parameters["cluster_max_distance"]:
                            # Add edge in graph only if two indels are close to each other (distance <= max_delta)
                            connection_graph.add_edge(i1, i2)
                    else:
                        if partition_sample[i1].span_loc_distance(partition_sample[i2], parameters["distance_normalizer"]) <= parameters["cluster_max_distance"]:
                            # Add edge in graph only if two indels are close to each other (distance <= max_delta)
                            connection_graph.add_edge(i1, i2)
        clusters_indices = nx.find_cliques(connection_graph)
        for cluster in clusters_indices:
            clusters_full.append([partition_sample[index] for index in cluster])
    return clusters_full


def calculate_score(cigar_evidences, suppl_evidences):
    num_evidences = min(20, cigar_evidences) + min(20, suppl_evidences)
    evidence_boost = 0
    if cigar_evidences > 0:
        evidence_boost += 10
    if suppl_evidences > 0:
        evidence_boost += 20
    return num_evidences + evidence_boost


def calculate_score_inversion(inv_cluster, inversion_length, parameters):
    directions = [ev.direction for ev in inv_cluster]
    direction_counts = [0, 0, 0, 0, 0]
    for direction in directions:
        if direction == "left_fwd":
            direction_counts[0] += 1
        if direction == "left_rev":
            direction_counts[1] += 1
        if direction == "right_fwd":
            direction_counts[2] += 1
        if direction == "right_rev":
            direction_counts[3] += 1
        if direction == "all":
            direction_counts[4] += 1
    left_evidences = direction_counts[0] + direction_counts[1]
    right_evidences = direction_counts[2] + direction_counts[3]
    valid_suppl_evidences = min(left_evidences, right_evidences) + direction_counts[4]
    if inversion_length > parameters["max_sv_size"]:
        return 0
    else:
        return valid_suppl_evidences


def consolidate_clusters_unilocal(clusters, parameters):
    """Consolidate clusters to a list of (type, contig, mean start, mean end, cluster size, members) tuples."""
    consolidated_clusters = []
    for cluster in clusters:
        average_start = sum([member.get_source()[1] for member in cluster]) / len(cluster)
        average_end = sum([member.get_source()[2] for member in cluster]) / len(cluster)
        if cluster[0].type == "inv":
            score = calculate_score_inversion(cluster, average_end - average_start, parameters)
        else:
            cigar_evidences = [member for member in cluster if member.evidence == "cigar"]
            suppl_evidences = [member for member in cluster if member.evidence == "suppl"]
            score = calculate_score(len(cigar_evidences), len(suppl_evidences))
        consolidated_clusters.append(EvidenceClusterUniLocal(cluster[0].get_source()[0], int(round(average_start)), int(round(average_end)), score, len(cluster), cluster, cluster[0].type))
    return consolidated_clusters


def consolidate_clusters_bilocal(clusters):
    """Consolidate clusters to a list of (type, contig, mean start, mean end, cluster size, members) tuples."""
    consolidated_clusters = []
    for cluster in clusters:
        cigar_evidences = [member for member in cluster if member.evidence == "cigar"]
        suppl_evidences = [member for member in cluster if member.evidence == "suppl"]
        score = calculate_score(len(cigar_evidences), len(suppl_evidences))
        
        #Source
        source_average_start = sum([member.get_source()[1] for member in cluster]) / len(cluster)
        source_average_end = sum([member.get_source()[2] for member in cluster]) / len(cluster)

        if cluster[0].type == "dup":
            max_copies = max([member.copies for member in cluster])
            consolidated_clusters.append(EvidenceClusterBiLocal(cluster[0].get_source()[0],
                                                                int(round(source_average_start)),
                                                                int(round(source_average_end)),
                                                                cluster[0].get_source()[0],
                                                                int(round(source_average_end)),
                                                                int(round(source_average_end)) + max_copies *
                                                                (int(round(source_average_end)) -
                                                                 int(round(source_average_start))),
                                                                score, len(cluster), cluster, cluster[0].type))
        else:
            #Destination
            destination_average_start = sum([member.get_destination()[1] for member in cluster]) / len(cluster)
            destination_average_end = sum([member.get_destination()[2] for member in cluster]) / len(cluster)
            consolidated_clusters.append(EvidenceClusterBiLocal(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                                cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), score, len(cluster), cluster, cluster[0].type))
    return consolidated_clusters


def partition_and_cluster_candidates(candidates, parameters):
    partitions = form_partitions(candidates, parameters["partition_max_distance"])
    clusters = clusters_from_partitions(partitions, parameters)
    logging.info("Cluster results: {0} partitions and {1} clusters".format(len(partitions), len(clusters)))

    final_candidates = []
    for cluster in clusters:
        combined_score = sum([candidate.score for candidate in cluster])
        combined_members = [member for candidate in cluster for member in candidate.members]

        #Source
        source_average_start = (sum([candidate.get_source()[1] for candidate in cluster]) / len(cluster))
        source_average_end = (sum([candidate.get_source()[2] for candidate in cluster]) / len(cluster))

        #Destination
        destination_average_start = (sum([candidate.get_destination()[1] for candidate in cluster]) / len(cluster))
        destination_average_end = (sum([candidate.get_destination()[2] for candidate in cluster]) / len(cluster))

        if cluster[0].type == "ins":
            final_candidates.append(CandidateInsertion(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                       cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), combined_members, combined_score))
        elif cluster[0].type == "dup_int":
            final_candidates.append(CandidateDuplicationInterspersed(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                       cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), combined_members, combined_score))
    return final_candidates


def partition_and_cluster_unilocal(evidences, parameters):
    partitions = form_partitions(evidences, parameters["partition_max_distance"])
    clusters = clusters_from_partitions(partitions, parameters)
    logging.info("Cluster results: {0} partitions and {1} clusters".format(len(partitions), len(clusters)))
    return sorted(consolidate_clusters_unilocal(clusters, parameters), key=lambda cluster: (cluster.contig, (cluster.end + cluster.start) / 2))


def partition_and_cluster_bilocal(evidences, parameters):
    partitions = form_partitions(evidences, parameters["partition_max_distance"])
    clusters = clusters_from_partitions(partitions, parameters)
    logging.info("Cluster results: {0} partitions and {1} clusters".format(len(partitions), len(clusters)))
    return consolidate_clusters_bilocal(clusters)
