from __future__ import print_function

import sys
import logging

import networkx as nx
from random import sample
from statistics import mean, stdev

from SVSignature import SignatureClusterUniLocal, SignatureClusterBiLocal
from SVCandidate import CandidateDuplicationInterspersed


def form_partitions(sv_signatures, max_delta):
    """Form partitions of signatures using mean distance."""
    sorted_signatures = sorted(sv_signatures, key=lambda evi: evi.get_key())
    partitions = []
    current_partition = []
    for signature in sorted_signatures:
        if len(current_partition) > 0 and current_partition[0].mean_distance_to(signature) > max_delta:
            partitions.append(current_partition[:])
            current_partition = []
        current_partition.append(signature)
    if len(current_partition) > 0:
        partitions.append(current_partition[:])
    return partitions


def clusters_from_partitions(partitions, options):
    """Form clusters in partitions using Gowda-Diday distance and clique finding in a distance graph."""
    clusters_full = []
    # Find clusters in each partition individually.
    for num, partition in enumerate(partitions):
        if len(partition) > 100:
            partition_sample = sample(partition, 100)
        else:
            partition_sample = partition
        largest_signature = sorted(partition_sample, key=lambda evi: (evi.get_source()[2] - evi.get_source()[1]))[-1]
        largest_indel_size = largest_signature.get_source()[2] - largest_signature.get_source()[1]
        connection_graph = nx.Graph()
        connection_graph.add_nodes_from(range(len(partition_sample)))
        for i1 in range(len(partition_sample)):
            for i2 in range(len(partition_sample)):
                if i1 != i2:
                    if options.distance_metric == "gd":
                        if partition_sample[i1].gowda_diday_distance(partition_sample[i2], largest_indel_size) <= options.cluster_max_distance:
                            # Add edge in graph only if two indels are close to each other (distance <= max_delta)
                            connection_graph.add_edge(i1, i2)
                    else:
                        if partition_sample[i1].span_loc_distance(partition_sample[i2], options.distance_normalizer) <= options.cluster_max_distance:
                            # Add edge in graph only if two indels are close to each other (distance <= max_delta)
                            connection_graph.add_edge(i1, i2)
        clusters_indices = nx.find_cliques(connection_graph)
        for cluster in clusters_indices:
            clusters_full.append([partition_sample[index] for index in cluster])
    return clusters_full


def calculate_score(cigar_signatures, suppl_signatures, std_span, std_pos, span):
    if std_span == None or std_pos == None:
        span_deviation_score = 0
        pos_deviation_score = 0
    else:
        span_deviation_score = 1 - min(1, std_span / span)
        pos_deviation_score = 1 - min(1, std_pos / span)

    num_signatures = min(20, cigar_signatures) + min(20, suppl_signatures)
    signature_boost = 0
    if cigar_signatures > 0:
        signature_boost += 10
    if suppl_signatures > 0:
        signature_boost += 20
    #return min(50, num_signatures) + span_deviation_score * 25 + pos_deviation_score * 25
    #return 2 * min(40, num_signatures) + span_deviation_score * 15 + pos_deviation_score * 5
    return num_signatures + signature_boost + span_deviation_score * 20 + pos_deviation_score * 10


def calculate_score_inversion(inv_cluster, std_span, std_pos, span):
    if std_span == None or std_pos == None:
        span_deviation_score = 0
        pos_deviation_score = 0
    else:
        span_deviation_score = 1 - min(1, std_span / span)
        pos_deviation_score = 1 - min(1, std_pos / span)

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
    left_signatures = direction_counts[0] + direction_counts[1]
    right_signatures = direction_counts[2] + direction_counts[3]
    valid_suppl_signatures = min(left_signatures, right_signatures) + direction_counts[4]

    return min(70, valid_suppl_signatures) + span_deviation_score * 20 + pos_deviation_score * 10


def consolidate_clusters_unilocal(clusters, options):
    """Consolidate clusters to a list of (type, contig, mean start, mean end, cluster size, members) tuples."""
    consolidated_clusters = []
    for cluster in clusters:
        average_start = sum([member.get_source()[1] for member in cluster]) / len(cluster)
        average_end = sum([member.get_source()[2] for member in cluster]) / len(cluster)
        if len(cluster) > 1:
            std_span = stdev([member.get_source()[2] - member.get_source()[1] for member in cluster])
            std_pos = stdev([(member.get_source()[2] + member.get_source()[1]) / 2 for member in cluster])
        else:
            std_span = None
            std_pos = None
        if cluster[0].type == "inv":
            score = calculate_score_inversion(cluster, std_span, std_pos, average_end - average_start)
        else:
            cigar_signatures = [member for member in cluster if member.signature == "cigar"]
            suppl_signatures = [member for member in cluster if member.signature == "suppl"]
            score = calculate_score(len(cigar_signatures), len(suppl_signatures), std_span, std_pos, average_end - average_start)
        consolidated_clusters.append(SignatureClusterUniLocal(cluster[0].get_source()[0], int(round(average_start)), int(round(average_end)), score, len(cluster), cluster, cluster[0].type, std_span, std_pos))
    return consolidated_clusters


def consolidate_clusters_bilocal(clusters):
    """Consolidate clusters to a list of (type, contig, mean start, mean end, cluster size, members) tuples."""
    consolidated_clusters = []
    for cluster in clusters:
        cigar_signatures = [member for member in cluster if member.signature == "cigar"]
        suppl_signatures = [member for member in cluster if member.signature == "suppl"]
        
        #Source
        source_average_start = sum([member.get_source()[1] for member in cluster]) / len(cluster)
        source_average_end = sum([member.get_source()[2] for member in cluster]) / len(cluster)
        if len(cluster) > 1:
            source_std_span = stdev([member.get_source()[2] - member.get_source()[1] for member in cluster])
            source_std_pos = stdev([(member.get_source()[2] + member.get_source()[1]) / 2 for member in cluster])
        else:
            source_std_span = None
            source_std_pos = None

        if cluster[0].type == "dup":
            max_copies = max([member.copies for member in cluster])
            score = calculate_score(len(cigar_signatures), len(suppl_signatures), source_std_span, source_std_pos, source_average_end - source_average_start)
            consolidated_clusters.append(SignatureClusterBiLocal(cluster[0].get_source()[0],
                                                                int(round(source_average_start)),
                                                                int(round(source_average_end)),
                                                                cluster[0].get_source()[0],
                                                                int(round(source_average_end)),
                                                                int(round(source_average_end)) + max_copies *
                                                                (int(round(source_average_end)) -
                                                                 int(round(source_average_start))),
                                                                score, len(cluster), cluster, cluster[0].type, source_std_span, source_std_pos))
        else:
            #Destination
            destination_average_start = sum([member.get_destination()[1] for member in cluster]) / len(cluster)
            destination_average_end = sum([member.get_destination()[2] for member in cluster]) / len(cluster)
            if len(cluster) > 1:
                destination_std_span = stdev([member.get_destination()[2] - member.get_destination()[1] for member in cluster])
                destination_std_pos = stdev([(member.get_destination()[2] + member.get_destination()[1]) / 2 for member in cluster])
            else:
                destination_std_span = None
                destination_std_pos = None
            if source_std_span == None or source_std_pos == None or destination_std_span == None or destination_std_pos == None:
                score = calculate_score(len(cigar_signatures), len(suppl_signatures), None, None, mean([source_average_end - source_average_start, destination_average_end - destination_average_start]))
                consolidated_clusters.append(SignatureClusterBiLocal(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                                cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), score, len(cluster), cluster, cluster[0].type, None, None))
            else:
                score = calculate_score(len(cigar_signatures), len(suppl_signatures), mean([source_std_span, destination_std_span]), mean([source_std_pos, destination_std_pos]), mean([source_average_end - source_average_start, destination_average_end - destination_average_start]))
                consolidated_clusters.append(SignatureClusterBiLocal(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                                cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), score, len(cluster), cluster, cluster[0].type, mean([source_std_span, destination_std_span]), mean([source_std_pos, destination_std_pos])))
    return consolidated_clusters


def partition_and_cluster_candidates(candidates, options, type):
    partitions = form_partitions(candidates, options.partition_max_distance)
    clusters = clusters_from_partitions(partitions, options)
    logging.info("Clustered {0}: {1} partitions and {2} clusters".format(type, len(partitions), len(clusters)))

    final_candidates = []
    for cluster in clusters:
        combined_score = max([candidate.score for candidate in cluster])
        combined_members = [member for candidate in cluster for member in candidate.members]
        try:
            combined_std_span = mean([candidate.std_span for candidate in cluster])
        except TypeError:
            combined_std_span = None
        try:
            combined_std_pos = mean([candidate.std_pos for candidate in cluster])
        except TypeError:
            combined_std_pos = None

        #Source
        source_average_start = (sum([candidate.get_source()[1] for candidate in cluster]) / len(cluster))
        source_average_end = (sum([candidate.get_source()[2] for candidate in cluster]) / len(cluster))

        #Destination
        destination_average_start = (sum([candidate.get_destination()[1] for candidate in cluster]) / len(cluster))
        destination_average_end = (sum([candidate.get_destination()[2] for candidate in cluster]) / len(cluster))

        if cluster[0].type == "dup_int":
            final_candidates.append(CandidateDuplicationInterspersed(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                       cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), combined_members, combined_score, combined_std_span, combined_std_pos))
    return final_candidates


def partition_and_cluster_unilocal(signatures, options, type):
    partitions = form_partitions(signatures, options.partition_max_distance)
    clusters = clusters_from_partitions(partitions, options)
    logging.info("Clustered {0}: {1} partitions and {2} clusters".format(type, len(partitions), len(clusters)))
    return sorted(consolidate_clusters_unilocal(clusters, options), key=lambda cluster: (cluster.contig, (cluster.end + cluster.start) / 2))


def partition_and_cluster_bilocal(signatures, options, type):
    partitions = form_partitions(signatures, options.partition_max_distance)
    clusters = clusters_from_partitions(partitions, options)
    logging.info("Clustered {0}: {1} partitions and {2} clusters".format(type, len(partitions), len(clusters)))
    return consolidate_clusters_bilocal(clusters)
