from __future__ import print_function

import sys
import logging

import networkx as nx
from random import sample

from SVEvidence import EvidenceClusterUniLocal, EvidenceClusterBiLocal
from SVCandidate import CandidateInsertion, CandidateDuplicationInterspersed


def form_partitions(sv_evidences, max_delta=1000):
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


def clusters_from_partitions(partitions, max_delta=1):
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
                if i1 == i2 or partition_sample[i1].gowda_diday_distance(partition_sample[i2], largest_indel_size) > max_delta:
                    pass
                else:
                    # Add edge in graph only if two indels are close to each other (distance <= max_delta)
                    connection_graph.add_edge(i1, i2)
        clusters_indices = nx.find_cliques(connection_graph)
        for cluster in clusters_indices:
            clusters_full.append([partition_sample[index] for index in cluster])
    return clusters_full


def calculate_score(cigar_evidences, kmer_evidences, suppl_evidences):
    num_evidences = cigar_evidences + kmer_evidences + suppl_evidences
    evidence_boost = 0
    if cigar_evidences > 0:
        evidence_boost += 6
    if kmer_evidences > 0:
        evidence_boost += 4
    if suppl_evidences > 0:
        evidence_boost += 10
    return num_evidences + evidence_boost


def consolidate_clusters_unilocal(clusters):
    """Consolidate clusters to a list of (type, contig, mean start, mean end, cluster size, members) tuples."""
    consolidated_clusters = []
    for cluster in clusters:
        cigar_evidences = [member for member in cluster if member.evidence == "cigar"]
        kmer_evidences = [member for member in cluster if member.evidence == "kmer"]
        suppl_evidences = [member for member in cluster if member.evidence == "suppl"]
        score = calculate_score(len(cigar_evidences), len(kmer_evidences), len(suppl_evidences))
        average_start = (2 * sum([member.get_source()[1] for member in cigar_evidences]) + sum([member.get_source()[1] for member in kmer_evidences]) + sum([member.get_source()[1] for member in suppl_evidences])) / float(2*len(cigar_evidences) + len(kmer_evidences) + len(suppl_evidences))
        average_end = (2 * sum([member.get_source()[2] for member in cigar_evidences]) + sum([member.get_source()[2] for member in kmer_evidences]) + sum([member.get_source()[2] for member in suppl_evidences])) / float(2*len(cigar_evidences) + len(kmer_evidences) + len(suppl_evidences))
        consolidated_clusters.append(EvidenceClusterUniLocal(cluster[0].get_source()[0], int(round(average_start)), int(round(average_end)), score, len(cluster), cluster, cluster[0].type))
    return consolidated_clusters


def consolidate_clusters_bilocal(clusters):
    """Consolidate clusters to a list of (type, contig, mean start, mean end, cluster size, members) tuples."""
    consolidated_clusters = []
    for cluster in clusters:
        cigar_evidences = [member for member in cluster if member.evidence == "cigar"]
        kmer_evidences = [member for member in cluster if member.evidence == "kmer"]
        suppl_evidences = [member for member in cluster if member.evidence == "suppl"]
        score = calculate_score(len(cigar_evidences), len(kmer_evidences), len(suppl_evidences))
        
        #Source
        source_average_start = (sum([member.get_source()[1] for member in cigar_evidences]) + sum([member.get_source()[1] for member in kmer_evidences]) + sum([member.get_source()[1] for member in suppl_evidences])) / float(len(cigar_evidences) + len(kmer_evidences) + len(suppl_evidences))
        source_average_end = (sum([member.get_source()[2] for member in cigar_evidences]) + sum([member.get_source()[2] for member in kmer_evidences]) + sum([member.get_source()[2] for member in suppl_evidences])) / float(len(cigar_evidences) + len(kmer_evidences) + len(suppl_evidences))

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
            destination_average_start = (sum([member.get_destination()[1] for member in cigar_evidences]) + sum([member.get_destination()[1] for member in kmer_evidences]) + sum([member.get_destination()[1] for member in suppl_evidences])) / float(len(cigar_evidences) + len(kmer_evidences) + len(suppl_evidences))
            destination_average_end = (sum([member.get_destination()[2] for member in cigar_evidences]) + sum([member.get_destination()[2] for member in kmer_evidences]) + sum([member.get_destination()[2] for member in suppl_evidences])) / float(len(cigar_evidences) + len(kmer_evidences) + len(suppl_evidences))

            consolidated_clusters.append(EvidenceClusterBiLocal(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                                cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), score, len(cluster), cluster, cluster[0].type))
    return consolidated_clusters


def partition_and_cluster_candidates(candidates):
    partitions = form_partitions(candidates)
    clusters = clusters_from_partitions(partitions)
    logging.info("Cluster results: {0} partitions and {1} clusters".format(len(partitions), len(clusters)))

    final_candidates = []
    for cluster in clusters:
        combined_score = sum([candidate.score for candidate in cluster])
        combined_members = [member for candidate in cluster for member in candidate.members]

        #Source
        source_average_start = (sum([candidate.get_source()[1] for candidate in cluster]) / float(len(cluster)))
        source_average_end = (sum([candidate.get_source()[2] for candidate in cluster]) / float(len(cluster)))

        #Destination
        destination_average_start = (sum([candidate.get_destination()[1] for candidate in cluster]) / float(len(cluster)))
        destination_average_end = (sum([candidate.get_destination()[2] for candidate in cluster]) / float(len(cluster)))

        if cluster[0].type == "ins":
            final_candidates.append(CandidateInsertion(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                       cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), combined_members, combined_score))
        elif cluster[0].type == "dup_int":
            final_candidates.append(CandidateDuplicationInterspersed(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                       cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), combined_members, combined_score))
    return final_candidates


def partition_and_cluster_unilocal(evidences):
    partitions = form_partitions(evidences)
    clusters = clusters_from_partitions(partitions)
    logging.info("Cluster results: {0} partitions and {1} clusters".format(len(partitions), len(clusters)))
    return sorted(consolidate_clusters_unilocal(clusters), key=lambda cluster: (cluster.contig, (cluster.end + cluster.start) / 2))


def partition_and_cluster_bilocal(evidences):
    partitions = form_partitions(evidences)
    clusters = clusters_from_partitions(partitions)
    logging.info("Cluster results: {0} partitions and {1} clusters".format(len(partitions), len(clusters)))
    return consolidate_clusters_bilocal(clusters)
