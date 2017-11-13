from __future__ import print_function

import sys

import networkx as nx

from SVEvidence import EvidenceClusterUniLocal, EvidenceClusterBiLocal


def form_partitions(sv_evidences, max_delta=1000):
    """Form partitions of evidences using mean distance."""
    sorted_evidences = sorted(sv_evidences, key=lambda evi: evi.get_key())
    partitions = []
    current_partition = []
    for evidence in sorted_evidences:
        if len(current_partition) < 1:
            current_partition.append(evidence)
            continue
        if current_partition[0].mean_distance_to(evidence) > max_delta:
            partitions.append(current_partition[:])
            while len(current_partition) > 0 and current_partition[0].mean_distance_to(evidence) > max_delta:
                current_partition.pop(0)
        current_partition.append(evidence)
    if len(current_partition) > 0:
        partitions.append(current_partition[:])
    return partitions


def clusters_from_partitions(partitions, max_delta=1):
    """Form clusters in partitions using Gowda-Diday distance and clique finding in a distance graph."""
    clusters_full = []
    # Find clusters in each partition individually.
    for num, partition in enumerate(partitions):
        largest_evidence = sorted(partition, key=lambda evi: (evi.end - evi.start))[-1]
        largest_indel_size = largest_evidence.end - largest_evidence.start
        connection_graph = nx.Graph()
        connection_graph.add_nodes_from(range(len(partition)))
        for i1 in range(len(partition)):
            for i2 in range(len(partition)):
                if i1 == i2 or partition[i1].gowda_diday_distance(partition[i2], largest_indel_size) > max_delta:
                    pass
                else:
                    # Add edge in graph only if two indels are close to each other (distance <= max_delta)
                    connection_graph.add_edge(i1, i2)
        clusters_indices = nx.find_cliques(connection_graph)
        for cluster in clusters_indices:
            clusters_full.append([partition[index] for index in cluster])
    return clusters_full


def consolidate_clusters_unilocal(clusters):
    """Consolidate clusters to a list of (type, contig, mean start, mean end, cluster size, members) tuples."""
    consolidated_clusters = []
    for cluster in clusters:
        cigar_evidences = [member for member in cluster if member.evidence == "cigar"]
        kmer_evidences = [member for member in cluster if member.evidence == "kmer"]
        suppl_evidences = [member for member in cluster if member.evidence == "suppl"]
        score = len(cigar_evidences) + len(kmer_evidences) + len(suppl_evidences)
        if len(cigar_evidences) > 0:
            score += 5
        if len(kmer_evidences) > 0:
            score += 5
        if len(suppl_evidences) > 0:
            score += 5
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
        score = len(cigar_evidences) + len(kmer_evidences) + len(suppl_evidences)
        if len(cigar_evidences) > 0:
            score += 5
        if len(kmer_evidences) > 0:
            score += 5
        if len(suppl_evidences) > 0:
            score += 5
        
        #Source
        source_average_start = (sum([member.get_source()[1] for member in cigar_evidences]) + sum([member.get_source()[1] for member in kmer_evidences]) + sum([member.get_source()[1] for member in suppl_evidences])) / float(len(cigar_evidences) + len(kmer_evidences) + len(suppl_evidences))
        source_average_end = (sum([member.get_source()[2] for member in cigar_evidences]) + sum([member.get_source()[2] for member in kmer_evidences]) + sum([member.get_source()[2] for member in suppl_evidences])) / float(len(cigar_evidences) + len(kmer_evidences) + len(suppl_evidences))
        
        #Destination
        destination_average_start = (sum([member.get_destination()[1] for member in cigar_evidences]) + sum([member.get_destination()[1] for member in kmer_evidences]) + sum([member.get_destination()[1] for member in suppl_evidences])) / float(len(cigar_evidences) + len(kmer_evidences) + len(suppl_evidences))
        destination_average_end = (sum([member.get_destination()[2] for member in cigar_evidences]) + sum([member.get_destination()[2] for member in kmer_evidences]) + sum([member.get_destination()[2] for member in suppl_evidences])) / float(len(cigar_evidences) + len(kmer_evidences) + len(suppl_evidences))
        
        consolidated_clusters.append(EvidenceClusterBiLocal(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                            cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), score, len(cluster), cluster, cluster[0].type))
    return consolidated_clusters


def partition_and_cluster_unilocal(evidences):
    partitions = form_partitions(evidences)
    print("Formed {0} partitions".format(len(partitions)), file=sys.stderr)
    clusters = clusters_from_partitions(partitions)
    print("Subdivided partitions into {0} clusters".format(len(clusters)), file=sys.stderr)
    return sorted(consolidate_clusters_unilocal(clusters), key=lambda cluster: (cluster.contig, (cluster.end + cluster.start) / 2))


def partition_and_cluster_bilocal(evidences):
    partitions = form_partitions(evidences)
    print("Formed {0} partitions".format(len(partitions)), file=sys.stderr)
    clusters = clusters_from_partitions(partitions)
    print("Subdivided partitions into {0} clusters".format(len(clusters)), file=sys.stderr)
    return consolidate_clusters_bilocal(clusters)
