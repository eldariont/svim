from __future__ import print_function

import sys

import networkx as nx

from SVEvidence import EvidenceCluster


def mean_distance(evidence1, evidence2):
    """Return distance between means of two evidences."""
    if evidence1.contig == evidence2.contig and evidence1.type == evidence2.type:
        return abs(((evidence1.start + evidence1.end) / 2) - ((evidence2.start + evidence2.end) / 2))
    else:
        return float("inf")


def gowda_diday_distance(evidence1, evidence2, largest_indel_size):
    """Return Gowda-Diday distance between two evidences."""
    # different chromosomes
    if evidence1.contig != evidence2.contig:
        return float("inf")
    # different SV type
    if evidence1.type != evidence2.type:
        return float("inf")
    # non-intersecting
    if evidence1.end <= evidence2.start or evidence2.end <= evidence1.start:
        return float("inf")
    dist_pos = abs(evidence1.start - evidence2.start) / float(largest_indel_size)
    span1 = abs(evidence1.end - evidence1.start)
    span2 = abs(evidence2.end - evidence2.start)
    span_total = abs(max(evidence1.end, evidence2.end) - min(evidence1.start, evidence2.start))
    dist_span = abs(span1 - span2) / float(span_total)
    inter = min(evidence1.end, evidence2.end) - max(evidence1.start, evidence2.start)
    dist_content = (span1 + span2 - 2 * inter) / float(span_total)
    return dist_pos + dist_span + dist_content


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
                if i1 == i2 or gowda_diday_distance(partition[i1], partition[i2], largest_indel_size) > max_delta:
                    pass
                else:
                    # Add edge in graph only if two indels are close to each other (distance <= max_delta)
                    connection_graph.add_edge(i1, i2)
        clusters_indices = nx.find_cliques(connection_graph)
        for cluster in clusters_indices:
            clusters_full.append([partition[index] for index in cluster])
    return clusters_full


def consolidate_clusters(clusters):
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
        average_start = (2 * sum([member.start for member in cigar_evidences]) + sum([member.start for member in kmer_evidences]) + sum([member.start for member in suppl_evidences])) / float(2*len(cigar_evidences) + len(kmer_evidences) + len(suppl_evidences))
        average_end = (2 * sum([member.end for member in cigar_evidences]) + sum([member.end for member in kmer_evidences]) + sum([member.end for member in suppl_evidences])) / float(2*len(cigar_evidences) + len(kmer_evidences) + len(suppl_evidences))
        consolidated_clusters.append(EvidenceCluster(cluster[0].contig, int(round(average_start)), int(round(average_end)), score, len(cluster), cluster, cluster[0].type))
    return consolidated_clusters


def partition_and_cluster(evidences):
    partitions = form_partitions(evidences)
    print("Formed {0} partitions".format(len(partitions)), file=sys.stderr)
    clusters = clusters_from_partitions(partitions)
    print("Subdivided partitions into {0} clusters".format(len(clusters)), file=sys.stderr)
    return consolidate_clusters(clusters)
