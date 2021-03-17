from __future__ import print_function

import sys
import logging

from random import seed, sample
from statistics import mean, stdev
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from edlib import align
from pysam import FastaFile

from svim.SVSignature import SignatureClusterUniLocal, SignatureClusterBiLocal
from svim.SVCandidate import CandidateDuplicationInterspersed


def form_partitions(sv_signatures, max_distance):
    """Form partitions of signatures using mean distance."""
    sorted_signatures = sorted(sv_signatures, key=lambda evi: evi.get_key())
    partitions = []
    current_partition = []
    for signature in sorted_signatures:
        if len(current_partition) > 0 and current_partition[-1].downstream_distance_to(signature) > max_distance:
            partitions.append(current_partition[:])
            current_partition = []
        current_partition.append(signature)
    if len(current_partition) > 0:
        partitions.append(current_partition[:])
    return partitions


def compute_haplotype_edit_distance(signature1, signature2, reference, window_padding = 100):
    window_start = min(signature1.start, signature2.start) - window_padding
    window_end = max(signature1.start, signature2.start) + window_padding

    #construct haplotype sequences for both signatures
    haplotype1 = reference.fetch(signature1.contig, max(0, window_start), max(0, signature1.start)).upper()
    haplotype1 += signature1.sequence.upper()
    haplotype1 += reference.fetch(signature1.contig, max(0, signature1.start), max(0, window_end)).upper()

    haplotype2 = reference.fetch(signature2.contig, max(0, window_start), max(0, signature2.start)).upper()
    haplotype2 += signature2.sequence.upper()
    haplotype2 += reference.fetch(signature2.contig, max(0, signature2.start), max(0, window_end)).upper()

    return align(haplotype1, haplotype2)["editDistance"]

def span_position_distance(signature1, signature2, signature_type, reference, distance_normalizer, cluster_max_distance):
    if signature_type == "DEL" or signature_type == "DUP_TAN":
        span1 = signature1.get_source()[2] - signature1.get_source()[1]
        span2 = signature2.get_source()[2] - signature2.get_source()[1]
        center1 = (signature1.get_source()[1] + signature1.get_source()[2]) // 2
        center2 = (signature2.get_source()[1] + signature2.get_source()[2]) // 2
        position_distance = abs(center1 - center2) / distance_normalizer
        span_distance = abs(span1 - span2) / max(span1, span2)
        return position_distance + span_distance
    elif signature_type == "INV": #two signatures from same read can be clustered together"
        span1 = signature1.get_source()[2] - signature1.get_source()[1]
        span2 = signature2.get_source()[2] - signature2.get_source()[1]
        center1 = (signature1.get_source()[1] + signature1.get_source()[2]) // 2
        center2 = (signature2.get_source()[1] + signature2.get_source()[2]) // 2
        position_distance = abs(center1 - center2) / distance_normalizer
        span_distance = abs(span1 - span2) / max(span1, span2)
        return position_distance + span_distance
    elif signature_type == "INS": #center is the insertion location
        span1 = signature1.get_source()[2] - signature1.get_source()[1]
        span2 = signature2.get_source()[2] - signature2.get_source()[1]
        center1 = signature1.get_source()[1]
        center2 = signature2.get_source()[1]
        position_distance = abs(center1 - center2) / distance_normalizer
        if position_distance > 2*cluster_max_distance:
            #do not compute edit distance if insertions are too distant
            span_distance = abs(span1 - span2) / max(span1, span2)
            return position_distance + span_distance
        else:
            edit_distance = compute_haplotype_edit_distance(signature1, signature2, reference)
            sequence_distance = edit_distance / max(span1, span2) / 2
            return position_distance + sequence_distance
    elif signature_type == "DUP_INT": #position distance is computed for source and destination
        span1 = signature1.get_source()[2] - signature1.get_source()[1]
        span2 = signature2.get_source()[2] - signature2.get_source()[1]
        source_center1 = (signature1.get_source()[1] + signature1.get_source()[2]) // 2
        source_center2 = (signature2.get_source()[1] + signature2.get_source()[2]) // 2
        position_distance_source = abs(source_center1 - source_center2) / distance_normalizer
        position_distance_destination = abs(signature1.get_destination()[1] - signature2.get_destination()[1]) / distance_normalizer
        span_distance = abs(span1 - span2) / max(span1, span2)
        return position_distance_source + position_distance_destination + span_distance
    elif signature_type == "BND": #only position distance is computed
        dist1 = abs(signature1.get_source()[1] - signature2.get_source()[1])
        dist2 = abs(signature1.get_destination()[1] - signature2.get_destination()[1])
        if signature1.direction1 == signature2.direction1 and signature1.direction2 == signature2.direction2:
            position_distance = (dist1 + dist2) / 3000
        else:
            position_distance = 99999
        return position_distance
    else:
        return None


def span_position_distance_clusters(cluster1, cluster2, distance_normalizer):
    "Span position distance function for merging clusters"
    span1 = cluster1.get_source()[2] - cluster1.get_source()[1]
    span2 = cluster2.get_source()[2] - cluster2.get_source()[1]
    center1 = (cluster1.get_source()[1] + cluster1.get_source()[2]) // 2
    center2 = (cluster2.get_source()[1] + cluster2.get_source()[2]) // 2
    position_distance = abs(center1 - center2) / distance_normalizer
    span_distance = abs(span1 - span2) / max(span1, span2)    
    return position_distance + span_distance


def span_position_distance_intdup_candidates(signature1, signature2, distance_normalizer):
    "Span position distance function for clustering candidates"
    span1 = signature1.get_source()[2] - signature1.get_source()[1]
    span2 = signature2.get_source()[2] - signature2.get_source()[1]
    source_center1 = (signature1.get_source()[1] + signature1.get_source()[2]) // 2
    source_center2 = (signature2.get_source()[1] + signature2.get_source()[2]) // 2
    position_distance_source = abs(source_center1 - source_center2) / distance_normalizer
    position_distance_destination = abs(signature1.get_destination()[1] - signature2.get_destination()[1]) / distance_normalizer
    span_distance = abs(span1 - span2) / max(span1, span2)
    return position_distance_source + position_distance_destination + span_distance


def clusters_from_partitions(partitions, reference, options):
    """Finds clusters in partitions using span-position distance and hierarchical clustering. 
    Assumes that all signatures in the given partition are of the same type and on the same contig"""
    clusters_final = []
    large_partitions = 0
    duplicate_signatures = 0
    #initialize random number generator with fixed number to produce same output from same input
    seed(1524)
    # Find clusters in each partition individually.
    for partition in partitions:
        if len(partition) > 100:
            partition_sample = sample(partition, 100)
            large_partitions += 1
        else:
            partition_sample = partition
        element_type = partition_sample[0].type
        assert(element_type in ["DEL", "DUP_TAN", "INV", "INS", "DUP_INT", "BND"])

        #remove similar signatures coming from the same read
        if element_type == "INV":
            #no duplication removal for inversions because they consist of two complementary signatures from the same read
            partition_sample_without_duplicates = partition_sample
        else:
            duplicates_from_same_read = set()
            for i in range(len(partition_sample)-1):
                for j in range(i+1, len(partition_sample)):
                    if partition_sample[i].read == partition_sample[j].read and span_position_distance(partition_sample[i], partition_sample[j], element_type, reference, options.distance_normalizer, options.cluster_max_distance) <= options.cluster_max_distance:
                        duplicates_from_same_read.add(j)
            duplicate_signatures += len(duplicates_from_same_read)
            partition_sample_without_duplicates = [partition_sample[i] for i in range(len(partition_sample)) if i not in duplicates_from_same_read]

        if len(partition_sample_without_duplicates) == 1:
            clusters_final.append([partition_sample_without_duplicates[0]])
            continue

        #compute pairwise distances
        distances = []
        if element_type == "INV":
            for i in range(len(partition_sample_without_duplicates)-1):
                for j in range(i+1, len(partition_sample_without_duplicates)):
                    distances.append(span_position_distance(partition_sample_without_duplicates[i], partition_sample_without_duplicates[j], element_type, reference, options.distance_normalizer, options.cluster_max_distance))
        else:
            for i in range(len(partition_sample_without_duplicates)-1):
                for j in range(i+1, len(partition_sample_without_duplicates)):
                    if partition_sample_without_duplicates[i].read == partition_sample_without_duplicates[j].read:
                        distances.append(99999)
                    else:
                        distances.append(span_position_distance(partition_sample_without_duplicates[i], partition_sample_without_duplicates[j], element_type, reference, options.distance_normalizer, options.cluster_max_distance))
        Z = linkage(np.array(distances), method = "average")
        cluster_indices = list(fcluster(Z, options.cluster_max_distance, criterion='distance'))
        new_clusters = [[] for i in range(max(cluster_indices))]
        for signature_index, cluster_index in enumerate(cluster_indices):
            new_clusters[cluster_index-1].append(partition_sample_without_duplicates[signature_index])
        clusters_final.extend(new_clusters)
    if len(partitions) > 0:
        if len(partitions[0]) > 0:
            logging.debug("%d out of %d partitions for %s exceeded 100 elements." % (large_partitions, len(partitions), partitions[0][0].type))
            logging.debug("%d %s signatures were removed due to similarity to another signature from the same read." % (duplicate_signatures, partitions[0][0].type))
    return clusters_final


def calculate_score(cluster, std_span, std_pos, span, type):
    if std_span == None or std_pos == None:
        span_deviation_score = 0
        pos_deviation_score = 0
    else:
        span_deviation_score = 1 - min(1, std_span / span)
        pos_deviation_score = 1 - min(1, std_pos / span)

    if type == "INV":
        directions = [signature.direction for signature in cluster]
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
        valid_signatures = min(left_signatures, right_signatures) + direction_counts[4]
        num_signatures = min(80, valid_signatures)
    else:
        num_signatures = min(80, len(cluster))
    return num_signatures + span_deviation_score * (num_signatures / 8) + pos_deviation_score * (num_signatures / 8)


def consolidate_clusters_unilocal(clusters):
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
        score = calculate_score(cluster, std_span, std_pos, average_end - average_start, cluster[0].type)
        consolidated_clusters.append(SignatureClusterUniLocal(cluster[0].get_source()[0], int(round(average_start)), int(round(average_end)), score, len(cluster), cluster, cluster[0].type, std_span, std_pos))
    return consolidated_clusters


def consolidate_clusters_bilocal(clusters):
    """Consolidate clusters to a list of (type, contig, mean start, mean end, cluster size, members) tuples."""
    consolidated_clusters = []
    for cluster in clusters:
        #Source
        source_average_start = sum([member.get_source()[1] for member in cluster]) / len(cluster)
        source_average_end = sum([member.get_source()[2] for member in cluster]) / len(cluster)
        if len(cluster) > 1:
            source_std_span = stdev([member.get_source()[2] - member.get_source()[1] for member in cluster])
            source_std_pos = stdev([(member.get_source()[2] + member.get_source()[1]) / 2 for member in cluster])
        else:
            source_std_span = None
            source_std_pos = None

        if cluster[0].type == "DUP_TAN":
            max_copies = max([member.copies for member in cluster])
            score = calculate_score(cluster, source_std_span, source_std_pos, source_average_end - source_average_start, cluster[0].type)
            consolidated_clusters.append(SignatureClusterBiLocal(cluster[0].get_source()[0],
                                                                int(round(source_average_start)),
                                                                int(round(source_average_end)),
                                                                cluster[0].get_source()[0],
                                                                int(round(source_average_end)),
                                                                int(round(source_average_end)) + max_copies *
                                                                (int(round(source_average_end)) -
                                                                 int(round(source_average_start))),
                                                                score, len(cluster), cluster, cluster[0].type, source_std_span, source_std_pos))
        elif cluster[0].type == "DUP_INT":
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
                score = calculate_score(cluster, None, None, mean([source_average_end - source_average_start, destination_average_end - destination_average_start]), cluster[0].type)
                consolidated_clusters.append(SignatureClusterBiLocal(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                                cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), score, len(cluster), cluster, cluster[0].type, None, None))
            else:
                score = calculate_score(cluster, mean([source_std_span, destination_std_span]), mean([source_std_pos, destination_std_pos]), mean([source_average_end - source_average_start, destination_average_end - destination_average_start]), cluster[0].type)
                consolidated_clusters.append(SignatureClusterBiLocal(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                                cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), score, len(cluster), cluster, cluster[0].type, mean([source_std_span, destination_std_span]), mean([source_std_pos, destination_std_pos])))
        elif cluster[0].type == "BND":
            #Destination
            destination_average_start = sum([member.get_destination()[1] for member in cluster]) / len(cluster)
            destination_average_end = sum([member.get_destination()[2] for member in cluster]) / len(cluster)
            #Directions
            directions1 = list(set([member.direction1 for member in cluster]))
            assert len(directions1) == 1
            direction1 = directions1[0]
            directions2 = list(set([member.direction2 for member in cluster]))
            assert len(directions2) == 1
            direction2 = directions2[0]

            if len(cluster) > 1:
                destination_std_pos = stdev([(member.get_destination()[2] + member.get_destination()[1]) / 2 for member in cluster])
            else:
                destination_std_span = None
                destination_std_pos = None
            if source_std_pos == None or destination_std_pos == None:
                score = calculate_score(cluster, None, None, 500, cluster[0].type)
                new_signature_cluster = SignatureClusterBiLocal(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                                cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), score, len(cluster), cluster, cluster[0].type, None, None)
            else:
                score = calculate_score(cluster, source_std_pos, destination_std_pos, 500, cluster[0].type)
                new_signature_cluster = SignatureClusterBiLocal(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                                cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), score, len(cluster), cluster, cluster[0].type, source_std_pos, destination_std_pos)
            new_signature_cluster.direction1 = direction1
            new_signature_cluster.direction2 = direction2
            consolidated_clusters.append(new_signature_cluster)
    return consolidated_clusters


def partition_and_cluster_candidates(candidates, options, type):
    partitions = form_partitions(candidates, options.partition_max_distance)
    clusters = []
    large_partitions = 0
    #initialize random number generator with fixed number to produce same output from same input
    seed(1524)
    # Find clusters in each partition individually.
    for partition in partitions:
        if len(partition) == 1:
            clusters.append([partition[0]])
            continue
        elif len(partition) > 100:
            partition_sample = sample(partition, 100)
            large_partitions += 1
        else:
            partition_sample = partition
        element_type = partition_sample[0].type
        distances = []
        for i in range(len(partition_sample)-1):
            for j in range(i+1, len(partition_sample)):
                distances.append(span_position_distance_intdup_candidates(partition_sample[i], partition_sample[j], options.distance_normalizer))
        Z = linkage(np.array(distances), method = "average")

        cluster_indices = list(fcluster(Z, options.cluster_max_distance, criterion='distance'))
        new_clusters = [[] for i in range(max(cluster_indices))]
        for signature_index, cluster_index in enumerate(cluster_indices):
            new_clusters[cluster_index-1].append(partition_sample[signature_index])
        clusters.extend(new_clusters)
    if len(partitions) > 0:
        if len(partitions[0]) > 0:
            logging.debug("%d out of %d partitions for %s exceeded 100 elements." % (large_partitions, len(partitions), partitions[0][0].type))
    logging.info("Clustered {0}: {1} partitions and {2} clusters".format(type, len(partitions), len(clusters)))

    final_candidates = []
    for cluster in clusters:
        combined_score = max([candidate.score for candidate in cluster])
        combined_members = [member for candidate in cluster for member in candidate.members]

        stds_span = [candidate.std_span for candidate in cluster if candidate.std_span != None]
        if len(stds_span) >= 1:
            combined_std_span = mean(stds_span)
        else:
            combined_std_span = None
        stds_pos = [candidate.std_pos for candidate in cluster if candidate.std_pos != None]
        if len(stds_pos) >= 1:
            combined_std_pos = mean(stds_pos)
        else:
            combined_std_pos = None

        #Source
        source_average_start = (sum([candidate.get_source()[1] for candidate in cluster]) / len(cluster))
        source_average_end = (sum([candidate.get_source()[2] for candidate in cluster]) / len(cluster))

        #Destination
        destination_average_start = (sum([candidate.get_destination()[1] for candidate in cluster]) / len(cluster))
        destination_average_end = (sum([candidate.get_destination()[2] for candidate in cluster]) / len(cluster))

        #Origin deleted?
        cutpaste = False
        for member in cluster:
            if member.cutpaste:
                cutpaste = True

        if cluster[0].type == "DUP_INT":
            final_candidates.append(CandidateDuplicationInterspersed(cluster[0].get_source()[0], int(round(source_average_start)), int(round(source_average_end)),
                                                       cluster[0].get_destination()[0], int(round(destination_average_start)), int(round(destination_average_end)), combined_members, combined_score, combined_std_span, combined_std_pos, cutpaste))
    return final_candidates


def partition_and_cluster(signatures, options, type):
    partitions = form_partitions(signatures, options.partition_max_distance)
    with FastaFile(options.genome) as reference:
        clusters = clusters_from_partitions(partitions, reference, options)
    logging.info("Clustered {0}: {1} partitions and {2} clusters".format(type, len(partitions), len(clusters)))
    if type == "deleted regions" or type == "inserted regions" or type == "inverted regions":
        return sorted(consolidate_clusters_unilocal(clusters), key=lambda cluster: (cluster.contig, (cluster.end + cluster.start) / 2))
    elif type == "tandem duplicated regions" or type == "inserted regions with detected region of origin" or type == "translocation breakpoints":
        return consolidate_clusters_bilocal(clusters)
    else:
        logging.error("Unknown parameter type={0} to function partition_and_cluster.")
