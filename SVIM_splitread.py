from __future__ import print_function

import sys

from SVEvidence import EvidenceDeletion, EvidenceInsertion, EvidenceInversion, EvidenceTranslocation, EvidenceDuplicationTandem, EvidenceInsertionFrom
from SVIM_clustering import consolidate_clusters_bilocal, clusters_from_partitions


def analyze_full_read_segments_two(primary_aln, supplementary_aln, full_bam, parameters):
    read_name = primary_aln.query_name
    primary_ref_chr = full_bam.getrname(primary_aln.reference_id)
    primary_ref_start = primary_aln.reference_start
    primary_ref_end = primary_aln.reference_end
    primary_q_start = primary_aln.query_alignment_start
    primary_q_end = primary_aln.query_alignment_end

    supplementary_ref_chr = full_bam.getrname(supplementary_aln.reference_id)
    supplementary_ref_start = supplementary_aln.reference_start
    supplementary_ref_end = supplementary_aln.reference_end
    supplementary_q_start = supplementary_aln.query_alignment_start
    supplementary_q_end = supplementary_aln.query_alignment_end

    if primary_ref_chr == supplementary_ref_chr:
        if (primary_aln.is_reverse and supplementary_aln.is_reverse) or (not primary_aln.is_reverse and not supplementary_aln.is_reverse):
            if supplementary_q_start - primary_q_end >= -parameters["segment_overlap_tolerance"]:
                individual_dist = supplementary_q_start - primary_q_end
                reference_dist = supplementary_ref_start - primary_ref_end
                if reference_dist >= 0:                    
                    deviation = individual_dist - reference_dist
                    if deviation >= parameters["min_length"]:
                        #INS candidate
                        if reference_dist <= parameters["max_segment_gap_tolerance"]:
                            #print("Insertion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_end, primary_ref_end + deviation, deviation), file=sys.stdout)
                            return EvidenceInsertion(primary_ref_chr, primary_ref_end, primary_ref_end + deviation, "suppl", read_name)
                        else:
                            pass
                            #print("Insertion detected (imprecise): {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_end, primary_ref_end + deviation, deviation), file=sys.stdout)
                    elif -parameters["max_deletion_size"] <= deviation <= -parameters["min_length"]:
                        #DEL candidate
                        if individual_dist <= parameters["max_segment_gap_tolerance"]:
                            #print("Deletion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_end, primary_ref_end - deviation, -deviation), file=sys.stdout)
                            return EvidenceDeletion(primary_ref_chr, primary_ref_end, primary_ref_end - deviation, "suppl", read_name)
                        else:
                            pass
                            #print("Deletion detected (imprecise): {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_end, primary_ref_end - deviation, -deviation), file=sys.stdout)
                    elif deviation < -parameters["max_deletion_size"]:
                        #Either very large DEL or TRANS
                        #print("Translocation breakpoint detected: {0}:{1} -> {2}:{3}".format(primary_ref_chr, primary_ref_end, supplementary_ref_chr, supplementary_ref_start), file=sys.stdout)
                        return EvidenceTranslocation(primary_ref_chr, primary_ref_end, supplementary_ref_chr, supplementary_ref_start, "suppl", read_name)
                elif reference_dist < -parameters["min_length"] and supplementary_ref_end > primary_ref_start:
                    #Tandem Duplication
                    #print("Tandem duplication detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_start, primary_ref_end, primary_ref_end - supplementary_ref_start), file=sys.stdout)
                    return EvidenceDuplicationTandem(primary_ref_chr, supplementary_ref_start, primary_ref_end, 1, "suppl", read_name)
            elif primary_q_start - supplementary_q_end >= -parameters["segment_overlap_tolerance"]:
                reference_dist = primary_ref_start - supplementary_ref_end
                individual_dist = primary_q_start - supplementary_q_end
                if reference_dist >= 0:
                    deviation = individual_dist - reference_dist
                    if deviation >= parameters["min_length"]:
                        #INS candidate
                        if reference_dist <= parameters["max_segment_gap_tolerance"]:
                            #print("Insertion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_end, supplementary_ref_end + deviation, deviation), file=sys.stdout)
                            return EvidenceInsertion(primary_ref_chr, supplementary_ref_end, supplementary_ref_end + deviation, "suppl", read_name)
                        else:
                            pass
                            #print("Insertion detected (imprecise): {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_end, supplementary_ref_end + deviation, deviation), file=sys.stdout)
                    elif -parameters["max_deletion_size"] <= deviation <= -parameters["min_length"]:
                        #DEL candidate
                        if individual_dist <= parameters["max_segment_gap_tolerance"]:
                            #print("Deletion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_end, supplementary_ref_end - deviation, -deviation), file=sys.stdout)
                            return EvidenceDeletion(primary_ref_chr, supplementary_ref_end, supplementary_ref_end - deviation, "suppl", read_name)
                        else:
                            pass
                            #print("Deletion detected (imprecise): {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_end, supplementary_ref_end - deviation, -deviation), file=sys.stdout)
                    elif deviation < -parameters["max_deletion_size"]:
                        #Either very large DEL or TRANS
                        #print("Translocation breakpoint detected: {0}:{1} -> {2}:{3}".format(primary_ref_chr, primary_ref_end, supplementary_ref_chr, supplementary_ref_start), file=sys.stdout)
                        return EvidenceTranslocation(supplementary_ref_chr, supplementary_ref_end, primary_ref_chr, primary_ref_start, "suppl", read_name)
                elif reference_dist < -parameters["min_length"] and primary_ref_end > supplementary_ref_start:
                    #Tandem Duplication
                    #print("Tandem duplication detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_start, supplementary_ref_end, supplementary_ref_start - primary_ref_end), file=sys.stdout)
                    return EvidenceDuplicationTandem(primary_ref_chr, primary_ref_start, supplementary_ref_end, 1, "suppl", read_name)
            else:
                pass
                #print("Overlapping read segments in read", read_name)
        elif not primary_aln.is_reverse and supplementary_aln.is_reverse:
            supplementary_q_end = primary_aln.infer_read_length() - supplementary_aln.query_alignment_start
            supplementary_q_start = primary_aln.infer_read_length() - supplementary_aln.query_alignment_end
            if -parameters["segment_overlap_tolerance"] <= supplementary_q_start - primary_q_end <= parameters["max_segment_gap_tolerance"]:
                if supplementary_ref_start - primary_ref_end >= -parameters["segment_overlap_tolerance"]: # Case 1
                    #print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_end, supplementary_ref_end, supplementary_ref_end - primary_ref_end), file=sys.stdout)
                    return EvidenceInversion(primary_ref_chr, primary_ref_end, supplementary_ref_end, "suppl", read_name, "left_fwd")
                elif primary_ref_start - supplementary_ref_end >= -parameters["segment_overlap_tolerance"]: # Case 3
                    #print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_end, primary_ref_end, primary_ref_end - supplementary_ref_end), file=sys.stdout)
                    return EvidenceInversion(primary_ref_chr, supplementary_ref_end, primary_ref_end, "suppl", read_name, "left_rev")
            elif -parameters["segment_overlap_tolerance"] <= primary_q_start - supplementary_q_end <= parameters["max_segment_gap_tolerance"]:
                if primary_ref_start - supplementary_ref_end >= -parameters["segment_overlap_tolerance"]: # Case 2
                    #print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_start, primary_ref_start, primary_ref_start - supplementary_ref_start), file=sys.stdout)
                    return EvidenceInversion(primary_ref_chr, supplementary_ref_start, primary_ref_start, "suppl", read_name, "right_fwd")
                elif supplementary_ref_start - primary_ref_end >= -parameters["segment_overlap_tolerance"]: # Case 4
                    #print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_start, supplementary_ref_start, supplementary_ref_start - primary_ref_start), file=sys.stdout)
                    return EvidenceInversion(primary_ref_chr, primary_ref_start, supplementary_ref_start, "suppl", read_name, "right_rev")
            else:
                pass
                #print("Overlapping read segments in read", read_name)
        elif primary_aln.is_reverse and not supplementary_aln.is_reverse:
            primary_q_end = primary_aln.infer_read_length() - primary_aln.query_alignment_start
            primary_q_start = primary_aln.infer_read_length() - primary_aln.query_alignment_end
            if -parameters["segment_overlap_tolerance"] <= supplementary_q_start - primary_q_end <= parameters["max_segment_gap_tolerance"]:
                if supplementary_ref_start - primary_ref_end >= -parameters["segment_overlap_tolerance"]: # Case 2
                    #print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_start, supplementary_ref_start, supplementary_ref_start - primary_ref_start), file=sys.stdout)
                    return EvidenceInversion(primary_ref_chr, primary_ref_start, supplementary_ref_start, "suppl", read_name, "right_fwd")
                elif primary_ref_start - supplementary_ref_end >= -parameters["segment_overlap_tolerance"]: # Case 4
                    #print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_start, primary_ref_start, primary_ref_start - supplementary_ref_start), file=sys.stdout)
                    return EvidenceInversion(primary_ref_chr, supplementary_ref_start, primary_ref_start, "suppl", read_name, "right_rev")
            elif -parameters["segment_overlap_tolerance"] <= primary_q_start - supplementary_q_end <= parameters["max_segment_gap_tolerance"]:
                if primary_ref_start - supplementary_ref_end >= -parameters["segment_overlap_tolerance"]: # Case 1
                    #print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_end, primary_ref_end, primary_ref_end - supplementary_ref_end), file=sys.stdout)
                    return EvidenceInversion(primary_ref_chr, supplementary_ref_end, primary_ref_end, "suppl", read_name, "left_fwd")
                elif supplementary_ref_start - primary_ref_end >= -parameters["segment_overlap_tolerance"]: # Case 3
                    #print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_end, supplementary_ref_end, supplementary_ref_end - primary_ref_end), file=sys.stdout)
                    return EvidenceInversion(primary_ref_chr, primary_ref_end, supplementary_ref_end, "suppl", read_name, "left_rev")
            else:
                pass
                #print("Overlapping read segments in read", read_name)
    else:
        if (primary_aln.is_reverse and supplementary_aln.is_reverse) or (not primary_aln.is_reverse and not supplementary_aln.is_reverse):
            if -parameters["segment_overlap_tolerance"] <= supplementary_q_start - primary_q_end:
                if supplementary_q_start - primary_q_end <= parameters["max_segment_gap_tolerance"]:
                    #print("Translocation breakpoint detected: {0}:{1} -> {2}:{3}".format(primary_ref_chr, primary_ref_end, supplementary_ref_chr, supplementary_ref_start), file=sys.stdout)
                    return EvidenceTranslocation(primary_ref_chr, primary_ref_end, supplementary_ref_chr, supplementary_ref_start, "suppl", read_name)
                else:
                    pass
                    #print("Translocation breakpoint detected (imprecise): {0}:{1} -> {2}:{3}".format(primary_ref_chr, primary_ref_end, supplementary_ref_chr, supplementary_ref_start), file=sys.stdout)
            elif -parameters["segment_overlap_tolerance"] <= primary_q_start - supplementary_q_end:
                if primary_q_start - supplementary_q_end <= parameters["max_segment_gap_tolerance"]:
                    #print("Translocation breakpoint detected: {0}:{1} -> {2}:{3}".format(supplementary_ref_chr, supplementary_ref_end, primary_ref_chr, primary_ref_start), file=sys.stdout)
                    return EvidenceTranslocation(supplementary_ref_chr, supplementary_ref_end, primary_ref_chr, primary_ref_start, "suppl", read_name)
                    #print("Translocation breakpoint detected (imprecise): {0}:{1} -> {2}:{3}".format(supplementary_ref_chr, supplementary_ref_end, primary_ref_chr, primary_ref_start), file=sys.stdout)
                else:
                    pass
            else:
                pass
                #print("Overlapping read segments in read", read_name)
        else:
            #INV + TRANS
            pass
    return None


def analyze_full_read_segments_three(primary_aln, supplementary_aln1, supplementary_aln2, full_bam, parameters):
    read_name = primary_aln.query_name
    alns = [primary_aln, supplementary_aln1, supplementary_aln2]

    results = []
    ordered_alns = sorted(alns, key = lambda aln: aln.infer_read_length() - aln.query_alignment_end if aln.is_reverse else aln.query_alignment_start)
    ordered_alns_query_limits = [(aln.infer_read_length() - aln.query_alignment_end, aln.infer_read_length() - aln.query_alignment_start) if aln.is_reverse else (aln.query_alignment_start, aln.query_alignment_end) for aln in ordered_alns]
    ordered_alns_reference_names = [full_bam.getrname(aln.reference_id) for aln in ordered_alns]

    #all segments on same contig
    if ordered_alns_reference_names[0] == ordered_alns_reference_names[1] == ordered_alns_reference_names[2]:
        query_order_nice = True
        for i in xrange(len(ordered_alns) - 1):
            if ordered_alns_query_limits[i+1][0] - ordered_alns_query_limits[i][1] < -parameters["segment_overlap_tolerance"] or ordered_alns_query_limits[i+1][0] - ordered_alns_query_limits[i][1] > parameters["max_segment_gap_tolerance"]:
                query_order_nice = False
        
        reference_012 = True
        for i in xrange(len(ordered_alns) - 1):
            if ordered_alns[i+1].reference_start - ordered_alns[i].reference_end < -parameters["segment_overlap_tolerance"] or ordered_alns[i+1].reference_start - ordered_alns[i].reference_end > parameters["max_segment_gap_tolerance"]:
                reference_012 = False

        reference_210 = True
        for i in xrange(len(ordered_alns) - 1):
            if ordered_alns[i].reference_start - ordered_alns[i+1].reference_end < -parameters["segment_overlap_tolerance"] or ordered_alns[i].reference_start - ordered_alns[i+1].reference_end > parameters["max_segment_gap_tolerance"]:
                reference_210 = False
        
        if -parameters["segment_overlap_tolerance"] <= ordered_alns[2].reference_start - ordered_alns[0].reference_end <= parameters["max_segment_gap_tolerance"]:
            reference_02 = True
        else:
            reference_02 = False

        if -parameters["segment_overlap_tolerance"] <= ordered_alns[0].reference_start - ordered_alns[2].reference_end <= parameters["max_segment_gap_tolerance"]:
            reference_20 = True
        else:
            reference_20 = False

        if ordered_alns[1].reference_start < ordered_alns[0].reference_end and \
                        ordered_alns[1].reference_end > ordered_alns[0].reference_start:
            reference_01_overlap = True
        else:
            reference_01_overlap = False

        if ordered_alns[0].reference_start < ordered_alns[1].reference_end and \
                        ordered_alns[0].reference_end > ordered_alns[1].reference_start:
            reference_10_overlap = True
        else:
            reference_10_overlap = False

        if ordered_alns[2].reference_start < ordered_alns[1].reference_end and \
                        ordered_alns[2].reference_end > ordered_alns[1].reference_start:
            reference_12_overlap = True
        else:
            reference_12_overlap = False

        if ordered_alns[1].reference_start < ordered_alns[2].reference_end and \
                        ordered_alns[1].reference_end > ordered_alns[2].reference_start:
            reference_21_overlap = True
        else:
            reference_21_overlap = False

        # all segments come right after another
        if query_order_nice:
            # all segments are mapped right after another
            if reference_012 or reference_210:
                if not ordered_alns[0].is_reverse and ordered_alns[1].is_reverse and not ordered_alns[2].is_reverse:
                    #print("Inversion detected: {0}:{1}-{2} (length {3}, 3 segments)".format(ordered_alns_reference_names[0], ordered_alns[1].reference_start, ordered_alns[1].reference_end, ordered_alns[1].reference_end - ordered_alns[1].reference_start), file=sys.stdout)
                    results.append(EvidenceInversion(ordered_alns_reference_names[0], ordered_alns[1].reference_start, ordered_alns[1].reference_end, "suppl", read_name, "all"))
                elif ordered_alns[0].is_reverse and not ordered_alns[1].is_reverse and ordered_alns[2].is_reverse:
                    #print("Inversion detected: {0}:{1}-{2} (length {3}, 3 segments)".format(ordered_alns_reference_names[0], ordered_alns[1].reference_start, ordered_alns[1].reference_end, ordered_alns[1].reference_end - ordered_alns[1].reference_start), file=sys.stdout)
                    results.append(EvidenceInversion(ordered_alns_reference_names[0], ordered_alns[1].reference_start, ordered_alns[1].reference_end, "suppl", read_name, "all"))
            # segments 0 and 2 are mapped right after another and are from the + strand
            elif reference_02 and not ordered_alns[0].is_reverse and not ordered_alns[2].is_reverse:
                if ordered_alns[1].reference_start >= ordered_alns[2].reference_end or ordered_alns[1].reference_end <= ordered_alns[0].reference_start:
                    #Duplication or Insertion
                    #print("Duplication/Insertion detected from {0}:{1}-{2} to {0}:{3}".format(ordered_alns_reference_names[1], ordered_alns[1].reference_start, ordered_alns[1].reference_end, ordered_alns[0].reference_end), file=sys.stdout)
                    results.append(EvidenceInsertionFrom(ordered_alns_reference_names[1], ordered_alns[1].reference_start, ordered_alns[1].reference_end, ordered_alns_reference_names[0], ordered_alns[0].reference_end, "suppl", read_name))
            # segments 2 and 0 are mapped right after another and are from the - strand
            elif reference_20 and ordered_alns[0].is_reverse and ordered_alns[2].is_reverse:
                if ordered_alns[1].reference_start >= ordered_alns[0].reference_end or ordered_alns[1].reference_end <= ordered_alns[2].reference_start:
                    #Duplication or Insertion
                    #print("Duplication/Insertion detected from {0}:{1}-{2} to {0}:{3}".format(ordered_alns_reference_names[1], ordered_alns[1].reference_start, ordered_alns[1].reference_end, ordered_alns[2].reference_end), file=sys.stdout)
                    results.append(EvidenceInsertionFrom(ordered_alns_reference_names[1], ordered_alns[1].reference_start, ordered_alns[1].reference_end, ordered_alns_reference_names[2], ordered_alns[2].reference_end, "suppl", read_name))
            # all segments map to + strand
            elif not ordered_alns[0].is_reverse and not ordered_alns[1].is_reverse and not ordered_alns[2].is_reverse:
                if reference_01_overlap and reference_12_overlap:
                    overlaps = [EvidenceDuplicationTandem(ordered_alns_reference_names[0], ordered_alns[1].reference_start,
                                                      ordered_alns[0].reference_end, 1, "suppl", read_name),
                                EvidenceDuplicationTandem(ordered_alns_reference_names[0], ordered_alns[2].reference_start,
                                                      ordered_alns[1].reference_end, 1, "suppl", read_name)]
                    clustered_overlaps = consolidate_clusters_bilocal(clusters_from_partitions([overlaps], parameters["cluster_max_distance"]))
                    if len(clustered_overlaps) > 1:
                        results.extend(overlaps)
                    else:
                        contig, start, end = clustered_overlaps[0].get_source()
                        results.append(EvidenceDuplicationTandem(contig, start, end, 2, "suppl", read_name))
                elif reference_01_overlap:
                    results.append(EvidenceDuplicationTandem(ordered_alns_reference_names[0], ordered_alns[1].reference_start,
                                                      ordered_alns[0].reference_end, 1, "suppl", read_name))
                elif reference_12_overlap:
                    results.append(
                        EvidenceDuplicationTandem(ordered_alns_reference_names[0], ordered_alns[2].reference_start,
                                                  ordered_alns[1].reference_end, 1, "suppl", read_name))
            # all segments map to - strand
            elif ordered_alns[0].is_reverse and ordered_alns[1].is_reverse and ordered_alns[2].is_reverse:
                if reference_10_overlap and reference_21_overlap:
                    overlaps = [EvidenceDuplicationTandem(ordered_alns_reference_names[0], ordered_alns[0].reference_start,
                                                      ordered_alns[1].reference_end, 1, "suppl", read_name),
                                EvidenceDuplicationTandem(ordered_alns_reference_names[0], ordered_alns[1].reference_start,
                                                      ordered_alns[2].reference_end, 1, "suppl", read_name)]
                    clustered_overlaps = consolidate_clusters_bilocal(clusters_from_partitions([overlaps], parameters["cluster_max_distance"]))
                    if len(clustered_overlaps) > 1:
                        results.extend(overlaps)
                    else:
                        contig, start, end = clustered_overlaps[0].get_source()
                        results.append(EvidenceDuplicationTandem(contig, start, end, 2, "suppl", read_name))
                elif reference_10_overlap:
                    results.append(EvidenceDuplicationTandem(ordered_alns_reference_names[0], ordered_alns[0].reference_start,
                                                      ordered_alns[1].reference_end, 1, "suppl", read_name))
                elif reference_21_overlap:
                    results.append(
                        EvidenceDuplicationTandem(ordered_alns_reference_names[0], ordered_alns[1].reference_start,
                                                  ordered_alns[2].reference_end, 1, "suppl", read_name))
            else:
                pass
                #print("group01", read_name, file=sys.stderr)
        else:
            pass
            #print("group02", file=sys.stderr)

    elif ordered_alns_reference_names[0] == ordered_alns_reference_names[2]:
        query_order_nice = True
        for i in xrange(len(ordered_alns) - 1):
            if ordered_alns_query_limits[i+1][0] - ordered_alns_query_limits[i][1] < -parameters["segment_overlap_tolerance"] or ordered_alns_query_limits[i+1][0] - ordered_alns_query_limits[i][1] > parameters["max_segment_gap_tolerance"]:
                query_order_nice = False

        if -parameters["segment_overlap_tolerance"] <= ordered_alns[2].reference_start - ordered_alns[0].reference_end <= parameters["max_segment_gap_tolerance"]:
            reference_02 = True
        else:
            reference_02 = False

        if -parameters["segment_overlap_tolerance"] <= ordered_alns[0].reference_start - ordered_alns[2].reference_end <= parameters["max_segment_gap_tolerance"]:
            reference_20 = True
        else:
            reference_20 = False

        if query_order_nice:
            if reference_02:
                #Duplication or Insertion
                #print("Duplication/Insertion detected from {0}:{1}-{2} to {0}:{3}".format(ordered_alns_reference_names[1], ordered_alns[1].reference_start, ordered_alns[1].reference_end, ordered_alns[0].reference_end), file=sys.stdout)
                results.append(EvidenceInsertionFrom(ordered_alns_reference_names[1], ordered_alns[1].reference_start, ordered_alns[1].reference_end, ordered_alns_reference_names[0], ordered_alns[0].reference_end, "suppl", read_name))
            elif reference_20:
                #Duplication or Insertion
                #print("Duplication/Insertion detected from {0}:{1}-{2} to {0}:{3}".format(ordered_alns_reference_names[1], ordered_alns[1].reference_start, ordered_alns[1].reference_end, ordered_alns[2].reference_end), file=sys.stdout)
                results.append(EvidenceInsertionFrom(ordered_alns_reference_names[1], ordered_alns[1].reference_start, ordered_alns[1].reference_end, ordered_alns_reference_names[2], ordered_alns[2].reference_end, "suppl", read_name))
        else:
            pass
            #print("group03", file=sys.stderr)
    else:
        pass
        #print("group04", file=sys.stderr)
    return results


def analyze_full_read_segments(full_iterator_object, full_bam, parameters):
    """Analyze the positions and orientations of primary and supplementary alignments of a full read for signatures of SVs."""
    full_read_name, full_prim, full_suppl, full_sec = full_iterator_object

    sv_evidences = []

    if len(full_prim) != 1 or full_prim[0].is_unmapped or full_prim[0].mapping_quality < parameters["min_mapq"]:
        return sv_evidences

    good_suppl_alns = [aln for aln in full_suppl if not aln.is_unmapped and aln.mapping_quality >= parameters["min_mapq"]]
    if len(good_suppl_alns) == 1:
        result = analyze_full_read_segments_two(full_prim[0], good_suppl_alns[0], full_bam, parameters)
        if result != None:
            sv_evidences.append(result)
    elif len(good_suppl_alns) == 2:
        results = analyze_full_read_segments_three(full_prim[0], good_suppl_alns[0], good_suppl_alns[1], full_bam, parameters)
        sv_evidences.extend(results)
    elif 2 < len(good_suppl_alns) < 6:
        pass
        #print("Read", full_read_name.split("_")[1], "has", len(good_suppl_alns), "good supplementary segments.") 
    
    return sv_evidences