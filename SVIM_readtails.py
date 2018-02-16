from __future__ import print_function

import sys
import logging
import itertools

from collections import defaultdict

from SVEvidence import EvidenceDeletion, EvidenceInsertion, EvidenceInversion, EvidenceTranslocation
from SVIM_kmercounting import find_svs


def check_indel_candidate_minus(left_tail, right_tail, contig, full_read, reference, parameters):
    left_ref_end = left_tail.reference_end
    left_q_end = left_tail.query_alignment_end
    right_ref_start = right_tail.reference_start
    right_q_start = right_tail.query_alignment_start

    read_snippet = str(full_read[parameters.tail_span - left_q_end : len(full_read) - right_q_start].upper())
    ref_snippet = str(reference[contig].seq[right_ref_start:left_ref_end].upper().reverse_complement())
    sv_results = find_svs(ref_snippet, read_snippet, parameters, debug = False)
    
    sv_evidences = []
    for typ, start, end in sv_results:
        if typ == "del":
            #print("Deletion detected: {0}:{1}-{2} (length {3})".format(contig, left_ref_end - end , left_ref_end - start, end - start), file=sys.stdout)
            sv_evidences.append(EvidenceDeletion(contig, left_ref_end - end, left_ref_end - start, "kmer", left_tail.query_name))
        if typ == "ins":
            #print("Insertion detected: {0}:{1}-{2} (length {3})".format(contig, left_ref_end - start, left_ref_end - start + (end - start), end - start), file=sys.stdout)
            sv_evidences.append(EvidenceInsertion(contig, left_ref_end - start, left_ref_end - start + (end - start), "kmer", left_tail.query_name))
    return sv_evidences


def check_indel_candidate_plus(left_tail, right_tail, contig, full_read, reference, parameters):
    left_ref_start = left_tail.reference_start
    left_q_start = left_tail.query_alignment_start
    right_ref_end = right_tail.reference_end
    right_q_end = right_tail.query_alignment_end
    
    read_snippet = str(full_read[left_q_start:len(full_read) - parameters.tail_span + right_q_end].upper())
    ref_snippet = str(reference[contig].seq[left_ref_start:right_ref_end].upper())
    sv_results = find_svs(ref_snippet, read_snippet, parameters, debug = False)

    sv_evidences = []
    for typ, start, end in sv_results:
        if typ == "del":
            #print("Deletion detected: {0}:{1}-{2} (length {3})".format(contig, left_ref_start + start, left_ref_start + end, end - start), file=sys.stdout)
            sv_evidences.append(EvidenceDeletion(contig, left_ref_start + start, left_ref_start + end, "kmer", left_tail.query_name))
        if typ == "ins":
            #print("Insertion detected: {0}:{1}-{2} (length {3})".format(contig, left_ref_start + start, left_ref_start + end, end - start), file=sys.stdout)
            sv_evidences.append(EvidenceInsertion(contig, left_ref_start + start, left_ref_start + end, "kmer", left_tail.query_name))
    return sv_evidences


def check_inv_1(left_tail, right_tail, contig, full_read, reference, parameters):
    left_ref_start = left_tail.reference_start
    left_q_start = left_tail.query_alignment_start
    right_ref_start = right_tail.reference_start
    right_q_start = right_tail.query_alignment_start
    
    read_snippet = str(full_read[left_q_start:len(full_read) - right_q_start].upper())
    ref_snippet_1 = str(reference[contig].seq[left_ref_start:left_ref_start+len(read_snippet)].upper())
    ref_snippet_2 = str(reference[contig].seq[right_ref_start:right_ref_start+len(read_snippet)].upper().reverse_complement())
    sv_results = find_svs(ref_snippet_1 + ref_snippet_2, read_snippet, parameters, debug = False)

    sv_evidences = []
    for typ, start, end in sv_results:
        if typ == "del" and start < len(ref_snippet_1) and end > len(ref_snippet_1):
            inv_start = left_ref_start + start
            inv_end = right_ref_start + (len(ref_snippet_1 + ref_snippet_2) - end)
            if inv_end > inv_start:
                #print("Inversion detected: {0}:{1}-{2} (length {3})".format(contig, inv_start, inv_end, inv_end - inv_start), file=sys.stdout)
                sv_evidences.append(EvidenceInversion(contig, inv_start, inv_end, "kmer", left_tail.query_name, "left_fwd"))
    return sv_evidences


def check_inv_2(left_tail, right_tail, contig, full_read, reference, parameters):
    left_ref_end = left_tail.reference_end
    left_q_end = left_tail.query_alignment_end
    right_ref_end = right_tail.reference_end
    right_q_end = right_tail.query_alignment_end
    
    read_snippet = str(full_read[parameters.tail_span - left_q_end : len(full_read) - parameters.tail_span + right_q_end].upper())
    ref_snippet_1 = str(reference[contig].seq[left_ref_end - len(read_snippet):left_ref_end].upper().reverse_complement())
    ref_snippet_2 = str(reference[contig].seq[right_ref_end - len(read_snippet):right_ref_end].upper())
    sv_results = find_svs(ref_snippet_1 + ref_snippet_2, read_snippet, parameters, debug = False)

    sv_evidences = []
    for typ, start, end in sv_results:
        if typ == "del" and start < len(ref_snippet_1) and end > len(ref_snippet_1):
            inv_start = left_ref_end - start
            inv_end = right_ref_end - (len(ref_snippet_1 + ref_snippet_2) - end)
            if inv_end > inv_start:
                #print("Inversion detected: {0}:{1}-{2} (length {3})".format(contig, inv_start, inv_end, inv_end - inv_start), file=sys.stdout)
                sv_evidences.append(EvidenceInversion(contig, inv_start, inv_end, "kmer", left_tail.query_name, "right_fwd"))
    return sv_evidences


def check_inv_3(left_tail, right_tail, contig, full_read, reference, parameters):
    left_ref_start = left_tail.reference_start
    left_q_start = left_tail.query_alignment_start
    right_ref_start = right_tail.reference_start
    right_q_start = right_tail.query_alignment_start
    
    read_snippet = str(full_read[left_q_start:len(full_read) - right_q_start].upper())
    ref_snippet_1 = str(reference[contig].seq[left_ref_start:left_ref_start+len(read_snippet)].upper())
    ref_snippet_2 = str(reference[contig].seq[right_ref_start:right_ref_start+len(read_snippet)].upper().reverse_complement())
    sv_results = find_svs(ref_snippet_1 + ref_snippet_2, read_snippet, parameters, debug = False)

    sv_evidences = []
    for typ, start, end in sv_results:
        if typ == "del" and start < len(ref_snippet_1) and end > len(ref_snippet_1):
            inv_start = right_ref_start + (len(ref_snippet_1 + ref_snippet_2) - end)
            inv_end = left_ref_start + start
            if inv_end > inv_start:
                #print("Inversion detected: {0}:{1}-{2} (length {3})".format(contig, inv_start, inv_end, inv_end - inv_start), file=sys.stdout)
                sv_evidences.append(EvidenceInversion(contig, inv_start, inv_end, "kmer", left_tail.query_name, "left_rev"))
    return sv_evidences


def check_inv_4(left_tail, right_tail, contig, full_read, reference, parameters):
    left_ref_end = left_tail.reference_end
    left_q_end = left_tail.query_alignment_end
    right_ref_end = right_tail.reference_end
    right_q_end = right_tail.query_alignment_end
    
    read_snippet = str(full_read[parameters.tail_span - left_q_end : len(full_read) - parameters.tail_span + right_q_end].upper())
    ref_snippet_1 = str(reference[contig].seq[left_ref_end - len(read_snippet):left_ref_end].upper().reverse_complement())
    ref_snippet_2 = str(reference[contig].seq[right_ref_end - len(read_snippet):right_ref_end].upper())
    sv_results = find_svs(ref_snippet_1 + ref_snippet_2, read_snippet, parameters, debug = False)
    sv_evidences = []
    for typ, start, end in sv_results:
        if typ == "del" and start < len(ref_snippet_1) and end > len(ref_snippet_1):
            inv_start = right_ref_end - (len(ref_snippet_1 + ref_snippet_2) - end)
            inv_end = left_ref_end - start
            if inv_end > inv_start:
                #print("Inversion detected: {0}:{1}-{2} (length {3})".format(contig, inv_start, inv_end, inv_end - inv_start), file=sys.stdout)
                sv_evidences.append(EvidenceInversion(contig, inv_start, inv_end, "kmer", left_tail.query_name, "right_rev"))
    return sv_evidences


def check_trans_candidate_1(left_tail, right_tail, left_contig, right_contig, full_read, reference, parameters):
    left_ref_start = left_tail.reference_start
    left_q_start = left_tail.query_alignment_start
    right_ref_end = right_tail.reference_end
    right_q_end = right_tail.query_alignment_end
    
    read_snippet = str(full_read[left_q_start:len(full_read) - parameters.tail_span + right_q_end].upper())
    ref_snippet_1 = str(reference[left_contig].seq[left_ref_start:left_ref_start+len(read_snippet)].upper())
    ref_snippet_2 = str(reference[right_contig].seq[right_ref_end - len(read_snippet):right_ref_end].upper())
    sv_results = find_svs(ref_snippet_1 + ref_snippet_2, read_snippet, parameters, debug = False)

    sv_evidences = []
    for typ, start, end in sv_results:
        if typ == "del" and start < len(ref_snippet_1) and end > len(ref_snippet_1):
            breakpoint1 = left_ref_start + start
            breakpoint2 = right_ref_end - (len(ref_snippet_1 + ref_snippet_2) - end)
            #print("Translocation breakpoint detected: {0}:{1} -> {2}:{3}".format(left_contig, breakpoint1, right_contig, breakpoint2), file=sys.stdout)
            sv_evidences.append(EvidenceTranslocation(left_contig, breakpoint1, right_contig, breakpoint2, "kmer", left_tail.query_name))
    return sv_evidences


def check_trans_candidate_2(left_tail, right_tail, left_contig, right_contig, full_read, reference, parameters):
    left_ref_end = left_tail.reference_end
    left_q_end = left_tail.query_alignment_end
    right_ref_start = right_tail.reference_start
    right_q_start = right_tail.query_alignment_start
    
    read_snippet = str(full_read[parameters.tail_span - left_q_end:len(full_read) - right_q_start].upper())
    ref_snippet_1 = str(reference[left_contig].seq[left_ref_end - len(read_snippet):left_ref_end].upper().reverse_complement())
    ref_snippet_2 = str(reference[right_contig].seq[right_ref_start:right_ref_start + len(read_snippet)].upper().reverse_complement())
    sv_results = find_svs(ref_snippet_1 + ref_snippet_2, read_snippet, parameters, debug = False)

    sv_evidences = []
    for typ, start, end in sv_results:
        if typ == "del" and start < len(ref_snippet_1) and end > len(ref_snippet_1):
            breakpoint1 = left_ref_end - start
            breakpoint2 = right_ref_start + (len(ref_snippet_1 + ref_snippet_2) - end)
            #print("Translocation breakpoint detected: {0}:{1} -> {2}:{3}".format(left_contig, breakpoint1, right_contig, breakpoint2), file=sys.stdout)
            sv_evidences.append(EvidenceTranslocation(left_contig, breakpoint1, right_contig, breakpoint2, "kmer", left_tail.query_name))
    return sv_evidences


def analyze_pair_of_read_tails(left_iterator_object, right_iterator_object, left_bam, right_bam, reads, reference, parameters):
    left_read_name, left_prim, left_suppl, left_sec = left_iterator_object
    right_read_name, right_prim, right_suppl, right_sec = right_iterator_object

    if len(left_prim) != 1 or left_prim[0].is_unmapped or left_prim[0].mapping_quality < parameters.min_mapq:
        return []
    if len(right_prim) != 1 or right_prim[0].is_unmapped or right_prim[0].mapping_quality < parameters.min_mapq:
        return []

    left_ref_chr = left_bam.getrname(left_prim[0].reference_id)
    left_ref_start = left_prim[0].reference_start
    left_ref_end = left_prim[0].reference_end
    left_q_start = left_prim[0].query_alignment_start
    left_q_end = left_prim[0].query_alignment_end

    right_ref_chr = right_bam.getrname(right_prim[0].reference_id)
    right_ref_start = right_prim[0].reference_start
    right_ref_end = right_prim[0].reference_end
    right_q_start = right_prim[0].query_alignment_start
    right_q_end = right_prim[0].query_alignment_end

    full_read = reads[left_read_name].seq
    read_length = len(full_read)

    if left_ref_chr == right_ref_chr:
        if left_prim[0].is_reverse and right_prim[0].is_reverse:
            reference_dist = left_ref_start - right_ref_end
            if reference_dist > 0:
                individual_dist = read_length - right_q_end - (parameters.tail_span - left_q_start)
                percent_shift = (individual_dist - reference_dist) / float(read_length)
                if percent_shift > parameters.tail_max_deviation or percent_shift < parameters.tail_min_deviation:
                    size_estimate = individual_dist - reference_dist - (0.04 * read_length)
                    if size_estimate > -parameters.max_deletion_size:
                        #INDEL candidate, check with k-mer counting
                        return check_indel_candidate_minus(left_prim[0], right_prim[0], left_ref_chr, full_read, reference, parameters)
                    else:
                        #Either very large DEL or TRANS
                        pass 
            else:
                #TRANS candidate
                return check_trans_candidate_2(left_prim[0], right_prim[0], left_ref_chr, right_ref_chr, full_read, reference, parameters)
        elif not left_prim[0].is_reverse and not right_prim[0].is_reverse:
            reference_dist = right_ref_start - left_ref_end
            if reference_dist > 0:
                individual_dist = read_length - left_q_end - (parameters.tail_span - right_q_start)
                percent_shift = (individual_dist - reference_dist) / float(read_length)
                if percent_shift > parameters.tail_max_deviation or percent_shift < parameters.tail_min_deviation:
                    size_estimate = individual_dist - reference_dist - (0.04 * read_length)
                    if size_estimate > -parameters.max_deletion_size:
                        #INDEL candidate, check with k-mer counting
                        return check_indel_candidate_plus(left_prim[0], right_prim[0], left_ref_chr, full_read, reference, parameters)
                    else:
                        #Either very large DEL or TRANS
                        pass
            else:
                #TRANS candidate
                return check_trans_candidate_1(left_prim[0], right_prim[0], left_ref_chr, right_ref_chr, full_read, reference, parameters)
        elif not left_prim[0].is_reverse and right_prim[0].is_reverse:
            if right_ref_start >= left_ref_start:
                #INV candidate, right tail in inverted region
                return check_inv_1(left_prim[0], right_prim[0], left_ref_chr, full_read, reference, parameters)
            else:
                #INV candidate, left tail in inverted region
                return check_inv_3(left_prim[0], right_prim[0], left_ref_chr, full_read, reference, parameters)
        elif left_prim[0].is_reverse and not right_prim[0].is_reverse:
            if right_ref_start >= left_ref_start:
                #INV candidate, left tail in inverted region
                return check_inv_2(left_prim[0], right_prim[0], left_ref_chr, full_read, reference, parameters)
            else:
                #INV candidate, right tail in inverted region
                return check_inv_4(left_prim[0], right_prim[0], left_ref_chr, full_read, reference, parameters)
    else:
        #TRANS candidate
        if left_prim[0].is_reverse and right_prim[0].is_reverse:
            return check_trans_candidate_2(left_prim[0], right_prim[0], left_ref_chr, right_ref_chr, full_read, reference, parameters)
        elif not left_prim[0].is_reverse and not right_prim[0].is_reverse:
            return check_trans_candidate_1(left_prim[0], right_prim[0], left_ref_chr, right_ref_chr, full_read, reference, parameters)
        elif not left_prim[0].is_reverse and right_prim[0].is_reverse:
            #TRANS + INV
            pass
        elif left_prim[0].is_reverse and not right_prim[0].is_reverse:
            #TRANS + INV
            pass
    return []


def confirm_del(left_bam, right_bam, evidence_cluster, reads, reference, parameters):
    contig, start, end = evidence_cluster.get_source()
    num_nearby_tails = left_bam.count(contig, start-10000, end+10000) + right_bam.count(contig, start-10000, end+10000)
    logging.info("Analyzed region {0}:{1}-{2} contains {3} tails".format(contig, start, end, num_nearby_tails))

    if num_nearby_tails > 1000:
        return

    tails = defaultdict(list)
    for left_tail in left_bam.fetch(contig, start-10000, end+10000):
        if not left_tail.is_unmapped and left_tail.mapping_quality >= parameters.min_mapq and not left_tail.is_secondary and not left_tail.is_supplementary:
            tails[left_tail.query_name].append(left_tail)
    for right_tail in right_bam.fetch(contig, start-10000, end+10000):
        if right_tail.query_name in tails:
            if not right_tail.is_unmapped and right_tail.mapping_quality >= parameters.min_mapq and not right_tail.is_secondary and not right_tail.is_supplementary:
                tails[right_tail.query_name].append(right_tail)

    evidences = []
    num_spanning_reads = 0
    for read_name, tails_pair in tails.iteritems():
        if len(tails_pair) != 2:
            continue
        left_tail = tails_pair[0]
        right_tail = tails_pair[1]

        left_ref_chr = left_bam.getrname(left_tail.reference_id)
        left_ref_start = left_tail.reference_start
        left_ref_end = left_tail.reference_end
        right_ref_start = right_tail.reference_start
        right_ref_end = right_tail.reference_end

        full_read = reads[read_name].seq

        if left_tail.is_reverse and right_tail.is_reverse and right_ref_end < start and left_ref_start > end:
            num_spanning_reads += 1
            evidences.extend(check_indel_candidate_minus(left_tail, right_tail, left_ref_chr, full_read, reference, parameters))
        elif not left_tail.is_reverse and not right_tail.is_reverse and left_ref_end < start and right_ref_start > end:
            num_spanning_reads += 1
            evidences.extend(check_indel_candidate_plus(left_tail, right_tail, left_ref_chr, full_read, reference, parameters))

    deletion_confirmations = [ev for ev in evidences if ev.type == "del" and evidence_cluster.gowda_diday_distance(ev, 10000) < 1]
    logging.info("Found {0} confirmations from {1} spanning reads".format(len(deletion_confirmations), num_spanning_reads))
    return (len(deletion_confirmations), num_spanning_reads)


def confirm_ins(left_bam, right_bam, evidence_cluster, reads, reference, parameters):
    contig, start, end = evidence_cluster.get_source()
    num_nearby_tails = left_bam.count(contig, start-10000, start+10000) + right_bam.count(contig, start-10000, start+10000)
    logging.info("Analyzed region {0}:{1}-{2} contains {3} tails".format(contig, start, end, num_nearby_tails))

    if num_nearby_tails > 1000:
        return

    tails = defaultdict(list)
    for left_tail in left_bam.fetch(contig, start-10000, start+10000):
        if not left_tail.is_unmapped and left_tail.mapping_quality >= parameters.min_mapq and not left_tail.is_secondary and not left_tail.is_supplementary:
            tails[left_tail.query_name].append(left_tail)
    for right_tail in right_bam.fetch(contig, start-10000, start+10000):
        if right_tail.query_name in tails:
            if not right_tail.is_unmapped and right_tail.mapping_quality >= parameters.min_mapq and not right_tail.is_secondary and not right_tail.is_supplementary:
                tails[right_tail.query_name].append(right_tail)

    evidences = []
    num_spanning_reads = 0
    for read_name, tails_pair in tails.iteritems():
        if len(tails_pair) != 2:
            continue
        left_tail = tails_pair[0]
        right_tail = tails_pair[1]

        left_ref_chr = left_bam.getrname(left_tail.reference_id)
        left_ref_start = left_tail.reference_start
        left_ref_end = left_tail.reference_end
        right_ref_start = right_tail.reference_start
        right_ref_end = right_tail.reference_end

        full_read = reads[read_name].seq

        if left_tail.is_reverse and right_tail.is_reverse and right_ref_end < start and left_ref_start > start:
            num_spanning_reads += 1
            evidences.extend(check_indel_candidate_minus(left_tail, right_tail, left_ref_chr, full_read, reference, parameters))
        elif not left_tail.is_reverse and not right_tail.is_reverse and left_ref_end < start and right_ref_start > start:
            num_spanning_reads += 1
            evidences.extend(check_indel_candidate_plus(left_tail, right_tail, left_ref_chr, full_read, reference, parameters))

    insertion_confirmations = [ev for ev in evidences if ev.type == "ins" and evidence_cluster.gowda_diday_distance(ev, 10000) < 1]
    logging.info("Found {0} confirmations from {1} spanning reads".format(len(insertion_confirmations), num_spanning_reads))
    return (len(insertion_confirmations), num_spanning_reads)


def confirm_inv(left_bam, right_bam, evidence_cluster, reads, reference, parameters):
    contig, start, end = evidence_cluster.get_source()
    num_nearby_tails = left_bam.count(contig, start-10000, end) + right_bam.count(contig, start-10000, start+10000)
    num_nearby_tails += left_bam.count(contig, start, end + 10000) + right_bam.count(contig, start, end+10000)
    logging.info("Analyzed region {0}:{1}-{2} contains {3} tails".format(contig, start, end, num_nearby_tails))

    if num_nearby_tails > 1000:
        return

    tails1 = defaultdict(list)
    for left_tail in left_bam.fetch(contig, start-10000, end):
        if not left_tail.is_unmapped and left_tail.mapping_quality >= parameters.min_mapq and not left_tail.is_secondary and not left_tail.is_supplementary:
            tails1[left_tail.query_name].append(left_tail)
    for right_tail in right_bam.fetch(contig, start-10000, end):
        if right_tail.query_name in tails1:
            if not right_tail.is_unmapped and right_tail.mapping_quality >= parameters.min_mapq and not right_tail.is_secondary and not right_tail.is_supplementary:
                tails1[right_tail.query_name].append(right_tail)
    tails2 = defaultdict(list)
    for left_tail in left_bam.fetch(contig, start, end + 10000):
        if not left_tail.is_unmapped and left_tail.mapping_quality >= parameters.min_mapq and not left_tail.is_secondary and not left_tail.is_supplementary:
            tails2[left_tail.query_name].append(left_tail)
    for right_tail in right_bam.fetch(contig, start, end + 10000):
        if right_tail.query_name in tails2:
            if not right_tail.is_unmapped and right_tail.mapping_quality >= parameters.min_mapq and not right_tail.is_secondary and not right_tail.is_supplementary:
                tails2[right_tail.query_name].append(right_tail)

    evidences = []
    num_spanning_reads = 0
    for read_name, tails_pair in itertools.chain(tails1.iteritems(), tails2.iteritems()):
        if len(tails_pair) != 2:
            continue
        left_tail = tails_pair[0]
        right_tail = tails_pair[1]

        left_ref_chr = left_bam.getrname(left_tail.reference_id)
        left_ref_start = left_tail.reference_start
        right_ref_chr = right_bam.getrname(right_tail.reference_id)
        right_ref_start = right_tail.reference_start

        full_read = reads[read_name].seq

        if not left_tail.is_reverse and right_tail.is_reverse:
            if right_ref_start >= left_ref_start:
                #INV candidate, right tail in inverted region
                num_spanning_reads += 1
                evidences.extend(check_inv_1(left_tail, right_tail, left_ref_chr, full_read, reference, parameters))
            else:
                #INV candidate, left tail in inverted region
                num_spanning_reads += 1
                evidences.extend(check_inv_3(left_tail, right_tail, left_ref_chr, full_read, reference, parameters))
        elif left_tail.is_reverse and not right_tail.is_reverse:
            if right_ref_start >= left_ref_start:
                #INV candidate, left tail in inverted region
                num_spanning_reads += 1
                evidences.extend(check_inv_2(left_tail, right_tail, left_ref_chr, full_read, reference, parameters))
            else:
                #INV candidate, right tail in inverted region
                num_spanning_reads += 1
                evidences.extend(check_inv_4(left_tail, right_tail, left_ref_chr, full_read, reference, parameters))

    inversion_confirmations = [ev for ev in evidences if ev.type == "inv" and evidence_cluster.gowda_diday_distance(ev, 10000) < 1]
    logging.info("Found {0} confirmations from {1} spanning reads".format(len(inversion_confirmations), num_spanning_reads))
    return (len(inversion_confirmations), num_spanning_reads)
