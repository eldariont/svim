from __future__ import print_function

import sys

from SVEvidence import SVEvidence


def find_indels_in_cigar_tuples(tuples, min_length=50):
    """Parses CIGAR tuples (op, len) and returns Indels with a length > minLength"""
    pos = 0
    indels = []
    for operation, length in tuples:
        if operation == 0:                     # alignment match
            pos += length
        elif operation == 1:                   # insertion
            if length >= min_length:
                indels.append((pos, length, 'ins'))
        elif operation == 2:                   # deletion
            if length >= min_length:
                indels.append((pos, length, 'del'))
            pos += length
        elif operation == 7 or operation == 8:        # match or mismatch
            pos += length
    return indels


def analyze_one_supplementary(primary_aln, supplementary_aln, full_bam):
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
            if supplementary_q_start >= primary_q_end:
                individual_dist = supplementary_q_start - primary_q_end
                reference_dist = supplementary_ref_start - primary_ref_end
                if reference_dist >= 0:                    
                    deviation = individual_dist - reference_dist
                    if deviation > 50:
                        #INS candidate
                        if reference_dist <= 10:
                            print("Insertion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_end, primary_ref_end + deviation, deviation), file=sys.stdout)
                            return SVEvidence(primary_ref_chr, primary_ref_end, primary_ref_end + deviation, "ins", "suppl", read_name)
                        else:
                            print("Insertion detected (imprecise): {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_end, primary_ref_end + deviation, deviation), file=sys.stdout)
                    elif -10000 < deviation < -50:
                        #DEL candidate
                        if individual_dist <= 10:
                            print("Deletion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_end, primary_ref_end - deviation, -deviation), file=sys.stdout)
                            return SVEvidence(primary_ref_chr, primary_ref_end, primary_ref_end - deviation, "del", "suppl", read_name)
                        else:
                            print("Deletion detected (imprecise): {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_end, primary_ref_end - deviation, -deviation), file=sys.stdout)
                    elif deviation <= -10000:
                        #Either very large DEL or TRANS
                        pass
                elif reference_dist < -50 and supplementary_ref_start >= primary_ref_start:
                    #Tandem Duplication
                    print("Tandem duplication detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_start, primary_ref_end, primary_ref_end - supplementary_ref_start), file=sys.stdout)
                    return SVEvidence(primary_ref_chr, supplementary_ref_start, primary_ref_end, "dup_tan", "suppl", read_name)    
            elif primary_q_start >= supplementary_q_end:
                reference_dist = primary_ref_start - supplementary_ref_end
                individual_dist = primary_q_start - supplementary_q_end
                if reference_dist >= 0:
                    deviation = individual_dist - reference_dist
                    if deviation > 50:
                        #INS candidate
                        if reference_dist <= 10:
                            print("Insertion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_end, supplementary_ref_end + deviation, deviation), file=sys.stdout)
                            return SVEvidence(primary_ref_chr, supplementary_ref_end, supplementary_ref_end + deviation, "ins", "suppl", read_name)
                        else:
                            print("Insertion detected (imprecise): {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_end, supplementary_ref_end + deviation, deviation), file=sys.stdout)
                    elif -10000 < deviation < -50:
                        #DEL candidate
                        if individual_dist <= 10:
                            print("Deletion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_end, supplementary_ref_end - deviation, -deviation), file=sys.stdout)
                            return SVEvidence(primary_ref_chr, supplementary_ref_end, supplementary_ref_end - deviation, "del", "suppl", read_name)
                        else:
                            print("Deletion detected (imprecise): {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_end, supplementary_ref_end - deviation, -deviation), file=sys.stdout)
                    elif deviation <= -10000:
                        #Either very large DEL or TRANS
                        pass
                elif reference_dist < -50 and primary_ref_start >= supplementary_ref_start:
                    #Tandem Duplication
                    print("Tandem duplication detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_start, supplementary_ref_end, supplementary_ref_start - primary_ref_end), file=sys.stdout)
                    return SVEvidence(primary_ref_chr, primary_ref_start, supplementary_ref_end, "dup_tan", "suppl", read_name)    
            else:
                print("Overlapping read segments in read", read_name)
        elif not primary_aln.is_reverse and supplementary_aln.is_reverse:
            supplementary_q_end = primary_aln.infer_read_length() - supplementary_aln.query_alignment_start
            supplementary_q_start = primary_aln.infer_read_length() - supplementary_aln.query_alignment_end
            if supplementary_q_start - primary_q_end >= -5:
                individual_dist = supplementary_q_start - primary_q_end
                if individual_dist <= 10:
                    if supplementary_ref_start >= primary_ref_end: # Case 1
                        print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_end, supplementary_ref_end, supplementary_ref_end - primary_ref_end), file=sys.stdout)
                        return SVEvidence(primary_ref_chr, primary_ref_end, supplementary_ref_end, "inv", "suppl", read_name)
                    else: # Case 3
                        print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_end, primary_ref_end, primary_ref_end - supplementary_ref_end), file=sys.stdout)
                        return SVEvidence(primary_ref_chr, supplementary_ref_end, primary_ref_end, "inv", "suppl", read_name)
            elif primary_q_start - supplementary_q_end >= -5:
                individual_dist = primary_q_start - supplementary_q_end
                if individual_dist <= 10:
                    if primary_ref_start >= supplementary_ref_end: # Case 2
                        print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_start, primary_ref_start, primary_ref_start - supplementary_ref_start), file=sys.stdout)
                        return SVEvidence(primary_ref_chr, supplementary_ref_start, primary_ref_start, "inv", "suppl", read_name)
                    else: # Case 4
                        print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_start, supplementary_ref_start, supplementary_ref_start - primary_ref_start), file=sys.stdout)
                        return SVEvidence(primary_ref_chr, primary_ref_start, supplementary_ref_start, "inv", "suppl", read_name)
            elif supplementary_q_start >= primary_q_start and supplementary_q_end <= primary_q_end:
                print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_start, supplementary_ref_end, supplementary_ref_end - supplementary_ref_start), file=sys.stdout)
                return SVEvidence(primary_ref_chr, supplementary_ref_start, supplementary_ref_end, "inv", "suppl", read_name)
            else:
                print("Overlapping read segments in read", read_name)
        elif primary_aln.is_reverse and not supplementary_aln.is_reverse:
            primary_q_end = primary_aln.infer_read_length() - primary_aln.query_alignment_start
            primary_q_start = primary_aln.infer_read_length() - primary_aln.query_alignment_end
            if supplementary_q_start - primary_q_end > -5:
                individual_dist = supplementary_q_start - primary_q_end
                if individual_dist <= 10:
                    if supplementary_ref_start >= primary_ref_end: # Case 2
                        print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_start, supplementary_ref_start, supplementary_ref_start - primary_ref_start), file=sys.stdout)
                        return SVEvidence(primary_ref_chr, primary_ref_start, supplementary_ref_start, "inv", "suppl", read_name)
                    else: # Case 4
                        print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_start, primary_ref_start, primary_ref_start - supplementary_ref_start), file=sys.stdout)
                        return SVEvidence(primary_ref_chr, supplementary_ref_start, primary_ref_start, "inv", "suppl", read_name)
            elif primary_q_start - supplementary_q_end >= -5:
                individual_dist = primary_q_start - supplementary_q_end
                if individual_dist <= 10:
                    if primary_ref_start >= supplementary_ref_end: # Case 1
                        print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_end, primary_ref_end, primary_ref_end - supplementary_ref_end), file=sys.stdout)
                        return SVEvidence(primary_ref_chr, supplementary_ref_end, primary_ref_end, "inv", "suppl", read_name)
                    else: # Case 3
                        print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, primary_ref_end, supplementary_ref_end, supplementary_ref_end - primary_ref_end), file=sys.stdout)
                        return SVEvidence(primary_ref_chr, primary_ref_end, supplementary_ref_end, "inv", "suppl", read_name)
            elif supplementary_q_start >= primary_q_start and supplementary_q_end <= primary_q_end:
                print("Inversion detected: {0}:{1}-{2} (length {3})".format(primary_ref_chr, supplementary_ref_start, supplementary_ref_end, supplementary_ref_end - supplementary_ref_start), file=sys.stdout)
                return SVEvidence(primary_ref_chr, supplementary_ref_start, supplementary_ref_end, "inv", "suppl", read_name)
            else:
                print("Overlapping read segments in read", read_name)
    else:
        if (primary_aln.is_reverse and supplementary_aln.is_reverse) or (not primary_aln.is_reverse and not supplementary_aln.is_reverse):
            if supplementary_q_start >= primary_q_end:
                individual_dist = supplementary_q_start - primary_q_end
                if individual_dist <= 10:
                    print("Translocation breakpoint detected: {0}:{1} -> {2}:{3}".format(primary_ref_chr, primary_ref_end, supplementary_ref_chr, supplementary_ref_start), file=sys.stdout)
                    return SVEvidence(primary_ref_chr, primary_ref_end, supplementary_ref_start, "trans", "suppl", read_name, contig2 = supplementary_ref_chr)
                else:
                    print("Translocation breakpoint detected (imprecise): {0}:{1} -> {2}:{3}".format(primary_ref_chr, primary_ref_end, supplementary_ref_chr, supplementary_ref_start), file=sys.stdout)
            elif primary_q_start >= supplementary_q_end:
                individual_dist = primary_q_start - supplementary_q_end
                if individual_dist <= 10:
                    print("Translocation breakpoint detected: {0}:{1} -> {2}:{3}".format(supplementary_ref_chr, supplementary_ref_end, primary_ref_chr, primary_ref_start), file=sys.stdout)
                    return SVEvidence(supplementary_ref_chr, supplementary_ref_end, primary_ref_start, "trans", "suppl", read_name, contig2 = primary_ref_chr)
                    print("Translocation breakpoint detected (imprecise): {0}:{1} -> {2}:{3}".format(supplementary_ref_chr, supplementary_ref_end, primary_ref_chr, primary_ref_start), file=sys.stdout)
            else:
                print("Overlapping read segments in read", read_name)
        else:
            #INV + TRANS
            pass
    return None


def analyze_two_supplementary(primary_aln, supplementary_aln1, supplementary_aln2, full_bam):
    return None


def analyze_full_read(full_iterator_object, full_bam, parameters):
    full_read_name, full_prim, full_suppl, full_sec = full_iterator_object

    if len(full_prim) != 1 or full_prim[0].is_unmapped or full_prim[0].mapping_quality < parameters.tail_min_mapq:
        return None

    sv_evidences = []

    full_ref_chr = full_bam.getrname(full_prim[0].reference_id)
    full_ref_start = full_prim[0].reference_start
    indels = find_indels_in_cigar_tuples(full_prim[0].cigartuples)
    for pos, length, typ in indels:
        if typ == "del":
            print("Deletion detected: {0}:{1}-{2} (length {3})".format(full_ref_start, pos, full_ref_start + pos + length, length), file=sys.stdout)
            sv_evidences.append(SVEvidence(full_ref_chr, full_ref_start + pos, full_ref_start + pos + length, typ, "cigar", full_read_name))
        elif typ == "ins":
            print("Insertion detected: {0}:{1}-{2} (length {3})".format(full_ref_chr, full_ref_start + pos, full_ref_start + pos + length, length), file=sys.stdout)
            sv_evidences.append(SVEvidence(full_ref_chr, full_ref_start + pos, full_ref_start + pos + length, typ, "cigar", full_read_name))

    good_suppl_alns = [aln for aln in full_suppl if not aln.is_unmapped and aln.mapping_quality >= parameters.tail_min_mapq]
    if len(good_suppl_alns) == 1:
        result = analyze_one_supplementary(full_prim[0], good_suppl_alns[0], full_bam)
        if result != None:
            sv_evidences.append(result)
    elif len(good_suppl_alns) == 2:
        result = analyze_two_supplementary(full_prim[0], good_suppl_alns[0], good_suppl_alns[1], full_bam)
        if result != None:
            sv_evidences.append(result)
    elif 2 < len(good_suppl_alns) < 6:
        print("Read", full_read_name.split("_")[1], "has", len(good_suppl_alns), "good supplementary segments.") 
    
    return sv_evidences