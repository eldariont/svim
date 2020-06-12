from __future__ import print_function

import sys
from statistics import mean

from svim.SVSignature import SignatureDeletion, SignatureInsertion, SignatureInversion, SignatureTranslocation, SignatureDuplicationTandem, SignatureInsertionFrom
from svim.SVIM_clustering import consolidate_clusters_bilocal, clusters_from_partitions


def is_similar(chr1, start1, end1, chr2, start2, end2):
    if chr1 == chr2 and abs(start1 - start2) < 20 and abs(end1 - end2) < 20:
        return True
    else:
        return False


def analyze_read_segments(primary, supplementaries, bam, options):
    read_name = primary.query_name
    alignments = [primary] + supplementaries
    alignment_list = []
    for alignment in alignments:
        #correct query coordinates for reversely mapped reads
        if alignment.is_reverse:
            q_start = alignment.infer_read_length() - alignment.query_alignment_end
            q_end = alignment.infer_read_length() - alignment.query_alignment_start
        else:
            q_start = alignment.query_alignment_start
            q_end = alignment.query_alignment_end

        new_alignment_dict = {  'q_start': q_start, 
                                'q_end': q_end, 
                                'ref_id': alignment.reference_id, 
                                'ref_start': alignment.reference_start, 
                                'ref_end': alignment.reference_end,
                                'is_reverse': alignment.is_reverse  }
        alignment_list.append(new_alignment_dict)

    sorted_alignment_list = sorted(alignment_list, key=lambda aln: (aln['q_start'], aln['q_end']))
    #inferred_read_length = alignments[0].infer_read_length()

    sv_signatures = []
    #Translocation signatures from other SV classes are stored separately for --all_bnd option
    translocation_signatures_all_bnds = []
    tandem_duplications = []
    translocations = []

    for index in range(len(sorted_alignment_list) - 1):
        alignment_current = sorted_alignment_list[index]
        alignment_next = sorted_alignment_list[index + 1]

        distance_on_read = alignment_next['q_start'] - alignment_current['q_end']

        #Same chromosome
        if alignment_current['ref_id'] == alignment_next['ref_id']:
            ref_chr = bam.getrname(alignment_current['ref_id'])
            #Same orientation
            if alignment_current['is_reverse'] == alignment_next['is_reverse']:
                #Compute distance on reference depending on orientation
                if alignment_current['is_reverse']:
                    distance_on_reference = alignment_current['ref_start'] - alignment_next['ref_end']
                else:
                    distance_on_reference = alignment_next['ref_start'] - alignment_current['ref_end']
                #No overlap on read
                if distance_on_read >= -options.segment_overlap_tolerance:
                    #No overlap on reference
                    if distance_on_reference >= -options.segment_overlap_tolerance:
                        deviation = distance_on_read - distance_on_reference
                        #INS candidate
                        if deviation >= options.min_sv_size:
                            #No gap on reference
                            if distance_on_reference <= options.segment_gap_tolerance:
                                if not alignment_current['is_reverse']:
                                    try:
                                        insertion_seq = primary.query_sequence[alignment_current['q_end']:alignment_current['q_end']+deviation]
                                    except TypeError:
                                        insertion_seq = ""
                                    sv_signatures.append(SignatureInsertion(ref_chr, alignment_current['ref_end'], alignment_current['ref_end'] + deviation, "suppl", read_name, insertion_seq))
                                else:
                                    try:
                                        insertion_seq = primary.query_sequence[primary.infer_read_length() - alignment_next['q_start']:primary.infer_read_length() - alignment_next['q_start'] + deviation]
                                    except TypeError:
                                        insertion_seq = ""
                                    sv_signatures.append(SignatureInsertion(ref_chr, alignment_current['ref_start'], alignment_current['ref_start'] + deviation, "suppl", read_name, insertion_seq))
                        #DEL candidate
                        elif -options.max_sv_size <= deviation <= -options.min_sv_size:
                            #No gap on read
                            if distance_on_read <= options.segment_gap_tolerance:
                                if not alignment_current['is_reverse']:
                                    sv_signatures.append(SignatureDeletion(ref_chr, alignment_current['ref_end'], alignment_current['ref_end'] - deviation, "suppl", read_name))
                                    if options.all_bnds:
                                        translocation_signatures_all_bnds.append(SignatureTranslocation(ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_current['ref_end'] - deviation, 'fwd', "suppl", read_name))
                                else:
                                    sv_signatures.append(SignatureDeletion(ref_chr, alignment_next['ref_end'], alignment_next['ref_end'] - deviation, "suppl", read_name))
                                    if options.all_bnds:
                                        translocation_signatures_all_bnds.append(SignatureTranslocation(ref_chr, alignment_next['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_end'] - deviation, 'fwd', "suppl", read_name))
                        #Either very large DEL or TRANS
                        elif deviation < -options.max_sv_size:
                            #No gap on read
                            if distance_on_read <= options.segment_gap_tolerance:
                                if not alignment_current['is_reverse']:
                                    sv_signatures.append(SignatureTranslocation(ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_start'], 'fwd', "suppl", read_name))
                                    translocations.append(('fwd', 'fwd', ref_chr, alignment_current['ref_end'] - 1, ref_chr, alignment_next['ref_start']))
                                else:
                                    sv_signatures.append(SignatureTranslocation(ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_end'] - 1, 'rev', "suppl", read_name))
                                    translocations.append(('rev', 'rev', ref_chr, alignment_current['ref_start'], ref_chr, alignment_next['ref_end'] - 1))
                    #overlap on reference
                    else:
                        #Tandem Duplication
                        if distance_on_reference <= -options.min_sv_size:
                            if not alignment_current['is_reverse']:
                                #Tandem Duplication
                                if alignment_next['ref_end'] > alignment_current['ref_start']:
                                    tandem_duplications.append((ref_chr, alignment_next['ref_start'], alignment_current['ref_end'], True, True))
                                    if options.all_bnds:
                                        translocation_signatures_all_bnds.append(SignatureTranslocation(ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_start'], 'fwd', "suppl", read_name))
                                #Large tandem duplication
                                elif distance_on_reference >= -options.max_sv_size:
                                    tandem_duplications.append((ref_chr, alignment_next['ref_start'], alignment_current['ref_end'], False, True))
                                    if options.all_bnds:
                                        translocation_signatures_all_bnds.append(SignatureTranslocation(ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_start'], 'fwd', "suppl", read_name))
                                #Either very large TANDEM or TRANS
                                else:
                                    sv_signatures.append(SignatureTranslocation(ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_start'], 'fwd', "suppl", read_name))
                                    translocations.append(('fwd', 'fwd', ref_chr, alignment_current['ref_end'] - 1, ref_chr, alignment_next['ref_start']))
                            else:
                                #Tandem Duplication
                                if alignment_next['ref_start'] < alignment_current['ref_end']:
                                    tandem_duplications.append((ref_chr, alignment_current['ref_start'], alignment_next['ref_end'], True, False))
                                    if options.all_bnds:
                                        translocation_signatures_all_bnds.append(SignatureTranslocation(ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_end'] - 1, 'rev', "suppl", read_name))
                                #Large tandem duplication
                                elif distance_on_reference >= -options.max_sv_size:
                                    tandem_duplications.append((ref_chr, alignment_current['ref_start'], alignment_next['ref_end'], False, False))
                                    if options.all_bnds:
                                        translocation_signatures_all_bnds.append(SignatureTranslocation(ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_end'] - 1, 'rev', "suppl", read_name))
                                #Either very large TANDEM or TRANS
                                else:
                                    sv_signatures.append(SignatureTranslocation(ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_end'] - 1, 'rev', "suppl", read_name))
                                    translocations.append(('rev', 'rev', ref_chr, alignment_current['ref_start'], ref_chr, alignment_next['ref_end'] - 1))
            #Different orientations
            else:
                #Normal to reverse
                if not alignment_current['is_reverse'] and alignment_next['is_reverse']:
                    if -options.segment_overlap_tolerance <= distance_on_read <= options.segment_gap_tolerance:
                        if alignment_next['ref_start'] - alignment_current['ref_end'] >= -options.segment_overlap_tolerance: # Case 1
                            #INV candidate
                            if options.min_sv_size <= alignment_next['ref_end'] - alignment_current['ref_end'] <= options.max_sv_size:
                                sv_signatures.append(SignatureInversion(ref_chr, alignment_current['ref_end'], alignment_next['ref_end'], "suppl", read_name, "left_fwd"))
                                if options.all_bnds:
                                    translocation_signatures_all_bnds.append(SignatureTranslocation(ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_end'] - 1, 'rev', "suppl", read_name))
                            #Either very large INV or TRANS
                            elif alignment_next['ref_end'] - alignment_current['ref_end'] > options.max_sv_size:
                                sv_signatures.append(SignatureTranslocation(ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_end'] - 1, 'rev', "suppl", read_name))
                                translocations.append(('fwd', 'rev', ref_chr, alignment_current['ref_end'] - 1, ref_chr, alignment_next['ref_end'] - 1))
                        elif alignment_current['ref_start'] - alignment_next['ref_end'] >= -options.segment_overlap_tolerance: # Case 3
                            #INV candidate
                            if options.min_sv_size <= alignment_current['ref_end'] - alignment_next['ref_end'] <= options.max_sv_size:
                                sv_signatures.append(SignatureInversion(ref_chr, alignment_next['ref_end'], alignment_current['ref_end'], "suppl", read_name, "left_rev"))
                                if options.all_bnds:
                                    translocation_signatures_all_bnds.append(SignatureTranslocation(ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_end'] - 1, 'rev', "suppl", read_name))
                            #Either very large INV or TRANS
                            elif alignment_current['ref_end'] - alignment_next['ref_end'] > options.max_sv_size:
                                sv_signatures.append(SignatureTranslocation(ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_end'] - 1, 'rev', "suppl", read_name))
                                translocations.append(('fwd', 'rev', ref_chr, alignment_current['ref_end'] - 1, ref_chr, alignment_next['ref_end'] - 1))
                    else:
                        pass
                        #print("Overlapping read segments in read", read_name)
                #Reverse to normal
                if alignment_current['is_reverse'] and not alignment_next['is_reverse']:
                    if -options.segment_overlap_tolerance <= distance_on_read <= options.segment_gap_tolerance:
                        if alignment_next['ref_start'] - alignment_current['ref_end'] >= -options.segment_overlap_tolerance: # Case 2
                            #INV candidate
                            if options.min_sv_size <= alignment_next['ref_start'] - alignment_current['ref_start'] <= options.max_sv_size:
                                sv_signatures.append(SignatureInversion(ref_chr, alignment_current['ref_start'], alignment_next['ref_start'], "suppl", read_name, "right_fwd"))
                                if options.all_bnds:
                                    translocation_signatures_all_bnds.append(SignatureTranslocation(ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_start'], 'fwd', "suppl", read_name))
                            #Either very large INV or TRANS
                            elif alignment_next['ref_start'] - alignment_current['ref_start'] > options.max_sv_size:
                                sv_signatures.append(SignatureTranslocation(ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_start'], 'fwd', "suppl", read_name))
                                translocations.append(('rev', 'fwd', ref_chr, alignment_current['ref_start'], ref_chr, alignment_next['ref_start']))
                        elif alignment_current['ref_start'] - alignment_next['ref_end'] >= -options.segment_overlap_tolerance: # Case 4
                            #INV candidate
                            if options.min_sv_size <= alignment_current['ref_start'] - alignment_next['ref_start'] <= options.max_sv_size:
                                sv_signatures.append(SignatureInversion(ref_chr, alignment_next['ref_start'], alignment_current['ref_start'], "suppl", read_name, "right_rev"))
                                if options.all_bnds:
                                    translocation_signatures_all_bnds.append(SignatureTranslocation(ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_start'], 'fwd', "suppl", read_name))
                            #Either very large INV or TRANS
                            elif alignment_current['ref_start'] - alignment_next['ref_start'] > options.max_sv_size:
                                sv_signatures.append(SignatureTranslocation(ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_start'], 'fwd', "suppl", read_name))
                                translocations.append(('rev', 'fwd', ref_chr, alignment_current['ref_start'], ref_chr, alignment_next['ref_start']))
                    else:
                        pass
                        #print("Overlapping read segments in read", read_name)
        #Different chromosomes
        else:
            ref_chr_current = bam.getrname(alignment_current['ref_id'])
            ref_chr_next = bam.getrname(alignment_next['ref_id'])
            #Same orientation
            if alignment_current['is_reverse'] == alignment_next['is_reverse']:
                #No overlap on read
                if distance_on_read >= -options.segment_overlap_tolerance:
                    #No gap on read
                    if distance_on_read <= options.segment_gap_tolerance:
                        if not alignment_current['is_reverse']:
                            sv_signatures.append(SignatureTranslocation(ref_chr_current, alignment_current['ref_end'] - 1, 'fwd', ref_chr_next, alignment_next['ref_start'], 'fwd', "suppl", read_name))
                            translocations.append(('fwd', 'fwd', ref_chr_current, alignment_current['ref_end'] - 1, ref_chr_next, alignment_next['ref_start']))
                        else:
                            sv_signatures.append(SignatureTranslocation(ref_chr_current, alignment_current['ref_start'], 'rev', ref_chr_next, alignment_next['ref_end'] - 1, 'rev', "suppl", read_name))
                            translocations.append(('rev', 'rev', ref_chr_current, alignment_current['ref_start'], ref_chr_next, alignment_next['ref_end'] - 1))
                #Overlap on read
                else:
                    pass
                    #print("Overlapping read segments in read", read_name)
            #Different orientation
            else:
                #No overlap on read
                if distance_on_read >= -options.segment_overlap_tolerance:
                    #No gap on read
                    if distance_on_read <= options.segment_gap_tolerance:
                        if not alignment_current['is_reverse']:
                            sv_signatures.append(SignatureTranslocation(ref_chr_current, alignment_current['ref_end'] - 1, 'fwd', ref_chr_next, alignment_next['ref_end'] - 1, 'rev', "suppl", read_name))
                            translocations.append(('fwd', 'rev', ref_chr_current, alignment_current['ref_end'] - 1, ref_chr_next, alignment_next['ref_end'] - 1))
                        else:
                            sv_signatures.append(SignatureTranslocation(ref_chr_current, alignment_current['ref_start'], 'rev', ref_chr_next, alignment_next['ref_start'], 'fwd', "suppl", read_name))
                            translocations.append(('rev', 'fwd', ref_chr_current, alignment_current['ref_start'], ref_chr_next, alignment_next['ref_start']))
                #Overlap on read
                else:
                    pass
                    #print("Overlapping read segments in read", read_name)

    #Handle tandem duplications
    current_chromosome = None
    current_starts = []
    current_ends = []
    current_copy_number = 0
    current_fully_covered = []
    for tandem_duplication in tandem_duplications:
        if current_chromosome == None:
            current_chromosome = tandem_duplication[0]
            current_starts.append(tandem_duplication[1])
            current_ends.append(tandem_duplication[2])
            current_copy_number = 1
            current_fully_covered.append(tandem_duplication[3])
            current_direction = tandem_duplication[4]
        else:
            if is_similar(current_chromosome, mean(current_starts), mean(current_ends), tandem_duplication[0], tandem_duplication[1], tandem_duplication[2]) and current_direction == tandem_duplication[4]:
                current_starts.append(tandem_duplication[1])
                current_ends.append(tandem_duplication[2])
                current_copy_number += 1
                current_fully_covered.append(tandem_duplication[3])
            else:
                fully_covered = True if sum(current_fully_covered) else False
                sv_signatures.append(SignatureDuplicationTandem(current_chromosome, int(mean(current_starts)), int(mean(current_ends)), current_copy_number, fully_covered, "suppl", read_name))
                current_chromosome = tandem_duplication[0]
                current_starts =[tandem_duplication[1]]
                current_ends =[tandem_duplication[2]]
                current_copy_number = 1
                current_fully_covered = [tandem_duplication[3]]
    if current_chromosome != None:
        fully_covered = True if sum(current_fully_covered) else False
        sv_signatures.append(SignatureDuplicationTandem(current_chromosome, int(mean(current_starts)), int(mean(current_ends)), current_copy_number, fully_covered, "suppl", read_name))

    #Handle insertions_from
    for this_index in range(len(translocations)):
        this_dir1 = translocations[this_index][0]
        this_dir2 = translocations[this_index][1]
        this_chr1 = translocations[this_index][2]
        this_pos1 = translocations[this_index][3]
        this_chr2 = translocations[this_index][4]
        this_pos2 = translocations[this_index][5]

        for before_dir1, before_dir2, before_chr1, before_pos1, before_chr2, before_pos2 in translocations[:this_index]:
            #Same direction at destination and origin
            if before_dir1 == this_dir2 and before_dir2 == this_dir1:
                #Same position at destination
                if is_similar(before_chr1, before_pos1, 0, this_chr2, this_pos2, 0):
                    #Same chromosome for origin
                    if before_chr2 == this_chr1:
                        #INS_DUP candidate
                        if before_dir2 == before_dir1:
                            if before_dir1 == 'fwd':
                                if options.min_sv_size <= this_pos1 - before_pos2 + 1 <= options.max_sv_size:
                                    sv_signatures.append(SignatureInsertionFrom(before_chr2, before_pos2, this_pos1 + 1, before_chr1, int(mean([before_pos1 + 1, this_pos2])), "suppl", read_name))
                            elif before_dir1 == 'rev':
                                if options.min_sv_size <= before_pos2 - this_pos1 <= options.max_sv_size:
                                    sv_signatures.append(SignatureInsertionFrom(before_chr2, this_pos1, before_pos2 + 1, before_chr1, int(mean([before_pos1, this_pos2 + 1])), "suppl", read_name))
                        #INV_INS_DUP candidate
                        else:
                            pass

    return sv_signatures, translocation_signatures_all_bnds
