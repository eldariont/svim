from __future__ import print_function

import sys

from SVEvidence import EvidenceDeletion, EvidenceInsertion


def analyze_cigar_indel(tuples, min_length=50):
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


def analyze_alignment_indel(alignment, full_bam, full_read_name):
    sv_evidences = []
    full_ref_chr = full_bam.getrname(alignment.reference_id)
    full_ref_start = alignment.reference_start
    indels = analyze_cigar_indel(alignment.cigartuples)
    for pos, length, typ in indels:
        if typ == "del":
            #print("Deletion detected: {0}:{1}-{2} (length {3})".format(full_ref_start, pos, full_ref_start + pos + length, length), file=sys.stdout)
            sv_evidences.append(EvidenceDeletion(full_ref_chr, full_ref_start + pos, full_ref_start + pos + length, "cigar", full_read_name))
        elif typ == "ins":
            #print("Insertion detected: {0}:{1}-{2} (length {3})".format(full_ref_chr, full_ref_start + pos, full_ref_start + pos + length, length), file=sys.stdout)
            sv_evidences.append(EvidenceInsertion(full_ref_chr, full_ref_start + pos, full_ref_start + pos + length, "cigar", full_read_name))
    return sv_evidences




def analyze_full_read_indel(full_iterator_object, full_bam, parameters):
    """Search the (primary and supplementary) alignments of a full read for contained insertions and deletions"""
    full_read_name, full_prim, full_suppl, full_sec = full_iterator_object

    sv_evidences = []

    if len(full_prim) != 1 or full_prim[0].is_unmapped or full_prim[0].mapping_quality < parameters.min_mapq:
        return sv_evidences

    # Search indels in primary alignment
    sv_evidences.extend(analyze_alignment_indel(full_prim[0], full_bam, full_read_name))

    # Search indels in good supplementary alignments
    good_suppl_alns = [aln for aln in full_suppl if not aln.is_unmapped and aln.mapping_quality >= parameters.min_mapq]
    for alignment in good_suppl_alns:
        sv_evidences.extend(analyze_alignment_indel(alignment, full_bam, full_read_name))    

    return sv_evidences


