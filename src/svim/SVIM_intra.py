from __future__ import print_function

import sys

from svim.SVSignature import SignatureDeletion, SignatureInsertion, SignatureTranslocation


def analyze_cigar_indel(tuples, min_length):
    """Parses CIGAR tuples (op, len) and returns Indels with a length > minLength"""
    pos_ref = 0
    pos_read = 0
    indels = []
    for operation, length in tuples:
        if operation == 0:                     # alignment match
            pos_ref += length
            pos_read += length
        elif operation == 1:                   # insertion
            if length >= min_length:
                indels.append((pos_ref, pos_read, length, "INS"))
            pos_read += length
        elif operation == 2:                   # deletion
            if length >= min_length:
                indels.append((pos_ref, pos_read, length, "DEL"))
            pos_ref += length
        elif operation == 4:                   # soft clip
            pos_read += length
        elif operation == 7 or operation == 8:        # match or mismatch
            pos_ref += length
            pos_read += length
    return indels


def analyze_alignment_indel(alignment, bam, query_name, options):
    sv_signatures = []
    #Translocation signatures from other SV classes are stored separately for --all_bnd option
    translocation_signatures_all_bnds = []
    ref_chr = bam.getrname(alignment.reference_id)
    ref_start = alignment.reference_start
    indels = analyze_cigar_indel(alignment.cigartuples, options.min_sv_size)
    for pos_ref, pos_read, length, typ in indels:
        if typ == "DEL":
            sv_signatures.append(SignatureDeletion(ref_chr, ref_start + pos_ref, ref_start + pos_ref + length, "cigar", query_name))
            if options.all_bnds:
                translocation_signatures_all_bnds.append(SignatureTranslocation(ref_chr, ref_start + pos_ref, 'fwd', ref_chr, ref_start + pos_ref + length, 'fwd', "cigar", query_name))
        elif typ == "INS":
            try:
                insertion_seq = alignment.query_sequence[pos_read:pos_read+length]
            except TypeError:
                insertion_seq = ""
            sv_signatures.append(SignatureInsertion(ref_chr, ref_start + pos_ref, ref_start + pos_ref + length, "cigar", query_name, insertion_seq))
    return sv_signatures, translocation_signatures_all_bnds


