from __future__ import print_function

import sys

from svim.SVSignature import SignatureDeletion, SignatureInsertion


def analyze_cigar_indel(tuples, min_length):
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


def analyze_alignment_indel(alignment, bam, query_name, options):
    sv_signatures = []
    ref_chr = bam.getrname(alignment.reference_id)
    ref_start = alignment.reference_start
    indels = analyze_cigar_indel(alignment.cigartuples, options.min_sv_size)
    for pos, length, typ in indels:
        if typ == "del":
            sv_signatures.append(SignatureDeletion(ref_chr, ref_start + pos, ref_start + pos + length, "cigar", query_name))
        elif typ == "ins":
            sv_signatures.append(SignatureInsertion(ref_chr, ref_start + pos, ref_start + pos + length, "cigar", query_name))
    return sv_signatures


