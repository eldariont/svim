import logging
import pysam

from svim.SVIM_intra import analyze_alignment_indel
from svim.SVIM_inter import analyze_read_segments


def bam_iterator(bam):
    """Returns an iterator for the given SAM/BAM file (must be query-sorted).
    In each call, the alignments of a single read are yielded as a 3-tuple: (list of primary pysam.AlignedSegment, list of supplementary pysam.AlignedSegment, list of secondary pysam.AlignedSegment)."""
    alignments = bam.fetch(until_eof=True)
    current_aln = next(alignments)
    current_read_name = current_aln.query_name
    current_prim = []
    current_suppl = []
    current_sec = []
    if current_aln.is_secondary:
        current_sec.append(current_aln)
    elif current_aln.is_supplementary:
        current_suppl.append(current_aln)
    else:
        current_prim.append(current_aln)
    while True:
        try:
            next_aln = next(alignments)
            next_read_name = next_aln.query_name
            if next_read_name != current_read_name:
                yield (current_prim, current_suppl, current_sec)
                current_read_name = next_read_name
                current_prim = []
                current_suppl = []
                current_sec = []
            if next_aln.is_secondary:
                current_sec.append(next_aln)
            elif next_aln.is_supplementary:
                current_suppl.append(next_aln)
            else:
                current_prim.append(next_aln)
        except StopIteration:
            break
    yield (current_prim, current_suppl, current_sec)


def retrieve_supplementary_alignments(primary_alignment, bam):
    """Reconstruct supplementary alignments for a given primary alignment from the SA tag"""
    try:
        sa_tag = primary_alignment.get_tag("SA").split(";")           
    except KeyError:
        return []
    supplementary_alignments = []
    # For each supplementary alignment encoded in the SA tag
    for element in sa_tag:
        # Read information from the tag
        fields = element.split(",")
        if len(fields) != 6:
            continue
        rname = fields[0]
        pos = int(fields[1])
        strand = fields[2]
        # CIGAR string encoded in SA tag is shortened
        cigar = fields[3]
        mapq = int(fields[4])
        nm = int(fields[5])

        # Generate an aligned segment from the information
        a = pysam.AlignedSegment()
        a.query_name = primary_alignment.query_name
        a.query_sequence= primary_alignment.query_sequence
        if strand == "+":
            a.flag = 2048
        else:
            a.flag = 2064
        a.reference_id = bam.get_tid(rname)
        a.reference_start = pos - 1
        try:
            a.mapping_quality = mapq
        except OverflowError:
            a.mapping_quality = 0
        a.cigarstring = cigar
        a.next_reference_id = -1
        a.next_reference_start = -1
        a.template_length = 0
        a.query_qualities = primary_alignment.query_qualities
        a.set_tags([("NM", nm, "i")])

        supplementary_alignments.append(a)
    return supplementary_alignments


def analyze_alignment_file_querysorted(bam, options):
    alignment_it = bam_iterator(bam)

    sv_signatures = []
    read_nr = 0

    while True:
        try:
            alignment_iterator_object = next(alignment_it)
            primary_aln, suppl_aln, sec_aln = alignment_iterator_object
            if len(primary_aln) != 1 or primary_aln[0].is_unmapped or primary_aln[0].mapping_quality < options.min_mapq:
                continue
            read_nr += 1
            if read_nr % 10000 == 0:
                logging.info("Processed read {0}".format(read_nr))
            good_suppl_alns = [aln for aln in suppl_aln if not aln.is_unmapped and aln.mapping_quality >= options.min_mapq]
            if not options.skip_indel:
                sv_signatures.extend(analyze_alignment_indel(primary_aln[0], bam, primary_aln[0].query_name, options))
                for alignment in good_suppl_alns:
                    sv_signatures.extend(analyze_alignment_indel(alignment, bam, alignment.query_name, options))
            if not options.skip_segment:
                sv_signatures.extend(analyze_read_segments(primary_aln[0], good_suppl_alns, bam, options))
        except StopIteration:
            break
        except KeyboardInterrupt:
            logging.warning('Execution interrupted by user. Stop detection and continue with next step..')
            break
    return sv_signatures


def analyze_alignment_file_coordsorted(bam, options):
    alignment_it = bam.fetch(until_eof=True, multiple_iterators=True)

    sv_signatures = []
    read_nr = 0

    while True:
        try:
            current_alignment = next(alignment_it)
            if current_alignment.is_unmapped or current_alignment.is_supplementary or current_alignment.is_secondary or current_alignment.mapping_quality < options.min_mapq:
                continue
            read_nr += 1
            if read_nr % 10000 == 0:
                logging.info("Processed read {0}".format(read_nr))
            supplementary_alignments = retrieve_supplementary_alignments(current_alignment, bam)
            good_suppl_alns = [aln for aln in supplementary_alignments if not aln.is_unmapped and aln.mapping_quality >= options.min_mapq]

            if not options.skip_indel:
                sv_signatures.extend(analyze_alignment_indel(current_alignment, bam, current_alignment.query_name, options))
                for alignment in good_suppl_alns:
                    sv_signatures.extend(analyze_alignment_indel(alignment, bam, alignment.query_name, options))
            if not options.skip_segment:
                sv_signatures.extend(analyze_read_segments(current_alignment, good_suppl_alns, bam, options))
        except StopIteration:
            break
        except KeyboardInterrupt:
            logging.warning('Execution interrupted by user. Stop detection and continue with next step..')
            break
    return sv_signatures