from svim.SVIM_intra import analyze_alignment_indel
from svim.SVIM_inter import analyze_read_segments
from svim.SVIM_COLLECT import retrieve_other_alignments


def span_position_distance(candidate, signature, distance_normalizer):
    if candidate.type == "ins" or candidate.type == "dup_int":
        c_contig, c_start, c_end = candidate.get_destination()
    else:
        c_contig, c_start, c_end = candidate.get_source()
    if signature.type == "dup_int":
        s_contig, s_start, s_end = signature.get_destination()
    else:
        s_contig, s_start, s_end = signature.get_source()
    #ins signatures can support dup_int candidates and vice versa
    if not (candidate.type == "ins" and signature.type == "dup_int") and \
       not (candidate.type == "dup_int" and signature.type == "ins") and \
       candidate.type != signature.type:
        return float("inf")
    if c_contig != s_contig:
        return float("inf")
    span1 = c_end - c_start
    span2 = s_end - s_start
    center1 = (c_start + c_end) // 2
    center2 = (s_start + s_end) // 2
    position_distance = min(abs(c_start - s_start), abs(c_end - s_end), abs(center1 - center2)) / distance_normalizer
    span_distance = abs(span1 - span2) / max(span1, span2)
    return position_distance + span_distance


def genotype(candidates, bam, type, options):
    for candidate in candidates:
        #Fetch alignments around variant locus
        if type == "ins" or type == "dup_int":
            contig, start, end = candidate.get_destination()
            #We need the insertion locus on the reference for which end is equal to start
            end = start
        else:
            contig, start, end = candidate.get_source()
        alignment_it = bam.fetch(contig=contig, start=start-1000, stop=end+1000)

        sv_signatures = []
        reads_covering_variant = set()
        #Loop through fetched alignments
        while True:
            try:
                current_alignment = next(alignment_it)
                if current_alignment.is_unmapped or current_alignment.is_secondary or current_alignment.mapping_quality < options.min_mapq:
                    continue
                other_alignments = retrieve_other_alignments(current_alignment, bam)
                good_other_alns = [aln for aln in other_alignments if not aln.is_unmapped and aln.mapping_quality >= options.min_mapq]

                if type == "del" or type == "inv" or type == "ins" or type == "dup_int":
                    if current_alignment.reference_start < start - 100 and current_alignment.reference_end > end + 100:
                        if not type == "inv" and not options.skip_indel:
                            sv_signatures.extend(analyze_alignment_indel(current_alignment, bam, current_alignment.query_name, options))
                        reads_covering_variant.add(current_alignment.query_name)
                    if not options.skip_segment:
                        sv_signatures.extend(analyze_read_segments(current_alignment, good_other_alns, bam, options))
            except StopIteration:
                break
        reads_supporting_variant = set([sig.read for sig in sv_signatures if span_position_distance(candidate, sig, options.distance_normalizer) < options.cluster_max_distance])
        reads_supporting_reference = reads_covering_variant - reads_supporting_variant
        if (len(reads_supporting_variant) + len(reads_supporting_reference)) > 3:
            candidate.support_fraction = len(reads_supporting_variant) / (len(reads_supporting_variant) + len(reads_supporting_reference))
            if candidate.support_fraction > 0.8:
                candidate.genotype = 2
            elif candidate.support_fraction > 0.3 and candidate.support_fraction < 0.7:
                candidate.genotype = 1
            elif candidate.support_fraction < 0.2:
                candidate.genotype = 0
            else:
                candidate.genotype = None
        elif (len(reads_supporting_variant) + len(reads_supporting_reference)) > 0:
            candidate.support_fraction = len(reads_supporting_variant) / (len(reads_supporting_variant) + len(reads_supporting_reference))
            candidate.genotype = None
        else:
            candidate.support_fraction = None
            candidate.genotype = None
        candidate.ref_reads = len(reads_supporting_reference)
        candidate.alt_reads = len(reads_supporting_variant)