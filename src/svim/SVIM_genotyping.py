import time
import logging

from svim.SVIM_intra import analyze_alignment_indel
from svim.SVIM_inter import analyze_read_segments
from svim.SVIM_COLLECT import retrieve_other_alignments


def span_position_distance(candidate, signature, distance_normalizer):
    if candidate.type == "INS" or candidate.type == "DUP_INT":
        c_contig, c_start, c_end = candidate.get_destination()
    else:
        c_contig, c_start, c_end = candidate.get_source()
    if signature.type == "DUP_INT":
        s_contig, s_start, s_end = signature.get_destination()
    else:
        s_contig, s_start, s_end = signature.get_source()
    #ins signatures can support dup_int candidates and vice versa
    if not (candidate.type == "INS" and signature.type == "DUP_INT") and \
       not (candidate.type == "DUP_INT" and signature.type == "INS") and \
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
    num_candidates = len(candidates)
    for nr, candidate in enumerate(candidates):
        if (nr+1) % 10000 == 0:
            logging.info("Processed {0} of {1} candidates".format(nr+1, num_candidates))
        if candidate.score < options.minimum_score:
            continue
        #Fetch alignments around variant locus
        if type == "INS" or type == "DUP_INT":
            contig, start, end = candidate.get_destination()
            #We need the insertion locus on the reference for which end is equal to start
            end = start
        else:
            contig, start, end = candidate.get_source()
        contig_length = bam.get_reference_length(contig)
        alignment_it = bam.fetch(contig=contig, start=max(0, start-1000), stop=min(contig_length, end+1000))

        reads_supporting_variant = set([sig.read for sig in candidate.members])
        #Count reads that overlap the locus and therefore support the reference
        reads_supporting_reference = set()        
        #Loop through fetched alignments
        aln_no = 0
        while aln_no < 500:
            try:
                current_alignment = next(alignment_it)
            except StopIteration:
                break
            
            if current_alignment.query_name in reads_supporting_variant:
                continue
            if current_alignment.is_unmapped or current_alignment.is_secondary or current_alignment.mapping_quality < options.min_mapq:
                continue
            aln_no += 1

            if type == "DEL" or type == "INV":
                minimum_overlap = min((end - start) / 2, 2000)
                if (current_alignment.reference_start < (end - minimum_overlap) and current_alignment.reference_end > (end + 100) or
                    current_alignment.reference_start < (start - 100) and current_alignment.reference_end > (start + minimum_overlap)):
                    reads_supporting_reference.add(current_alignment.query_name)
            if type == "INS" or type == "DUP_INT":
                if current_alignment.reference_start < (start - 100) and current_alignment.reference_end > (end + 100):
                    reads_supporting_reference.add(current_alignment.query_name)

        if (len(reads_supporting_variant) + len(reads_supporting_reference)) >= options.minimum_depth:
            candidate.support_fraction = len(reads_supporting_variant) / (len(reads_supporting_variant) + len(reads_supporting_reference))
            if candidate.support_fraction >= options.homozygous_threshold:
                candidate.genotype = "1/1"
            elif candidate.support_fraction >= options.heterozygous_threshold and candidate.support_fraction < options.homozygous_threshold:
                candidate.genotype = "0/1"
            elif candidate.support_fraction < options.heterozygous_threshold:
                candidate.genotype = "0/0"
            else:
                candidate.genotype = "./."
        elif (len(reads_supporting_variant) + len(reads_supporting_reference)) > 0:
            candidate.support_fraction = len(reads_supporting_variant) / (len(reads_supporting_variant) + len(reads_supporting_reference))
            candidate.genotype = "./."
        else:
            candidate.support_fraction = "."
            candidate.genotype = "./."
        candidate.ref_reads = len(reads_supporting_reference)
        candidate.alt_reads = len(reads_supporting_variant)