import os
import logging
import re

from collections import defaultdict
from math import pow, sqrt
import time
from statistics import mean, stdev
from pysam import FastaFile
from spoa import poa
from cpuinfo import get_cpu_info

from svim.SVIM_clustering import form_partitions, partition_and_cluster_candidates, calculate_score
from svim.SVCandidate import CandidateInversion, CandidateDuplicationTandem, CandidateDeletion, CandidateNovelInsertion, CandidateBreakend
from svim.SVIM_merging import flag_cutpaste_candidates, merge_translocations_at_insertions


def write_candidates(working_dir, candidates):
    int_duplication_candidates, inversion_candidates, tan_duplication_candidates, deletion_candidates, novel_insertion_candidates, breakend_candidates = candidates

    if not os.path.exists(working_dir + '/candidates'):
        os.mkdir(working_dir + '/candidates')
    deletion_candidate_output = open(working_dir + '/candidates/candidates_deletions.bed', 'w')
    inversion_candidate_output = open(working_dir + '/candidates/candidates_inversions.bed', 'w')
    tandem_duplication_candidate_source_output = open(working_dir + '/candidates/candidates_tan_duplications_source.bed', 'w')
    tandem_duplication_candidate_dest_output = open(working_dir + '/candidates/candidates_tan_duplications_dest.bed', 'w')
    interspersed_duplication_candidate_source_output = open(working_dir + '/candidates/candidates_int_duplications_source.bed', 'w')
    interspersed_duplication_candidate_dest_output = open(working_dir + '/candidates/candidates_int_duplications_dest.bed', 'w')
    novel_insertion_candidate_output = open(working_dir + '/candidates/candidates_novel_insertions.bed', 'w')
    breakend_candidate_output = open(working_dir + '/candidates/candidates_breakends.bed', 'w')

    for candidate in deletion_candidates:
        print(candidate.get_bed_entry(), file=deletion_candidate_output)
    for candidate in int_duplication_candidates:
        bed_entries = candidate.get_bed_entries()
        print(bed_entries[0], file=interspersed_duplication_candidate_source_output)
        print(bed_entries[1], file=interspersed_duplication_candidate_dest_output)
    for candidate in inversion_candidates:
        print(candidate.get_bed_entry(), file=inversion_candidate_output)
    for candidate in tan_duplication_candidates:
        bed_entries = candidate.get_bed_entries()
        print(bed_entries[0], file=tandem_duplication_candidate_source_output)
        print(bed_entries[1], file=tandem_duplication_candidate_dest_output)
    for candidate in novel_insertion_candidates:
        print(candidate.get_bed_entry(), file=novel_insertion_candidate_output)
    for candidate in breakend_candidates:
        bed_entries = candidate.get_bed_entries()
        print(bed_entries[0], file=breakend_candidate_output)
        print(bed_entries[1], file=breakend_candidate_output)

    deletion_candidate_output.close()
    inversion_candidate_output.close()
    interspersed_duplication_candidate_source_output.close()
    interspersed_duplication_candidate_dest_output.close()
    tandem_duplication_candidate_source_output.close()
    tandem_duplication_candidate_dest_output.close()
    novel_insertion_candidate_output.close()
    breakend_candidate_output.close()


def sorted_nicely(vcf_entries):
    """ Sort the given vcf entries (in the form ((contig, start, end), vcf_string, sv_type)) in the way that humans expect.
        e.g. chr10 comes after chr2
        Algorithm adapted from https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/"""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    tuple_key = lambda entry: ( alphanum_key(str(entry[0][0])), entry[0][1], entry[0][2] )
    return sorted(vcf_entries, key = tuple_key)


def write_final_vcf(int_duplication_candidates,
                    inversion_candidates, 
                    tandem_duplication_candidates, 
                    deletion_candidates, 
                    novel_insertion_candidates, 
                    breakend_candidates, 
                    version, 
                    contig_names, 
                    contig_lengths,
                    types_to_output,
                    options):
    vcf_output = open(options.working_dir + '/variants.vcf', 'w')

    # Write header lines
    print("##fileformat=VCFv4.2", file=vcf_output)
    print("##fileDate={0}".format(time.strftime("%Y-%m-%d|%I:%M:%S%p|%Z|%z")), file=vcf_output)
    print("##source=SVIM-v{0}".format(version), file=vcf_output)
    for contig_name, contig_length in zip(contig_names, contig_lengths):
        print("##contig=<ID={0},length={1}>".format(contig_name, contig_length), file=vcf_output)
    if "DEL" in types_to_output:
        print("##ALT=<ID=DEL,Description=\"Deletion\">", file=vcf_output)
    if "INV" in types_to_output:
        print("##ALT=<ID=INV,Description=\"Inversion\">", file=vcf_output)
    if (not options.tandem_duplications_as_insertions and "DUP:TANDEM" in types_to_output) or \
       (not options.interspersed_duplications_as_insertions and "DUP:INT" in types_to_output):
        print("##ALT=<ID=DUP,Description=\"Duplication\">", file=vcf_output)
    if not options.tandem_duplications_as_insertions and "DUP:TANDEM" in types_to_output:
        print("##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">", file=vcf_output)
    if not options.interspersed_duplications_as_insertions and "DUP:INT" in types_to_output:
        print("##ALT=<ID=DUP:INT,Description=\"Interspersed Duplication\">", file=vcf_output)
    if "INS" in types_to_output:
        print("##ALT=<ID=INS,Description=\"Insertion\">", file=vcf_output)
    if "BND" in types_to_output:
        print("##ALT=<ID=BND,Description=\"Breakend\">", file=vcf_output)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">", file=vcf_output)
    print("##INFO=<ID=CUTPASTE,Number=0,Type=Flag,Description=\"Genomic origin of interspersed duplication seems to be deleted\">", file=vcf_output)
    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">", file=vcf_output)
    print("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">", file=vcf_output)
    print("##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Number of reads supporting this variant\">", file=vcf_output)
    print("##INFO=<ID=STD_SPAN,Number=1,Type=Float,Description=\"Standard deviation in span of merged SV signatures\">", file=vcf_output)
    print("##INFO=<ID=STD_POS,Number=1,Type=Float,Description=\"Standard deviation in position of merged SV signatures\">", file=vcf_output)
    print("##INFO=<ID=STD_POS1,Number=1,Type=Float,Description=\"Standard deviation of breakend 1 position\">", file=vcf_output)
    print("##INFO=<ID=STD_POS2,Number=1,Type=Float,Description=\"Standard deviation of breakend 2 position\">", file=vcf_output)
    if options.insertion_sequences:
        print("##INFO=<ID=SEQS,Number=.,Type=String,Description=\"Insertion sequences from all supporting reads\">", file=vcf_output)
    if options.read_names:
        print("##INFO=<ID=READS,Number=.,Type=String,Description=\"Names of all supporting reads\">", file=vcf_output)
    if options.zmws:
        print("##INFO=<ID=ZMWS,Number=1,Type=Integer,Description=\"Number of supporting ZMWs (PacBio only)\">", file=vcf_output)
    print("##FILTER=<ID=hom_ref,Description=\"Genotype is homozygous reference\">", file=vcf_output)
    print("##FILTER=<ID=not_fully_covered,Description=\"Tandem duplication is not fully covered by a single read\">", file=vcf_output)
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file=vcf_output)
    print("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">", file=vcf_output)
    print("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth for each allele\">", file=vcf_output)
    if not options.tandem_duplications_as_insertions and "DUP:TANDEM" in types_to_output:
        print("##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of tandem duplication (e.g. 2 for one additional copy)\">", file=vcf_output)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + options.sample, file=vcf_output)

    # Open reference genome sequence file
    sequence_alleles = not options.symbolic_alleles
    if sequence_alleles:
        try:
            reference = FastaFile(options.genome)
        except ValueError:
            logging.warning("The given reference genome is missing an index file ({path}.fai). Sequence alleles cannot be retrieved.".format(options.genome))
            sequence_alleles = False
        except IOError:
            logging.warning("The given reference genome is missing ({path}). Sequence alleles cannot be retrieved.".format(options.genome))
            sequence_alleles = False
    else:
        reference = None

    # Prepare VCF entries depending on command-line parameters
    vcf_entries = []
    if "DEL" in types_to_output:
        for candidate in deletion_candidates:
            vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry(sequence_alleles, reference, options.read_names, options.zmws), "DEL"))
    if "INV" in types_to_output:
        for candidate in inversion_candidates:
            vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry(sequence_alleles, reference, options.read_names, options.zmws), "INV"))
    if "INS" in types_to_output:
        for candidate in novel_insertion_candidates:
            vcf_entries.append((candidate.get_destination(), candidate.get_vcf_entry(sequence_alleles, reference, options.insertion_sequences, options.read_names, options.zmws), "INS"))
    if options.tandem_duplications_as_insertions:
        if "INS" in types_to_output:
            for candidate in tandem_duplication_candidates:
                vcf_entries.append((candidate.get_destination(), candidate.get_vcf_entry_as_ins(sequence_alleles, reference, options.read_names, options.zmws), "INS"))
    else:
        if "DUP:TANDEM" in types_to_output:
            for candidate in tandem_duplication_candidates:
                vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry_as_dup(options.read_names, options.zmws), "DUP_TANDEM"))
    if options.interspersed_duplications_as_insertions:
        if "INS" in types_to_output:
            for candidate in int_duplication_candidates:
                vcf_entries.append((candidate.get_destination(), candidate.get_vcf_entry_as_ins(sequence_alleles, reference, options.read_names, options.zmws), "INS"))
    else:
        if "DUP:INT" in types_to_output:
            for candidate in int_duplication_candidates:
                vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry_as_dup(options.read_names, options.zmws), "DUP_INT"))
    if "BND" in types_to_output:
        for candidate in breakend_candidates:
            vcf_entries.append(((candidate.get_source()[0], candidate.get_source()[1], candidate.get_source()[1] + 1), candidate.get_vcf_entry(options.read_names, options.zmws), "BND"))
            vcf_entries.append(((candidate.get_destination()[0], candidate.get_destination()[1], candidate.get_destination()[1] + 1), candidate.get_vcf_entry_reverse(options.read_names, options.zmws), "BND"))

    if sequence_alleles:
        reference.close()

    # Sort and write entries to VCF
    svtype_counter = defaultdict(int)
    for source, entry, svtype in sorted_nicely(vcf_entries):
        variant_id = "svim.{svtype}.{number}".format(svtype = svtype, number = svtype_counter[svtype] + 1)
        entry_with_id = entry.replace("PLACEHOLDERFORID", variant_id, 1)
        svtype_counter[svtype] += 1
        print(entry_with_id, file=vcf_output)

    vcf_output.close()

def generate_insertion_consensus(ins_cluster, reference, window_padding = 100, maximum_haplotype_length = 10000, allowed_size_deviation = 2.0):
    #compute window containing all members of cluster
    member_pos = [member.start for member in ins_cluster.members]
    window_start = min(member_pos) - window_padding
    window_end = max(member_pos) + window_padding

    #construct haplotype sequences from all reads
    haplotypes = []
    for member in ins_cluster.members:
        haplotype_sequence = reference.fetch(ins_cluster.contig, max(0, window_start), max(0, member.start)).upper()
        haplotype_sequence += member.sequence.upper()
        haplotype_sequence += reference.fetch(ins_cluster.contig, max(0, member.start), max(0, window_end)).upper()
        haplotypes.append(haplotype_sequence)
    largest_haplotype_length = max([len(h) for h in haplotypes])
    if largest_haplotype_length > maximum_haplotype_length:
        logging.info("Skipping consensus computation for insertion with haplotypes exceeding maximum length ({0} > {1})".format(largest_haplotype_length, maximum_haplotype_length))
        return (1, ())

    #compute consensus of all reads using SPOA
    try:
        consensus_reads, msa_reads = poa(haplotypes, algorithm=1, m=2, n=-4, g=-4, e=-2, q=-24, c=-1)
    except MemoryError:
        logging.warning("Error: spoa ran out of memory while computing a consensus sequence for a cluster of insertion signatures (size = {0}, maximum haplotype length = {1}).".format(ins_cluster.size, largest_haplotype_length))
        return (2, ())
    except ValueError as ve:
        logging.warning("Error: spoa failed while computing a consensus sequence for a cluster of insertion signatures (size = {0}, maximum haplotype length = {1}): {2}".format(ins_cluster.size, largest_haplotype_length, ve))
        return (2, ())

    #re-align consensus sequence to reference sequence in the window
    ref_sequence = reference.fetch(ins_cluster.contig, max(0, window_start), max(0, window_end)).upper()
    try:
        consensus_reads_ref, msa_reads_ref = poa([consensus_reads, ref_sequence], algorithm=1, m=2, n=-4, g=-4, e=-2, q=-24, c=-1)
    except MemoryError:
        logging.warning("Error: spoa ran out of memory while aligning the consensus sequence for a cluster of insertion signatures (size = {0}, maximum haplotype length = {1}).".format(ins_cluster.size, largest_haplotype_length))
        return (2, ())
    except ValueError as ve:
        logging.warning("Error: spoa failed while aligning the consensus sequence for a cluster of insertion signatures (size = {0}, maximum haplotype length = {1}): {2}".format(ins_cluster.size, largest_haplotype_length, ve))
        return (2, ())

    #locate insertion relative to reference and check whether size is close to expected size
    expected_size = ins_cluster.end - ins_cluster.start
    matches = []
    for match in re.finditer(r'-+', msa_reads_ref[1]):
        match_size = match.end() - match.start()
        size_ratio = max(match_size, expected_size) / min(match_size, expected_size)
        matches.append((match.start(), match_size, size_ratio))
    good_matches = [m for m in matches if m[2] < allowed_size_deviation]
    if len(good_matches) == 0:
        logging.info("Consensus failure (no suitable insertion found in realignment step). \
                      Expected size: {0}; Match sizes: {1}; Ref: {2}; Consensus: {3}".format(expected_size,
                                                                                             "/".join([str(m[1]) for m in matches]),
                                                                                             msa_reads_ref[1],
                                                                                             msa_reads_ref[0]))
        return (3, ())
    elif len(good_matches) == 1:
        realigned_insertion_start = max(0, window_start) + good_matches[0][0]
        realigned_insertion_size = good_matches[0][1]
        insertion_consensus = msa_reads_ref[0][good_matches[0][0]:good_matches[0][0]+good_matches[0][1]]
        logging.debug("Consensus success. Expected size: {0}; Consensus size: {1}; Consensus insertion: {2}".format(expected_size, realigned_insertion_size, insertion_consensus))
        return (0, (realigned_insertion_start, realigned_insertion_size, insertion_consensus))
    else:
        logging.info("Consensus failure (multiple suitable insertions found in realignment step). \
                      Expected size: {0}; Match sizes: {1}; Ref: {2}; Consensus: {3}".format(expected_size,
                                                                                             "/".join([str(m[1]) for m in matches]),
                                                                                             msa_reads_ref[1],
                                                                                             msa_reads_ref[0]))
        return (4, ())


def prepare_insertion_candidates(insertion_signature_clusters, options):
    novel_insertion_candidates = []
    if options.skip_consensus or ("sse4_1" not in get_cpu_info()["flags"]):
        if options.skip_consensus:
            logging.info("Skipping computation of insertion consensus sequences because of --skip_consensus flag.")
        else:
            logging.warning("Skipping computation of insertion consensus sequences "
                            "because CPU does not support SSE 4.1 instruction set.")
        for ins_cluster in insertion_signature_clusters:
            if ins_cluster.score > 0:
                novel_insertion_candidates.append(CandidateNovelInsertion(ins_cluster.contig,
                                                                              ins_cluster.start,
                                                                              ins_cluster.end,
                                                                              "",
                                                                              ins_cluster.members,
                                                                              ins_cluster.score,
                                                                              ins_cluster.std_span,
                                                                              ins_cluster.std_pos))
        return novel_insertion_candidates
    logging.info("Generating and realigning consensus sequence for insertions..")
    with FastaFile(options.genome) as reference:
        # 0-successful, 1-skipped, 2-failed, 3-no consensus, 4-multiple consensuses
        status_counter = [0, 0, 0, 0, 0]
        for ins_cluster in insertion_signature_clusters:
            if ins_cluster.score > 0:
                if len(ins_cluster.members) < 3:
                    novel_insertion_candidates.append(CandidateNovelInsertion(ins_cluster.contig,
                                                                              ins_cluster.start,
                                                                              ins_cluster.end,
                                                                              ins_cluster.members[0].sequence,
                                                                              ins_cluster.members,
                                                                              ins_cluster.score,
                                                                              ins_cluster.std_span,
                                                                              ins_cluster.std_pos))
                    continue
                status, consensus_result = generate_insertion_consensus(ins_cluster,
                                                                        reference,
                                                                        maximum_haplotype_length=options.max_consensus_length)
                try:
                    status_counter[status] += 1
                except KeyError:
                    logging.error("Unknown status {0} returned by consensus computation.".format(status))
                    continue
                # Successful
                if status == 0:
                    realigned_insertion_start, realigned_insertion_size, insertion_consensus = consensus_result
                    novel_insertion_candidates.append(CandidateNovelInsertion(ins_cluster.contig,
                                                                              realigned_insertion_start,
                                                                              realigned_insertion_start + realigned_insertion_size,
                                                                              insertion_consensus,
                                                                              ins_cluster.members,
                                                                              ins_cluster.score,
                                                                              ins_cluster.std_span,
                                                                              ins_cluster.std_pos))
                #Unsuccessful
                else:
                    novel_insertion_candidates.append(CandidateNovelInsertion(ins_cluster.contig,
                                                                              ins_cluster.start,
                                                                              ins_cluster.end,
                                                                              "",
                                                                              ins_cluster.members,
                                                                              ins_cluster.score,
                                                                              ins_cluster.std_span,
                                                                              ins_cluster.std_pos))
    message = ("Generated and realigned consensus sequences for {0} insertions "
               "({1} skipped, {2} failed with an error, {3} failed with no "
               "consensus, {4} failed with multiple consensuses).")
    logging.info(message.format(status_counter[0],
                                status_counter[1],
                                status_counter[2],
                                status_counter[3],
                                status_counter[4]))
    return novel_insertion_candidates


def combine_clusters(signature_clusters, options):
    deletion_signature_clusters, insertion_signature_clusters, inversion_signature_clusters, tandem_duplication_signature_clusters, insertion_from_signature_clusters, translocation_signature_clusters = signature_clusters

    ###############################
    # Create inversion candidates #
    ###############################
    inversion_candidates = []
    for inv_cluster in inversion_signature_clusters:
        inversion_candidates.append(CandidateInversion(inv_cluster.contig, inv_cluster.start, inv_cluster.end, inv_cluster.members, inv_cluster.score, inv_cluster.std_span, inv_cluster.std_pos))

    ########################################
    # Create tandem duplication candidates #
    ########################################
    tan_dup_candidates = []
    for tan_dup_cluster in tandem_duplication_signature_clusters:
        source_contig, source_start, source_end = tan_dup_cluster.get_source()
        dest_contig, dest_start, dest_end = tan_dup_cluster.get_destination()
        num_copies = int(round((dest_end - dest_start) / (source_end - source_start)))
        fully_covered = True if sum([sig.fully_covered for sig in tan_dup_cluster.members]) else False
        tan_dup_candidates.append(CandidateDuplicationTandem(tan_dup_cluster.source_contig, tan_dup_cluster.source_start, tan_dup_cluster.source_end, num_copies, fully_covered, tan_dup_cluster.members, tan_dup_cluster.score, tan_dup_cluster.std_span, tan_dup_cluster.std_pos))

    #####################################
    # Cluster translocation breakpoints #
    #####################################

    # Cluster translocations by contig and pos1
    # logging.info("Cluster translocation breakpoints..")
    # translocations_fwdfwd = [tra for tra in translocation_signature_clusters if tra.direction1 == "fwd" and tra.direction2 == "fwd"]
    # translocations_revrev = [tra for tra in translocation_signature_clusters if tra.direction1 == "rev" and tra.direction2 == "rev"]
    # translocations_fwdrev = [tra for tra in translocation_signature_clusters if tra.direction1 == "fwd" and tra.direction2 == "rev"]
    # translocations_revfwd = [tra for tra in translocation_signature_clusters if tra.direction1 == "rev" and tra.direction2 == "fwd"]
    # translocation_partitions_fwdfwd = form_partitions(translocations_fwdfwd, options.trans_partition_max_distance)
    # translocation_partitions_revrev = form_partitions(translocations_revrev, options.trans_partition_max_distance)
    # translocation_partitions_fwdrev = form_partitions(translocations_fwdrev, options.trans_partition_max_distance)
    # translocation_partitions_revfwd = form_partitions(translocations_revfwd, options.trans_partition_max_distance)

    ##############################
    # Create breakend candidates #
    ##############################

    breakend_candidates = []
    for tra_cluster in translocation_signature_clusters:
        breakend_candidates.append(CandidateBreakend(tra_cluster.source_contig, 
                                                     tra_cluster.source_start, 
                                                     tra_cluster.direction1, 
                                                     tra_cluster.dest_contig, 
                                                     tra_cluster.dest_start, 
                                                     tra_cluster.direction2, 
                                                     tra_cluster.members, 
                                                     tra_cluster.score, 
                                                     tra_cluster.std_span, 
                                                     tra_cluster.std_pos))

    ###################################################
    # Merge translocation breakpoints with insertions #
    ###################################################

    logging.info("Combine inserted regions with translocation breakpoints..")
    new_insertion_from_clusters, inserted_regions_to_remove_1 = merge_translocations_at_insertions(translocation_signature_clusters, insertion_signature_clusters, options)
    insertion_from_signature_clusters.extend(new_insertion_from_clusters)

    ############################################################################
    # Create interspersed duplication candidates and flag cut&paste insertions #
    ############################################################################

    logging.info("Create interspersed duplication candidates and flag cut&paste insertions..")
    int_duplication_candidates = flag_cutpaste_candidates(insertion_from_signature_clusters, deletion_signature_clusters, options)

    ###################################
    # Remove inserted region clusters #
    ###################################

    #find all inserted regions overlapping interspersed duplication or tandem duplication candidates
    int_duplication_iterator = iter(sorted(int_duplication_candidates, key=lambda cand: cand.get_destination()))
    tan_duplication_iterator = iter(sorted(tan_dup_candidates, key=lambda cand: cand.get_destination()))
    int_duplications_end = False
    tan_duplications_end = False
    inserted_regions_to_remove_2 = []

    try:
        current_int_duplication = next(int_duplication_iterator)
    except StopIteration:
        int_duplications_end = True

    try:
        current_tan_duplication = next(tan_duplication_iterator)
    except StopIteration:
        tan_duplications_end = True

    for inserted_region_index, inserted_region in enumerate(insertion_signature_clusters):
        contig1, start1, end1 = inserted_region.get_source()
        length1 = end1 - start1
        if not int_duplications_end:
            contig2, start2, end2 = current_int_duplication.get_destination()
            while contig2 < contig1 or (contig2 == contig1 and end2 < start1):
                try:
                    current_int_duplication = next(int_duplication_iterator)
                    contig2, start2, end2 = current_int_duplication.get_destination()
                except StopIteration:
                    int_duplications_end = True
                    break
        if not int_duplications_end:
            length2 = end2 - start2
            #if overlapping interspersed duplication of similar length
            if contig2 == contig1 and start2 < end1 and (length1 - length2) / max(length1, length2) < 0.2:
                inserted_regions_to_remove_2.append(inserted_region_index)
        else:
            if not tan_duplications_end:
                contig2, start2, end2 = current_tan_duplication.get_destination()
                while contig2 < contig1 or (contig2 == contig1 and end2 < start1):
                    try:
                        current_tan_duplication = next(tan_duplication_iterator)
                        contig2, start2, end2 = current_tan_duplication.get_destination()
                    except StopIteration:
                        tan_duplications_end = True
                        break
            if not tan_duplications_end:
                length2 = end2 - start2
                #if overlapping tandem duplication of similar length
                if contig2 == contig1 and start2 < end1 and (length1 - length2) / max(length1, length2) < 0.2:
                    inserted_regions_to_remove_2.append(inserted_region_index)

    # remove found inserted regions
    all_inserted_regions_to_remove = sorted(list(set(inserted_regions_to_remove_1 + inserted_regions_to_remove_2)), reverse=True)
    for ins_index in all_inserted_regions_to_remove:
        del(insertion_signature_clusters[ins_index])

    ##############################
    # Create deletion candidates #
    ##############################
    deletion_candidates = []
    for del_cluster in deletion_signature_clusters:
        if del_cluster.score > 0:
            deletion_candidates.append(CandidateDeletion(del_cluster.contig, del_cluster.start, del_cluster.end, del_cluster.members, del_cluster.score, del_cluster.std_span, del_cluster.std_pos))

    #####################################
    # Create novel insertion candidates #
    #####################################
    novel_insertion_candidates = prepare_insertion_candidates(insertion_signature_clusters, options)

    ######################
    # Cluster candidates #
    ######################
    logging.info("Cluster interspersed duplication candidates one more time..")
    final_int_duplication_candidates = partition_and_cluster_candidates(int_duplication_candidates, options, "interspersed duplication candidates")

    return (deletion_candidates, inversion_candidates, final_int_duplication_candidates, tan_dup_candidates, novel_insertion_candidates, breakend_candidates)
