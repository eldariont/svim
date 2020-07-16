import os
import logging

from svim.SVIM_clustering import partition_and_cluster


def cluster_sv_signatures(sv_signatures, options):
    """Takes a list of SVSignatures and splits them up by type. The SVSignatures of each type are clustered and returned as a tuple of
    (deletion_signature_clusters, insertion_signature_clusters, inversion_signature_clusters, tandem_duplication_signature_clusters, insertion_from_signature_clusters, completed_translocation_signatures)."""

    deletion_signatures = [ev for ev in sv_signatures if ev.type == "DEL"]
    insertion_signatures = [ev for ev in sv_signatures if ev.type == "INS"]
    inversion_signatures = [ev for ev in sv_signatures if ev.type == "INV"]
    tandem_duplication_signatures = [ev for ev in sv_signatures if ev.type == "DUP_TAN"]
    translocation_signatures = [ev for ev in sv_signatures if ev.type == "BND"]
    insertion_from_signatures = [ev for ev in sv_signatures if ev.type == "DUP_INT"]

    # Cluster SV signatures
    deletion_signature_clusters = partition_and_cluster(deletion_signatures, options, "deleted regions")
    insertion_signature_clusters = partition_and_cluster(insertion_signatures, options, "inserted regions")
    inversion_signature_clusters = partition_and_cluster(inversion_signatures, options, "inverted regions")
    tandem_duplication_signature_clusters = partition_and_cluster(tandem_duplication_signatures, options, "tandem duplicated regions")
    translocation_signature_clusters = partition_and_cluster(translocation_signatures, options, "translocation breakpoints")
    insertion_from_signature_clusters = partition_and_cluster(insertion_from_signatures, options, "inserted regions with detected region of origin")

    return (deletion_signature_clusters, insertion_signature_clusters, inversion_signature_clusters, tandem_duplication_signature_clusters, insertion_from_signature_clusters, translocation_signature_clusters)


def write_signature_clusters_bed(working_dir, clusters):
    """Write signature clusters into working directory in BED format."""
    deletion_signature_clusters, insertion_signature_clusters, inversion_signature_clusters, tandem_duplication_signature_clusters, insertion_from_signature_clusters, translocation_signature_clusters = clusters

    # Print SV signature clusters
    if not os.path.exists(working_dir + '/signatures'):
        os.mkdir(working_dir + '/signatures')
    deletion_signature_output = open(working_dir + '/signatures/del.bed', 'w')
    insertion_signature_output = open(working_dir + '/signatures/ins.bed', 'w')
    inversion_signature_output = open(working_dir + '/signatures/inv.bed', 'w')
    tandem_duplication_signature_source_output = open(working_dir + '/signatures/dup_tan_source.bed', 'w')
    tandem_duplication_signature_dest_output = open(working_dir + '/signatures/dup_tan_dest.bed', 'w')
    translocation_signature_output = open(working_dir + '/signatures/trans.bed', 'w')
    insertion_from_signature_output = open(working_dir + '/signatures/dup_int.bed', 'w')

    for cluster in deletion_signature_clusters:
        print(cluster.get_bed_entry(), file=deletion_signature_output)
    for cluster in insertion_signature_clusters:
        print(cluster.get_bed_entry(), file=insertion_signature_output)
    for cluster in inversion_signature_clusters:
        print(cluster.get_bed_entry(), file=inversion_signature_output)
    for cluster in tandem_duplication_signature_clusters:
        bed_entries = cluster.get_bed_entries()
        print(bed_entries[0], file=tandem_duplication_signature_source_output)
        print(bed_entries[1], file=tandem_duplication_signature_dest_output)
    for cluster in insertion_from_signature_clusters:
        bed_entries = cluster.get_bed_entries()
        print(bed_entries[0], file=insertion_from_signature_output)
        print(bed_entries[1], file=insertion_from_signature_output)
    for cluster in translocation_signature_clusters:
        bed_entries = cluster.get_bed_entries()
        print(bed_entries[0], file=translocation_signature_output)
        print(bed_entries[1], file=translocation_signature_output)

    deletion_signature_output.close()
    insertion_signature_output.close()
    inversion_signature_output.close()
    tandem_duplication_signature_source_output.close()
    tandem_duplication_signature_dest_output.close()
    translocation_signature_output.close()
    insertion_from_signature_output.close()


def write_signature_clusters_vcf(working_dir, clusters, version):
    """Write signature clusters into working directory in VCF format."""
    deletion_signature_clusters, insertion_signature_clusters, inversion_signature_clusters, tandem_duplication_signature_clusters, insertion_from_signature_clusters, translocation_signature_clusters = clusters

    if not os.path.exists(working_dir + '/signatures'):
        os.mkdir(working_dir + '/signatures')
    vcf_output = open(working_dir + '/signatures/all.vcf', 'w')

    # Write header lines
    print("##fileformat=VCFv4.3", file=vcf_output)
    print("##source=SVIMV{0}".format(version), file=vcf_output)
    print("##ALT=<ID=DEL,Description=\"Deletion\">", file=vcf_output)
    print("##ALT=<ID=INV,Description=\"Inversion\">", file=vcf_output)
    print("##ALT=<ID=DUP,Description=\"Duplication\">", file=vcf_output)
    print("##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">", file=vcf_output)
    print("##ALT=<ID=INS,Description=\"Insertion\">", file=vcf_output)
    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">", file=vcf_output)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">", file=vcf_output)
    print("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">", file=vcf_output)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=vcf_output)

    vcf_entries = []
    for cluster in deletion_signature_clusters:
        vcf_entries.append((cluster.get_source(), cluster.get_vcf_entry()))
    for cluster in insertion_signature_clusters:
        vcf_entries.append((cluster.get_source(), cluster.get_vcf_entry()))
    for cluster in inversion_signature_clusters:
        vcf_entries.append((cluster.get_source(), cluster.get_vcf_entry()))
    for cluster in tandem_duplication_signature_clusters:
        vcf_entries.append((cluster.get_source(), cluster.get_vcf_entry()))

    # Sort and write entries to VCF
    for source, entry in sorted(vcf_entries, key=lambda pair: pair[0]):
        print(entry, file=vcf_output)

    vcf_output.close()