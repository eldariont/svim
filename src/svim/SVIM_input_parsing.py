import sys
import os
import logging
import argparse


def parse_arguments(program_version, arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""SVIM (pronounced SWIM) is a structural variant caller for long reads. 
It discriminates five different variant classes: deletions, tandem and interspersed duplications, 
inversions and insertions. SVIM is unique in its capability of extracting both the genomic origin and 
destination of duplications.

SVIM consists of four major steps:
- COLLECT detects signatures for SVs in long read alignments
- CLUSTER merges signatures that come from the same SV
- COMBINE combines clusters from different genomic regions and classifies them into distinct SV types
- GENOTYPE uses alignments spanning SVs to determine their genotype

SVIM can process two types of input. Firstly, it can detect SVs from raw reads by aligning them to a given reference genome first ("SVIM.py reads [options] working_dir reads genome").
Alternatively, it can detect SVs from existing reads alignments in SAM/BAM format ("SVIM.py alignment [options] working_dir bam_file").
""")

    subparsers = parser.add_subparsers(help='modes', dest='sub')
    parser.add_argument('--version',
                        '-v',
                        action='version',
                        version='%(prog)s {version}'.format(version=program_version))
    parser_fasta = subparsers.add_parser('reads',
                                          help='Detect SVs from raw reads. Align reads to given reference genome first.')
    parser_fasta.add_argument('working_dir',
                               type=str,
                               help='Working and output directory. \
                                     Existing files in the directory are overwritten. \
                                     If the directory does not exist, it is created.')
    parser_fasta.add_argument('reads',
                               type=str,
                               help='Read file (FASTA, FASTQ, gzipped FASTA, gzipped FASTQ or file list). \
                                     The read file has to have one of the following supported file endings: \
                                     FASTA: .fa, .fasta, .FA, .fa.gz, .fa.gzip, .fasta.gz, .fasta.gzip \
                                     FASTQ: .fq, .fastq, .FQ, .fq.gz, .fq.gzip, .fastq.gz, .fastq.gzip \
                                     FILE LIST: .fa.fn, fq.fn')
    parser_fasta.add_argument('genome',
                               type=str,
                               help='Reference genome file (FASTA)')
    parser_fasta.add_argument('--verbose',
                              action='store_true',
                              help='Enable more verbose logging (default: %(default)s)')
    group_fasta_align = parser_fasta.add_argument_group('ALIGN')
    group_fasta_align.add_argument('--cores',
                                        type=int,
                                        default=1,
                                        help='CPU cores to use for the alignment (default: %(default)s)')
    group_fasta_align.add_argument('--aligner',
                                        type=str,
                                        default="ngmlr",
                                        choices=["ngmlr", "minimap2"],
                                        help='Tool for read alignment: ngmlr or minimap2 (default: %(default)s)')
    group_fasta_align.add_argument('--nanopore',
                                        action='store_true',
                                        help='Use Nanopore settings for read alignment (default: %(default)s)')
    group_fasta_collect = parser_fasta.add_argument_group('COLLECT')
    group_fasta_collect.add_argument('--min_mapq',
                                        type=int,
                                        default=20,
                                        help='Minimum mapping quality of reads to consider (default: %(default)s). \
                                              Reads with a lower mapping quality are ignored.')
    group_fasta_collect.add_argument('--min_sv_size',
                                        type=int,
                                        default=40,
                                        help='Minimum SV size to detect (default: %(default)s). \
                                              SVIM can potentially detect events of any size but is limited by the \
                                              signal-to-noise ratio in the input alignments. That means that more \
                                              accurate reads and alignments enable the detection of smaller events. \
                                              For current PacBio or Nanopore data, we would recommend a minimum size \
                                              of 40bp or larger.')
    group_fasta_collect.add_argument('--max_sv_size',
                                        type=int,
                                        default=100000,
                                        help='Maximum SV size to detect (default: %(default)s). \
                                              This parameter is used to distinguish long deletions (and inversions) from \
                                              translocations which cannot be distinguished from the alignment alone. \
                                              Split read segments mapping far apart on the reference could either \
                                              indicate a very long deletion (inversion) or a translocation breakpoint. \
                                              SVIM calls a translocation breakpoint if the mapping distance is larger \
                                              than this parameter and a deletion (or inversion) if it is smaller or equal.')
    group_fasta_collect.add_argument('--segment_gap_tolerance',
                                        type=int,
                                        default=10,
                                        help='Maximum tolerated gap between adjacent alignment segments (default: %(default)s). \
                                              This parameter applies to gaps on the reference and the read. Example: \
                                              Deletions are detected from two subsequent segments of a split read that are mapped \
                                              far apart from each other on the reference. The segment gap tolerance determines \
                                              the maximum tolerated length of the read gap between both segments. If there is an \
                                              unaligned read segment larger than this value between the two segments, no deletion is called.')
    group_fasta_collect.add_argument('--segment_overlap_tolerance',
                                        type=int,
                                        default=5,
                                        help='Maximum tolerated overlap between adjacent alignment segments (default: %(default)s). \
                                              This parameter applies to overlaps on the reference and the read. Example: \
                                              Deletions are detected from two subsequent segments of a split read that are mapped \
                                              far apart from each other on the reference. The segment overlap tolerance determines \
                                              the maximum tolerated length of an overlap between both segments on the read. If the \
                                              overlap between the two segments on the read is larger than this value, no deletion is called.')
    group_fasta_collect.add_argument('--all_bnds',
                                        action='store_true',
                                        help='Output all rearrangements additionally in BND notation (default: %(default)s). \
                                              By default, SV signatures from the read alignments are used to detect complete SVs, \
                                              such as deletions, insertions and inversions. When this option is enabled, all SVs \
                                              are also output in breakend (BND) notation as defined in the VCF specs. For instance, \
                                              a deletion gets two records in the VCF output: 1. the normal <DEL> record and 2. \
                                              a <BND> record representing the novel adjacency between the deletion\'s start and \
                                              end coordinate in the sample genome.')
    group_fasta_cluster = parser_fasta.add_argument_group('CLUSTER')
    group_fasta_cluster.add_argument('--partition_max_distance',
                                        type=int,
                                        default=1000,
                                        help='Maximum distance in bp between SVs in a partition (default: %(default)s). \
                                              Before clustering, the SV signatures are divided into coarse partitions. This parameter \
                                              determines the maximum distance between two subsequent signatures in the same partition. \
                                              If the distance between two subsequent signatures \
                                              is larger than this parameter, they are distributed into separate partitions.')
    group_fasta_cluster.add_argument('--distance_normalizer',
                                        type=int,
                                        default=900,
                                        help='Distance normalizer used for span-position distance (default: %(default)s). \
                                              SVIM clusters the SV signatures using an hierarchical clustering approach and a \
                                              novel distance metric called \"span-position distance\". Span-position distance \
                                              is the sum of two components, span distance and position distance. \
                                              The span distance is the difference in lengths between signatures normalized \
                                              by the greater length and always lies in the interval [0,1]. \
                                              The position distance is the difference in position between signatures \
                                              normalized by the distance normalizer (this parameter). For a position difference \
                                              of 1.8kb and a distance normalizer of 900, the position distance will be 2. \
                                              A smaller distance normalizer leads to a higher position distance and as a \
                                              consequence increases the importance of the position distance in the \
                                              span-position distance relative to the span distance.')
    group_fasta_cluster.add_argument('--cluster_max_distance',
                                        type=float,
                                        default=0.3,
                                        help='Maximum span-position distance between SVs in a cluster (default: %(default)s). \
                                              This is the most important parameter because it determines the strictness \
                                              of clustering. Choosing a large value leads to fewer but larger clusters with larger \
                                              distances between its members. Choosing a small value leads to more but smaller \
                                              clusters with smaller distances between its members. \
                                              This parameter determines the height of the cut-off in the hierarchical clustering \
                                              dendrogram.')
    group_fasta_combine = parser_fasta.add_argument_group('COMBINE')
    group_fasta_combine.add_argument('--del_ins_dup_max_distance',
                                        type=float,
                                        default=1.0,
                                        help='Maximum span-position distance between the origin of an insertion and a deletion to be flagged as a potential cut&paste insertion (default: %(default)s)')
    group_fasta_combine.add_argument('--trans_sv_max_distance',
                                        type=int,
                                        default=500,
                                        help='Maximum distance in bp between a translocation breakpoint and an SV signature to be combined (default: %(default)s)')
    group_fasta_genotype = parser_fasta.add_argument_group('GENOTYPE')
    group_fasta_genotype.add_argument('--skip_genotyping',
                                        action='store_true',
                                        help='Disable genotyping (default: %(default)s)')
    group_fasta_genotype.add_argument('--minimum_score',
                                        type=int,
                                        default=3,
                                        help='Minimum score for genotyping (default: %(default)s). \
                                              Only SV candidates with a higher or equal score are genotyped. Depending on the \
                                              score distribution among the SV candidates, decreasing this value increases the \
                                              runtime. We recommend to choose a value close to the score threshold used \
                                              for filtering the SV candidates.')
    group_fasta_genotype.add_argument('--homozygous_threshold',
                                        type=float,
                                        default=0.8,
                                        help='Minimum variant allele frequency to be called as homozygous (default: %(default)s). \
                                              Allele frequency is computed as the fraction of reads supporting the variant over the \
                                              total number of reads covering the variant. Variants with an allele frequence greater \
                                              than or equal to this threshold are called as homozygous alternative.')
    group_fasta_genotype.add_argument('--heterozygous_threshold',
                                        type=float,
                                        default=0.2,
                                        help='Minimum variant allele frequency to be called as heterozygous (default: %(default)s). \
                                              Allele frequency is computed as the fraction of reads supporting the variant over the \
                                              total number of reads covering the variant. Variants with an allele frequence greater \
                                              than or equal to this threshold but lower than the homozygous threshold are called as \
                                              heterozygous alternative. Variants with an allele frequence lower than this threshold \
                                              are called as homozygous reference.')
    group_fasta_genotype.add_argument('--minimum_depth',
                                        type=int,
                                        default=4,
                                        help='Minimum total read depth for genotyping (default: %(default)s). \
                                              Variants covered by a total number of reads lower than this value are not assigned \
                                              a genotype (./. in the output VCF file).')
    group_fasta_output = parser_fasta.add_argument_group('OUTPUT')
    group_fasta_output.add_argument('--sample',
                                        type=str,
                                        default="Sample",
                                        help='Sample ID to include in output vcf file (default: %(default)s)')
    group_fasta_output.add_argument('--types',
                                        type=str,
                                        default="DEL,INS,INV,DUP:TANDEM,DUP:INT,BND",
                                        help='SV types to include in output VCF (default: %(default)s). \
                                              Give a comma-separated list of SV types. The possible SV types are: DEL (deletions), \
                                              INS (novel insertions), INV (inversions), DUP:TANDEM (tandem duplications), \
                                              DUP:INT (interspersed duplications), BND (breakends).')
    group_fasta_output.add_argument('--sequence_alleles',
                                        action='store_true',
                                        help='Use nucleotide sequences for alleles of deletions, inversions and insertions in output VCF (default: %(default)s). \
                                              By default, all SVs are represented by symbolic alleles, such as <DEL>, <INV> or <INS>. \
                                              If enabled, ALT alleles of insertions are obtained from the sequence of a random read that supports the variant.')
    group_fasta_output.add_argument('--insertion_sequences',
                                        action='store_true',
                                        help='Output insertion sequences in INFO tag of VCF (default: %(default)s). \
                                              If enabled, the INFO/SEQS tag contains a list of insertion sequences from the supporting reads. \
                                              However, the insertion sequences are not combined into a consensus sequence.')
    group_fasta_output.add_argument('--tandem_duplications_as_insertions',
                                        action='store_true',
                                        help='Represent tandem duplications as insertions in output VCF (default: %(default)s). \
                                              By default, tandem duplications are represented by the SVTYPE=DUP:TANDEM and the genomic source is given by the \
                                              POS and END tags. When enabling this option, duplications are instead represented by the SVTYPE=INS \
                                              and POS and END both give the insertion point of the duplication.')
    group_fasta_output.add_argument('--interspersed_duplications_as_insertions',
                                        action='store_true',
                                        help='Represent interspersed duplications as insertions in output VCF (default: %(default)s). \
                                              By default, interspersed duplications are represented by the SVTYPE=DUP:INT and the genomic source is given by the \
                                              POS and END tags. When enabling this option, duplications are instead represented by the SVTYPE=INS \
                                              and POS and END both give the insertion point of the duplication.')
    group_fasta_output.add_argument('--read_names',
                                        action='store_true',
                                        help='Output names of supporting reads in INFO tag of VCF (default: %(default)s). \
                                              If enabled, the INFO/READS tag contains the list of names of the supporting reads.')
    group_fasta_output.add_argument('--zmws',
                                        action='store_true',
                                        help='look for information on ZMWs in PacBio read names (default: %(default)s). \
                                              If enabled, the INFO/ZMWS tag contains the number of ZMWs that produced supporting reads.')

    parser_bam = subparsers.add_parser('alignment',
                                        help='Detect SVs from an existing alignment')
    parser_bam.add_argument('working_dir',
                             type=os.path.abspath,
                             help='Working and output directory. \
                                   Existing files in the directory are overwritten. \
                                   If the directory does not exist, it is created.')
    parser_bam.add_argument('bam_file',
                             type=str,
                             help='Coordinate-sorted and indexed BAM file with aligned long reads')
    parser_bam.add_argument('genome',
                               type=str,
                               help='Reference genome file that the long reads were aligned to (FASTA)')
    parser_bam.add_argument('--verbose',
                              action='store_true',
                              help='Enable more verbose logging (default: %(default)s)')
    group_bam_collect = parser_bam.add_argument_group('COLLECT')
    group_bam_collect.add_argument('--min_mapq',
                                      type=int,
                                      default=20,
                                      help='Minimum mapping quality of reads to consider (default: %(default)s). \
                                            Reads with a lower mapping quality are ignored.')
    group_bam_collect.add_argument('--min_sv_size',
                                      type=int,
                                      default=40,
                                      help='Minimum SV size to detect (default: %(default)s). \
                                            SVIM can potentially detect events of any size but is limited by the \
                                            signal-to-noise ratio in the input alignments. That means that more \
                                            accurate reads and alignments enable the detection of smaller events. \
                                            For current PacBio or Nanopore data, we would recommend a minimum size \
                                            of 40bp or larger.')
    group_bam_collect.add_argument('--max_sv_size',
                                      type=int,
                                      default=100000,
                                      help='Maximum SV size to detect (default: %(default)s). \
                                              This parameter is used to distinguish long deletions (and inversions) from \
                                              translocations which cannot be distinguished from the alignment alone. \
                                              Split read segments mapping far apart on the reference could either \
                                              indicate a very long deletion (inversion) or a translocation breakpoint. \
                                              SVIM calls a translocation breakpoint if the mapping distance is larger \
                                              than this parameter and a deletion (or inversion) if it is smaller or equal.')
    group_bam_collect.add_argument('--segment_gap_tolerance',
                                      type=int,
                                      default=10,
                                      help='Maximum tolerated gap between adjacent alignment segments (default: %(default)s). \
                                            This parameter applies to gaps on the reference and the read. Example: \
                                            Deletions are detected from two subsequent segments of a split read that are mapped \
                                            far apart from each other on the reference. The segment gap tolerance determines \
                                            the maximum tolerated length of the read gap between both segments. If there is an \
                                            unaligned read segment larger than this value between the two segments, no deletion is called.')
    group_bam_collect.add_argument('--segment_overlap_tolerance',
                                      type=int,
                                      default=5,
                                      help='Maximum tolerated overlap between adjacent alignment segments (default: %(default)s). \
                                            This parameter applies to overlaps on the reference and the read. Example: \
                                            Deletions are detected from two subsequent segments of a split read that are mapped \
                                            far apart from each other on the reference. The segment overlap tolerance determines \
                                            the maximum tolerated length of an overlap between both segments on the read. If the \
                                            overlap between the two segments on the read is larger than this value, no deletion is called.')
    group_bam_cluster = parser_bam.add_argument_group('CLUSTER')
    group_bam_cluster.add_argument('--partition_max_distance',
                                      type=int,
                                      default=1000,
                                      help='Maximum distance in bp between SVs in a partition (default: %(default)s). \
                                            Before clustering, the SV signatures are divided into coarse partitions. This parameter \
                                            determines the maximum distance between two subsequent signatures in the same partition. \
                                            If the distance between two subsequent signatures \
                                            is larger than this parameter, they are distributed into separate partitions.')
    group_bam_cluster.add_argument('--distance_normalizer',
                                      type=int,
                                      default=900,
                                      help='Distance normalizer used for span-position distance (default: %(default)s). \
                                            SVIM clusters the SV signatures using an hierarchical clustering approach and a \
                                            novel distance metric called \"span-position distance\". Span-position distance \
                                            is the sum of two components, span distance and position distance. \
                                            The span distance is the difference in lengths between signatures normalized \
                                            by the greater length and always lies in the interval [0,1]. \
                                            The position distance is the difference in position between signatures \
                                            normalized by the distance normalizer (this parameter). For a position difference \
                                            of 1.8kb and a distance normalizer of 900, the position distance will be 2. \
                                            A smaller distance normalizer leads to a higher position distance and as a \
                                            consequence increases the importance of the position distance in the \
                                            span-position distance relative to the span distance.')
    group_bam_cluster.add_argument('--cluster_max_distance',
                                      type=float,
                                      default=0.3,
                                      help='Maximum span-position distance between SVs in a cluster (default: %(default)s). \
                                            This is the most important parameter because it determines the strictness \
                                            of clustering. Choosing a large value leads to fewer but larger clusters with larger \
                                            distances between its members. Choosing a small value leads to more but smaller \
                                            clusters with smaller distances between its members. \
                                            This parameter determines the height of the cut-off in the hierarchical clustering \
                                            dendrogram.')
    group_bam_cluster.add_argument('--all_bnds',
                                        action='store_true',
                                        help='Output all rearrangements additionally in BND notation (default: %(default)s). \
                                              By default, SV signatures from the read alignments are used to detect complete SVs, \
                                              such as deletions, insertions and inversions. When this option is enabled, all SVs \
                                              are also output in breakend (BND) notation as defined in the VCF specs. For instance, \
                                              a deletion gets two records in the VCF output: 1. the normal <DEL> record and 2. \
                                              a <BND> record representing the novel adjacency between the deletion\'s start and \
                                              end coordinate in the sample genome.')
    group_bam_combine = parser_bam.add_argument_group('COMBINE')
    group_bam_combine.add_argument('--del_ins_dup_max_distance',
                                      type=float,
                                      default=1.0,
                                      help='Maximum span-position distance between the origin of an insertion and a deletion to be flagged as a potential cut&paste insertion (default: %(default)s)')
    group_bam_combine.add_argument('--trans_sv_max_distance',
                                      type=int,
                                      default=500,
                                      help='Maximum distance in bp between a translocation breakpoint and an SV signature to be combined (default: %(default)s)')
    group_bam_genotype = parser_bam.add_argument_group('GENOTYPE')
    group_bam_genotype.add_argument('--skip_genotyping',
                                        action='store_true',
                                        help='Disable genotyping (default: %(default)s)')
    group_bam_genotype.add_argument('--minimum_score',
                                      type=int,
                                      default=3,
                                      help='Minimum score for genotyping (default: %(default)s). \
                                            Only SV candidates with a higher or equal score are genotyped. Depending on the \
                                            score distribution among the SV candidates, decreasing this value increases the \
                                            runtime. We recommend to choose a value close to the score threshold used \
                                            for filtering the SV candidates.')
    group_bam_genotype.add_argument('--homozygous_threshold',
                                      type=float,
                                      default=0.8,
                                      help='Minimum variant allele frequency to be called as homozygous (default: %(default)s). \
                                            Allele frequency is computed as the fraction of reads supporting the variant over the \
                                            total number of reads covering the variant. Variants with an allele frequence greater \
                                            than or equal to this threshold are called as homozygous alternative.')
    group_bam_genotype.add_argument('--heterozygous_threshold',
                                      type=float,
                                      default=0.2,
                                      help='Minimum variant allele frequency to be called as heterozygous (default: %(default)s). \
                                            Allele frequency is computed as the fraction of reads supporting the variant over the \
                                            total number of reads covering the variant. Variants with an allele frequence greater \
                                            than or equal to this threshold but lower than the homozygous threshold are called as \
                                            heterozygous alternative. Variants with an allele frequence lower than this threshold \
                                            are called as homozygous reference.')
    group_bam_genotype.add_argument('--minimum_depth',
                                      type=int,
                                      default=4,
                                      help='Minimum total read depth for genotyping (default: %(default)s). \
                                            Variants covered by a total number of reads lower than this value are not assigned \
                                            a genotype (./. in the output VCF file).')
    group_bam_output = parser_bam.add_argument_group('OUTPUT')
    group_bam_output.add_argument('--sample',
                                        type=str,
                                        default="Sample",
                                        help='Sample ID to include in output vcf file (default: %(default)s)')
    group_bam_output.add_argument('--types',
                                        type=str,
                                        default="DEL,INS,INV,DUP:TANDEM,DUP:INT,BND",
                                        help='SV types to include in output VCF (default: %(default)s). \
                                              Give a comma-separated list of SV types. The possible SV types are: DEL (deletions), \
                                              INS (novel insertions), INV (inversions), DUP:TANDEM (tandem duplications), \
                                              DUP:INT (interspersed duplications), BND (breakends).')
    group_bam_output.add_argument('--sequence_alleles',
                                        action='store_true',
                                        help='Use nucleotide sequences for alleles of deletions, inversions and insertions in output VCF (default: %(default)s). \
                                              By default, all SVs are represented by symbolic alleles, such as <DEL>, <INV> or <INS>. \
                                              If enabled, ALT alleles of insertions are obtained from the sequence of a random read that supports the variant.')
    group_bam_output.add_argument('--insertion_sequences',
                                        action='store_true',
                                        help='Output insertion sequences in INFO tag of VCF (default: %(default)s). \
                                              If enabled, the INFO/SEQS tag contains a list of insertion sequences from the supporting reads. \
                                              However, the insertion sequences are not combined into a consensus sequence.')
    group_bam_output.add_argument('--tandem_duplications_as_insertions',
                                        action='store_true',
                                        help='Represent tandem duplications as insertions in output VCF (default: %(default)s). \
                                              By default, tandem duplications are represented by the SVTYPE=DUP:TANDEM and the genomic source is given by the \
                                              POS and END tags. When enabling this option, duplications are instead represented by the SVTYPE=INS \
                                              and POS and END both give the insertion point of the duplication.')
    group_bam_output.add_argument('--interspersed_duplications_as_insertions',
                                        action='store_true',
                                        help='Represent interspersed duplications as insertions in output VCF (default: %(default)s). \
                                              By default, interspersed duplications are represented by the SVTYPE=DUP:INT and the genomic source is given by the \
                                              POS and END tags. When enabling this option, duplications are instead represented by the SVTYPE=INS \
                                              and POS and END both give the insertion point of the duplication.')
    group_bam_output.add_argument('--read_names',
                                        action='store_true',
                                        help='Output names of supporting reads in INFO tag of VCF (default: %(default)s). \
                                              If enabled, the INFO/READS tag contains the list of names of the supporting reads.')
    group_bam_output.add_argument('--zmws',
                                        action='store_true',
                                        help='look for information on ZMWs in PacBio read names (default: %(default)s). \
                                              If enabled, the INFO/ZMWS tag contains the number of ZMWs that produced supporting reads.')

    return parser.parse_args(arguments)


def guess_file_type(reads_path):
    if reads_path.endswith(".fa") or reads_path.endswith(".fasta") or reads_path.endswith(".FA"):
        logging.info("Recognized reads file as FASTA format.")
        return "fasta"
    elif reads_path.endswith(".fq") or reads_path.endswith(".fastq") or reads_path.endswith(".FQ"):
        logging.info("Recognized reads file as FASTQ format.")
        return "fastq"
    elif reads_path.endswith(".fa.gz") or reads_path.endswith(".fasta.gz") or reads_path.endswith(".FA.gz") or reads_path.endswith(".fa.gzip") or reads_path.endswith(".fasta.gzip") or reads_path.endswith(".FA.gzip"):
        logging.info("Recognized reads file as gzipped FASTA format.")
        return "fasta_gzip"
    elif reads_path.endswith(".fq.gz") or reads_path.endswith(".fastq.gz") or reads_path.endswith(".FQ.gz") or reads_path.endswith(".fq.gzip") or reads_path.endswith(".fastq.gzip") or reads_path.endswith(".FQ.gzip"):
        logging.info("Recognized reads file as gzipped FASTQ format.")
        return "fastq_gzip"
    elif reads_path.endswith(".fa.fn") or reads_path.endswith(".fasta.fn") or reads_path.endswith(".FA.fn") or reads_path.endswith(".fq.fn") or reads_path.endswith(".fastq.fn") or reads_path.endswith(".FQ.fn"):
        logging.info("Recognized reads file as file list format.")
        return "list"
    else:
        logging.error("Unknown file ending of file {0}. See github.com/eldariont/svim/wiki/ for supported file endings. Exiting.".format(reads_path))
        return "unknown"


def read_file_list(path):
    file_list = open(path, "r")
    for line in file_list:
        yield line.strip()
    file_list.close()