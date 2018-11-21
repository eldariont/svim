import logging
import os

from subprocess import Popen, PIPE

from svim.SVIM_fullread import analyze_full_read_indel
from svim.SVIM_splitread import analyze_full_read_segments

def guess_file_type(reads_path):
    if reads_path.endswith(".fa") or reads_path.endswith(".fasta") or reads_path.endswith(".FA"):
        logging.info("Recognized reads file as FASTA format.")
        return "fasta"
    elif reads_path.endswith(".fq") or reads_path.endswith(".fastq") or reads_path.endswith(".FQ"):
        logging.info("Recognized reads file as FASTQ format.")
        return "fastq"
    elif reads_path.endswith(".fa.gz") or reads_path.endswith(".fasta.gz") or reads_path.endswith(".fa.gzip") or reads_path.endswith(".fasta.gzip"):
        logging.info("Recognized reads file as gzipped FASTA format.")
        return "fasta_gzip"
    elif reads_path.endswith(".fq.gz") or reads_path.endswith(".fastq.gz") or reads_path.endswith(".fq.gzip") or reads_path.endswith(".fastq.gzip"):
        logging.info("Recognized reads file as gzipped FASTQ format.")
        return "fastq_gzip"
    elif reads_path.endswith(".fa.fn"):
        logging.info("Recognized reads file as FASTA file list format.")
        return "list" 
    else:
        logging.error("Unknown file ending of file {0}. Exiting.".format(reads_path))
        return "unknown"


def create_full_file(working_dir, reads_path, reads_type):
    """Create FASTA file with full reads if it does not exist."""
    if not os.path.exists(working_dir):
        logging.error("Given working directory does not exist")
        sys.exit()

    reads_file_prefix = os.path.splitext(os.path.basename(reads_path))[0]

    if reads_type == "fasta_gzip" or reads_type == "fastq_gzip":
        full_reads_path = "{0}/{1}.fa".format(working_dir, reads_file_prefix)
        if not os.path.exists(full_reads_path):
            full_file = open(full_reads_path, 'w')
            write_full = True
        else:
            write_full = False
            logging.warning("FASTA file for full reads exists. Skip")
    else:
        full_reads_path = reads_path
        write_full = False

    if write_full:
        if reads_type == "fasta" or reads_type == "fastq":
            reads_file = open(reads_path, "r")
        elif reads_type == "fasta_gzip" or reads_type == "fastq_gzip":
            reads_file = gzip.open(reads_path, "rt")

        if reads_type == "fasta" or reads_type == "fasta_gzip":
            sequence = ""
            for line in reads_file:
                if line.startswith('>'):
                    if sequence != "":
                        if write_full:
                            print(">" + read_name, file=full_file)
                            print(sequence, file=full_file)
                    read_name = line.strip()[1:]
                    sequence = ""
                else:
                    sequence += line.strip()
            if sequence != "":
                if write_full:
                    print(">" + read_name, file=full_file)
                    print(sequence, file=full_file)
            reads_file.close()
        elif reads_type == "fastq" or reads_type == "fastq_gzip":
            sequence_line_is_next = False
            for line in reads_file:
                if line.startswith('@'):
                    read_name = line.strip()[1:]
                    sequence_line_is_next = True
                elif sequence_line_is_next:
                    sequence_line_is_next = False
                    sequence = line.strip()
                    if write_full:
                        print(">" + read_name, file=full_file)
                        print(sequence, file=full_file)
            reads_file.close()

    if write_full:
        full_file.close()

    logging.info("Full read files written")
    return full_reads_path


def run_full_alignment(working_dir, genome, reads_path, cores, aligner, nanopore):
    """Align full reads with NGM-LR."""
    reads_file_prefix = os.path.splitext(os.path.basename(reads_path))[0]
    full_aln = "{0}/{1}_aln.querysorted.bam".format(working_dir, reads_file_prefix)

    if not os.path.exists(full_aln):
        if aligner == "ngmlr":
            if nanopore:
                alignment = Popen(['ngmlr',
                           '-t', str(cores), '-r', genome, '-q', os.path.realpath(reads_path), '-x', 'ont'], stdout=PIPE)
            else:
                alignment = Popen(['ngmlr',
                           '-t', str(cores), '-r', genome, '-q', os.path.realpath(reads_path)], stdout=PIPE)
        elif aligner == "minimap2":
            if nanopore:
                alignment = Popen(['minimap2',
                           '-t', str(cores), '-x', 'map-ont', '-a', genome, os.path.realpath(reads_path)], stdout=PIPE)
            else:
                alignment = Popen(['minimap2',
                           '-t', str(cores), '-x', 'map-pb', '-a', genome, os.path.realpath(reads_path)], stdout=PIPE)
        view = Popen(['samtools',
                      'view', '-b', '-@', str(cores)], stdin=alignment.stdout, stdout=PIPE)
        sort = Popen(['samtools',
                      'sort', '-n', '-@', str(cores), '-o', full_aln],
                     stdin=view.stdout)
        sort.wait()
        logging.info("Alignment finished")
    else:
        logging.warning("Alignment for full sequences exists. Skip")


# def parse_sam_file(sam, contig):
#     """Parses a SAM file and returns a dict of reads (list of alignments for each read) for a given reference contig"""
#     alns = sam.fetch(reference=contig)
#     aln_dict = defaultdict(list)
#     for aln in alns:
#         aln_dict[aln.query_name].append(aln)
#     return aln_dict


# def natural_representation(qname): 
#     """Splits a read name into a tuple of strings and numbers. This facilitates the sort order applied by samtools -n
#        See https://www.biostars.org/p/102735/"""
#     return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', qname)]


def bam_iterator(bam):
    """Returns an iterator for the given SAM/BAM file (must be query-sorted). 
    In each call, the alignments of a single read are yielded as a 4-tuple: (read_name, list of primary pysam.AlignedSegment, list of supplementary pysam.AlignedSegment, list of secondary pysam.AlignedSegment)."""
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
                yield (current_read_name, current_prim, current_suppl, current_sec)
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
    yield (current_read_name, current_prim, current_suppl, current_sec)


def analyze_alignment(full_bam, options):
    full_it = bam_iterator(full_bam)

    sv_signatures = []
    read_nr = 0

    while True:
        try:
            full_iterator_object = next(full_it)
            read_nr += 1
            if read_nr % 10000 == 0:
                logging.info("Processed read {0}".format(read_nr))
            if not options.skip_indel:
                sv_signatures.extend(analyze_full_read_indel(full_iterator_object, full_bam, options))
            if not options.skip_segment:
                sv_signatures.extend(analyze_full_read_segments(full_iterator_object, full_bam, options))
        except StopIteration:
            break
        except KeyboardInterrupt:
            logging.warning('Execution interrupted by user. Stop detection and continue with next step..')
            break
    return sv_signatures


def read_file_list(path):
    file_list = open(path, "r")
    for line in file_list:
        yield line.strip()
    file_list.close()