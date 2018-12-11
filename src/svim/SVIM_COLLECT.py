import logging
import os

from subprocess import run, CalledProcessError

from svim.SVIM_fullread import analyze_full_read_indel
from svim.SVIM_splitread import analyze_full_read_segments


class ToolMissingError(Exception): pass

class AlignmentPipelineError(Exception): pass

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
    elif reads_path.endswith(".fa.fn") or reads_path.endswith(".fq.fn"):
        logging.info("Recognized reads file as file list format.")
        return "list"
    else:
        logging.error("Unknown file ending of file {0}. See github.com/eldariont/svim/wiki/ for supported file endings. Exiting.".format(reads_path))
        return "unknown"


def check_prereqisites(aligner):
    devnull = open(os.devnull, 'w')
    try:
        run(['gunzip', '--help'], stdout=devnull, stderr=devnull, check=True)
        run([aligner, '--help'], stdout=devnull, stderr=devnull, check=True)
        run(['samtools', '--help'], stdout=devnull, stderr=devnull, check=True)
    except FileNotFoundError as e:
        raise ToolMissingError('The alignment pipeline cannot be started because {0} was not found. Is it installed and in the PATH?'.format(e.filename)) from e
    except CalledProcessError as e:
        raise ToolMissingError('The alignment pipeline cannot be started because {0} failed.'.format(" ".join(e.cmd))) from e


def run_alignment(working_dir, genome, reads_path, reads_type, cores, aligner, nanopore):
    """Align full reads with NGMLR or minimap2."""
    check_prereqisites(aligner)
    reads_file_prefix = os.path.splitext(os.path.basename(reads_path))[0]
    full_aln = "{0}/{1}.{2}.querysorted.bam".format(working_dir, reads_file_prefix, aligner)
    if not os.path.exists(full_aln):
        try:
            command = ['set', '-o', 'pipefail', '&&']
            if aligner == "ngmlr":
                # We need to uncompress gzipped files for NGMLR first
                if reads_type == "fasta_gzip" or reads_type == "fastq_gzip":
                    command += ['gunzip', '-c', os.path.realpath(reads_path)]
                    command += ['|', 'ngmlr', '-t', str(cores), '-r', genome]
                    if nanopore:
                        command += ['-x', 'ont']
                else:
                    command += ['ngmlr', '-t', str(cores), '-r', genome, '-q', os.path.realpath(reads_path)]
                    if nanopore:
                        command += ['-x', 'ont']
            elif aligner == "minimap2":
                if nanopore:
                    command += ['minimap2', '-t', str(cores), '-x', 'map-ont', '-a', genome, os.path.realpath(reads_path)]
                else:
                    command += ['minimap2', '-t', str(cores), '-x', 'map-pb', '-a', genome, os.path.realpath(reads_path)]
            command += ['|', 'samtools', 'view', '-b', '-@', str(cores)]
            command += ['|', 'samtools', 'sort', '-n', '-@', str(cores), '-o', full_aln]
            logging.info("Starting alignment pipeline..")
            run(" ".join(command), shell=True, check=True)
        except CalledProcessError as e:
            raise AlignmentPipelineError('The alignment pipeline failed with exit code {0}. Command was: {1}'.format(e.returncode, e.cmd)) from e
        logging.info("Alignment pipeline finished")
        return full_aln
    else:
        logging.warning("Alignment output file {0} already exists. Skip alignment and use the existing file.".format(full_aln))
        return full_aln


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