SVIM - Structural variant identification using long reads
=========================================================

.. image:: https://badge.fury.io/py/svim.svg
    :target: https://badge.fury.io/py/svim

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
    :target: http://bioconda.github.io

SVIM (pronounced *SWIM*) is a structural variant caller for long reads.
It is able to detect, classify and genotype five different classes of structural variants.
Unlike existing methods, SVIM integrates information from across the genome to precisely distinguish similar events, such as tandem and interspersed duplications and insertions.
In our experiments on simulated data and real datasets from PacBio and Nanopore sequencing machines, SVIM reached consistently better results than competing methods.
Furthermore, it is unique in its capability of extracting both the genomic origin and destination of duplications.

Background on Structural Variants and Long Reads
------------------------------------------------

.. image:: https://raw.githubusercontent.com/eldariont/svim/master/docs/SVclasses.png
    :align: center

Structural variants (SVs) are typically defined as genomic variants larger than 50bps (e.g. deletions, duplications, inversions).
Studies have shown that they affect more bases in any given genome than SNPs and small Indels taken together.
Consequently, they have a large impact on genes and regulatory regions.
This is reflected in the large number of genetic diseases that are caused by SVs.

Common sequencing technologies by providers such as Illumina generate short reads with high accuracy.
However, they exhibit weaknesses in repeat and low-complexity regions.
This negatively affects SV detection because SVs are associated to such regions.
Single molecule long-read sequencing technologies from Pacific Biotechnologies and Oxford Nanopore produce reads with error rates of up to 15% but with lengths of several kb.
The high read lengths enable them to cover entire repeats and SVs which facilitates SV detection.

Installation
------------

.. code-block:: bash

    #Install via conda into a new environment (recommended): installs all dependencies including read alignment dependencies
    conda create -n svim_env --channel bioconda svim

    #Install via conda into existing (active) environment: installs all dependencies including read alignment dependencies
    conda install --channel bioconda svim

    #Install via pip (requires Python 3.6.*): installs all dependencies except those necessary for read alignment (ngmlr, minimap2, samtools)
    pip3 install svim

    #Install from github (requires Python 3.6.*): installs all dependencies except those necessary for read alignment (ngmlr, minimap2, samtools)
    git clone https://github.com/eldariont/svim.git
    cd svim
    pip3 install .

Changelog
---------
- **v1.4.0**: improve partitioning and clustering of insertions, update default value of --partition_max_distance to avoid very large partitions, add --verbose option
- **v1.3.1**: small changes to partitioning and clustering algorithm, add two new command-line options to output duplications as INS records in VCF, remove limit on number of supplementary alignments, remove q5 filter, bugfixes
- **v1.3.0**: improve BND detection, add INFO:ZMWS tag with number of supporting PacBio wells, add sequence alleles for INS, add FORMAT:CN tag for tandem duplications, bugfixes
- **v1.2.0**: add 3 more VCF output options: output sequence instead of symbolic alleles in VCF, output names of supporting reads, output insertion sequences of supporting reads
- **v1.1.0**: outputs BNDs in VCF, detects large tandem duplications, allows skipping genotyping, makes VCF output more flexible, adds genotype scatter plot
- **v1.0.0**: adds genotyping of deletions, inversions, insertions and interspersed duplications, produces plots of SV length distribution, improves help descriptions
- **v0.5.0**: replaces graph-based clustering with hierarchical clustering, modifies scoring function, improves partitioning prior to clustering, improves calling from coordinate-sorted SAM/BAM files, improves VCF output
- **v0.4.4**: includes exception message into log files, bug fixes, adds tests and sets up Travis
- **v0.4.3**: adds support for coordinate-sorted SAM/BAM files, improves VCF output and increases compatibility with IGV and truvari, bug fixes
    
Input
-----

SVIM analyzes long reads given as a FASTA/FASTQ file (uncompressed or gzipped) or a file list.
Alternatively, it can analyze an alignment file in BAM format.
SVIM was tested on both PacBio and Nanopore data.
It works best for alignment files produced by `NGMLR <https://github.com/philres/ngmlr>`_ but also supports the faster read mapper `minimap2 <https://github.com/lh3/minimap2>`_.

Output
------

SVIM's main output file called `variants.vcf` (formerly final_results.vcf) is placed into the given working directory.
For each of the five detected SV classes, SVIM also produces a BED file with the SV coordinates in the `candidates` subdirectory.

Usage
----------------------

Please see our `wiki <https://github.com/eldariont/svim/wiki>`_.

Contact
-------

If you experience problems or have suggestions please create an issue or a pull request or contact heller_d@molgen.mpg.de.

Citation
---------

Feel free to read and cite our paper in Bioinformatics: https://doi.org/10.1093/bioinformatics/btz041

License
-------

The project is licensed under the GNU General Public License.
