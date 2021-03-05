SVIM - Structural variant identification using long reads
=========================================================

.. image:: https://img.shields.io/pypi/v/svim?style=flat
    :target: https://pypi.org/project/svim/

.. image:: https://img.shields.io/conda/vn/bioconda/svim?style=flat
    :target: https://anaconda.org/bioconda/svim

.. image:: https://img.shields.io/conda/dn/bioconda/svim?label=bioconda%20downloads&style=flat
    :target: https://anaconda.org/bioconda/svim

.. image:: https://img.shields.io/badge/published%20in-Bioinformatics-blue.svg
    :target: https://doi.org/10.1093/bioinformatics/btz041

SVIM (pronounced *SWIM*) is a structural variant caller for long sequencing reads.
It is able to detect and classify the following six classes of structural variation: deletions, insertions, inversions, tandem duplications, interspersed duplications and translocations.
SVIM also estimates the genotypes of deletions, insertions, inversions and interspersed duplications.
Unlike other methods, SVIM integrates information from across the genome to precisely distinguish similar events, such as tandem and interspersed duplications and simple insertions.
In our experiments on simulated data and real datasets from PacBio and Nanopore sequencing machines, SVIM reached consistently better results than competing methods.

**Note!** To analyze haploid or diploid genome assemblies or contigs, please use our other method `SVIM-asm <https://github.com/eldariont/svim-asm>`_.

Background on Structural Variants and Long Reads
------------------------------------------------

.. image:: https://raw.githubusercontent.com/eldariont/svim/master/docs/SVclasses.png
    :align: center

Structural variants (SVs) are typically defined as genomic variants larger than 50bps (e.g. deletions, duplications, inversions).
Studies have shown that they affect more bases in an average genome than SNPs or small Indels.
Consequently, they have a large impact on genes and regulatory regions.
This is reflected in the large number of genetic disorders and other disease that are associated to SVs.

Common sequencing technologies by providers such as Illumina generate short reads with high accuracy.
However, they exhibit weaknesses in repeat and low-complexity regions where SVs are particularly common.
Single molecule long-read sequencing technologies from Pacific Biotechnologies and Oxford Nanopore produce reads with error rates of up to 15% but with lengths of several kbps.
The high read lengths enable them to cover entire repeats and SVs which facilitates SV detection.

Installation
------------

.. code-block:: bash

    #Install via conda into a new environment (recommended): installs all dependencies including read alignment dependencies
    conda create -n svim_env --channel bioconda svim

    #Install via conda into existing (active) environment: installs all dependencies including read alignment dependencies
    conda install --channel bioconda svim

    #Install via pip (requires Python 3.6.* or newer): installs all dependencies except those necessary for read alignment (ngmlr, minimap2, samtools)
    pip install svim

    #Install from github (requires Python 3.6.* or newer): installs all dependencies except those necessary for read alignment (ngmlr, minimap2, samtools)
    git clone https://github.com/eldariont/svim.git
    cd svim
    pip install .

Changelog
---------
- **v1.5.0**: prevent signatures from same read to be clustered together, bugfixes
- **v1.4.2**: fix invalid start coordinates in VCF output, issue warning for invalid characters in contig names 
- **v1.4.1**: improve clustering of translocation breakpoints (BNDs), improve --all_bnds mode, bugfixes
- **v1.4.0**: fix and improve clustering of insertions, add option --all_bnds to output all SV classes in breakend notation, update default value of --partition_max_distance to avoid very large partitions, bugfixes
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
SVIM has been successfully tested on PacBio CLR, PacBio CCS (HiFi) and Oxford Nanopore data.
It works best for alignment files produced by `NGMLR <https://github.com/philres/ngmlr>`_ but also supports the faster read mapper `minimap2 <https://github.com/lh3/minimap2>`_.

Output
------

SVIM's main output file called `variants.vcf` is placed into the given working directory.

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
