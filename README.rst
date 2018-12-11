SVIM - Structural variant identification using long reads
=========================================================

.. image:: https://badge.fury.io/py/svim.svg
    :target: https://badge.fury.io/py/svim

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
    :target: http://bioconda.github.io

SVIM (pronounced *SWIM*) is a structural variant caller for long reads. It is able to detect and classify five different classes of structural variants.  Unlike existing methods, SVIM integrates information from across the genome to precisely distinguish similar events, such as tandem and interspersed duplications and novel element insertions. In our experiments on simulated data and real datasets from PacBio and Nanopore sequencing machines, SVIM reached consistently better results than competing methods. Furthermore, it is unique in its capability of extracting both the genomic origin and destination of duplications.

Background on Structural Variants and Long Reads
------------------------------------------------

.. image:: https://raw.githubusercontent.com/eldariont/svim/master/docs/SVclasses.png
    :align: center

Structural variants (SVs) are typically defined as genomic variants larger than 50bps (e.g. deletions, duplications, inversions). Studies have shown that they affect more bases in any given genome than SNPs and small Indels taken together. Consequently, they have a large impact on genes and regulatory regions. This is reflected in the large number of genetic diseases that are caused by SVs.

Common sequencing technologies by providers such as Illumina generate short reads with high accuracy. However, they exhibit weaknesses in repeat and low-complexity regions. This negatively affects SV detection because SVs are associated to such regions. Single molecule long-read sequencing technologies from Pacific Biotechnologies and Oxford Nanopore produce reads with error rates of up to 15% but with lengths of several kb. The high read lengths enable them to cover entire repeats and SVs which facilitates SV detection.

Installation
------------

.. code-block:: bash

    #Install via conda: easiest option, installs all dependencies including read alignment dependencies
    conda install --channel bioconda svim

    #Install via pip (requires Python 3.6.*): installs all dependencies except those necessary for read alignment (ngmlr, minimap2, samtools)
    pip3 install svim

    #Install from github (requires Python 3.6.*): installs all dependencies except those necessary for read alignment (ngmlr, minimap2, samtools)
    git clone https://github.com/eldariont/svim.git
    cd svim
    pip3 install .


Input
-----

SVIM analyzes long reads given as a FASTA/FASTQ file (uncompressed or gzipped) or a file list. Alternatively, it can analyze an alignment file in BAM format. SVIM was tested on both PacBio and Nanopore data. It works best for alignment files produced by `NGMLR <https://github.com/philres/ngmlr>`_ but also supports the faster read mapper `minimap2 <https://github.com/lh3/minimap2>`_.

Output
------

SVIM distinguishes five different SV classes (see above schema): deletions, inversions, tandem and interspersed duplications and novel insertions. Additionally, SVIM indicates for detected interspersed duplications whether the genomic origin location seems to be deleted in at least one haplotype (indicating a cut&paste insertion) or not (indicating a canonic interspersed duplication). For each of these SV classes, SVIM produces a BED file with the SV coordinates. Additionally, a VCF file is produced containing all found SVs.

Usage
----------------------

Please see our `wiki <https://github.com/eldariont/svim/wiki>`_.

Contact
-------

If you experience problems or have suggestions please create an issue or a pull request or contact heller_d@molgen.mpg.de.

License
-------

The project is licensed under the GNU General Public License.
