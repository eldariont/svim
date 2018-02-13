## SVIM - Structural variant identification using long reads

SVIM (pronounced *SWIM*) is a structural variant caller for long reads. It combines read-tail mapping, full alignment analysis, and split-read mapping to distinguish five classes of structural variants. SVIM discriminates between similar SV classes such as interspersed duplications and cut&paste insertions and is unique in its capability of extracting both the genomic origin and destination of insertions and duplications.

### Background

Structural variants (SVs) are defined as genomic variants larger than 50 bps. They have been shown to affect more bases in any given genome than SNPs and small indels. Additionally, they have great impact on human phenotype and diversity and have been linked to numerous diseases. Due to their size and association with repeats, they are difficult to detect in shotgun sequencing, especially when based on short reads. Long read, single molecule sequencing technologies like those offered by Pacific Biosciences and Oxford Nanopore produce reads with lengths of several thousand base pairs. Despite their higher error rate, long read sequencing offers many advantages for the detection of structural variants, but available software tools still do not fully exploit the possibilities.

### Input

SVIM analyzes long reads contained in a FASTA file. It was tested on PacBio data only but might work with Nanopore reads as well. Alternatively, SVIM can analyze an alignment file in BAM format. It works best for alignment files produces by [NGM-LR](https://github.com/philres/ngmlr "NGM-LR repository").  

### Output

<img src="https://raw.githubusercontent.com/eldariont/svim/master/docs/SVclasses.png" align="right" width="400px">

SVIM distinguishes five different SV classes: deletions, inversions, cut&paste insertions, interspersed and tandem duplications. For each of these SV classes, it produces a BED file with the SV coordinates. Additonally, a VCF file is produced containing all found SVs.


### Contact
If you experience problems or have suggestions please create issues or pull requests or contact heller_d@molgen.mpg.de.

### License

The project is licensed under the GNU General Public License.
