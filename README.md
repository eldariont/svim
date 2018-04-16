## SVIM - Structural variant identification using long reads

SVIM (pronounced *SWIM*) is a structural variant caller for long reads. It is able to detect and classify six different classes of structural variants.  Unlike existing methods, SVIM integrates information from across the genome to precisely distinguish similar events, such as duplications and cut&paste insertions. In our experiments on simulated and real PacBio data, SVIM reached consistently better results than competing methods, particularly on low-coverage datasets. Furthermore, it is unique in its capability of extracting both the genomic origin and destination of insertions and duplications.

### Background on Structural Variants and Long Reads

Structural variants (SVs) are typically defined as genomic variants larger than 50bps (e.g. deletions, duplications, inversions). Studies have shown that they affect more bases in any given genome than SNPs and small Indels taken together. Consequently, they have a large impact on genes and regulatory regions. This is reflected in the large number of genetic diseases that are caused by SVs.

Common sequencing technologies by providers such as Illumina generate short reads with high accuracy. However, they exhibit weaknesses in repeat and low-complexity regions. This negatively affects SV detection because SVs are associated to such regions. Single molecule long-read sequencing technologies from Pacific Biotechnologies and Oxford Nanopore produce reads with error rates of up to 15% but with lengths of several kb. The high read lengths enable them to cover entire repeats and SVs which facilitates SV detection.

### Input

SVIM analyzes long reads contained in a FASTA file. It was tested on PacBio data only but might work with Nanopore reads as well. Alternatively, SVIM can analyze an alignment file in BAM format. It works best for alignment files produces by [NGM-LR](https://github.com/philres/ngmlr "NGM-LR repository").  

### Output

<img src="https://raw.githubusercontent.com/eldariont/svim/master/docs/SVclasses.png" align="right" width="400px">

SVIM distinguishes six different SV classes: deletions, inversions, cut&paste insertions, novel insertions, interspersed and tandem duplications. For each of these SV classes, it produces a BED file with the SV coordinates. Additonally, a VCF file is produced containing all found SVs.


### Contact
If you experience problems or have suggestions please create an issue or a pull request or contact heller_d@molgen.mpg.de.

### License

The project is licensed under the GNU General Public License.
