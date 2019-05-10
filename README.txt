Written by Keren Carss, 19th March 2018
kjn29@cam.ac.uk

Updated by Alba Sanchis, 8th Aug 2018
as2635@cam.ac.uk

Updated by Olga Shamardina, 10th May 2019
os360@cam.ac.uk

——————————————————————————————————————————
Flagship dataset (n = 13,037)

——————————————————————————————————————————
CONTENTS OF THIS README FILE

1) INTRODUCTION
2) STRUCTURAL VARIANT CALLING ALGORITHMS
3) FILES
4) FILTERING STRICT
5) FILTERING LENIENT
6) COLUMN DEFINITIONS
7) WARNINGS


——————————————————————————————————————————
1) INTRODUCTION

The files in this directory contain deletion data for the flagship dataset (n = 13,037) from all WGS10K projects.

Each row is a deletion in an individual (so if two individuals have the same deletion there will be two rows).

All coordinates and annotations are GRCh37.

——————————————————————————————————————————
2) STRUCTURAL VARIANT CALLING ALGORITHMS

Structural variants are called by Illumina using two Illumina algorithms: Canvas and Manta.

Canvas uses read depth information and is able to call copy number losses and copy number gains. It is optimised for variants >10kb. 

Canvas reference: Roller et al., Canvas: versatile and scalable detection of copy number variants. Bioinformatics. 2016 Aug 1;32(15):2375-7. doi: 10.1093/bioinformatics/btw163. 

Manta uses both paired-read fragment size and split read information and is able to call deletions, insertions, tandem duplications, inversions and translocations. For deletions it is optimised for 50bp-10kb. It cannot detect non-tandem repeats/amplifications, large insertions, or small inversions.

Manta reference: Chen et al., Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. Bioinformatics. 2016 Apr 15;32(8):1220-2. doi: 10.1093/bioinformatics/btv710.

For more information see also the Illumina Whole-Genome Sequencing Services User Guide which can be found here:

https://support.illumina.com/downloads/whole_genome_sequencing_service_user_guide.html

——————————————————————————————————————————
3) FILES

Files contain deletions called by Canvas only, by Manta only, and by both (defined as >0.7 reciprocal overlap between Canvas and Manta deletions in the same individual)

- Text files of deletions split per chromosome strict filtering (recommended for association analyses)
- Text files of deletions split per chromosome lenient filtering (recommended for pertinent findings analysis)
- As above, non GEL subsets for analysts without permission to see GEL data.

——————————————————————————————————————————
4) FILTERING STRICT

The following types of deletions have been removed:

>50Mb
Heterozygous deletion on X (not PAR) in a sample with one X chromosome
MT and Y
Called by Canvas only and CANVAS_QUAL<10
Called by Manta only and MANTA_FILTER != “PASS”

>0.7 of the deletion overlaps a flagged region*

Called by Canvas only, SNPS_DENSITY_HET_PASS > 0.5
Called by Canvas only in a sample with excessive number of Canvas calls (excessive is defined as CNV count >5SDs>mean) OR excessive total bps deleted (excessive is defined as >99.7th percentile) with more than 50 Canvas deletions in the lenient dataset (non-uniform coverage leading to excess of false positive calls).

Called by Manta only and MANTA_QUAL below 20th centile, calculated per chemistry
Called by Manta only and MANTA_GQ below 20th centile, calculated per chemistry
Called by Manta only and MANTA_IMPRECISE==TRUE
Called by Manta only and >1Mb
Called by Manta only and DUKE0_START = TRUE or DUKE0_STOP = TRUE
Called by Manta only and SNPS_N_HET_PASS > 0
Called by Manta only and heterozygous and READS_MEAN_MAPQ < 45


* Flagged regions are:
	2 x sets low mappability regions from Encode
	Centromeres
	Telomeres
	IG loci
	HLA loci
	Segmental duplications
	Duke non-unique regions

——————————————————————————————————————————
5) FILTERING LENIENT

The following types of deletions have been removed:

>50Mb
Heterozygous deletion on X (not PAR) in a sample with one X chromosome
MT and Y
Called by Canvas only and CANVAS_QUAL<10
Called by Manta only and MANTA_FILTER != “PASS”

——————————————————————————————————————————
6) COLUMN DEFINITIONS

CHROM				Chromosome
START				Start (if deletion is called by both Manta and Canvas, Manta coordinate is used)
END				End (if deletion is called by both Manta and Canvas, Manta coordinate is used)
variant_caller			Variant caller
size				Size (bps)
MANTA_CHROM			Chromosome (Manta, if applicable)
MANTA_START			Start (Manta, if applicable)
MANTA_END			End (Manta, if applicable)
MANTA_DEL_SIZE			Size (Manta, if applicable)
CANVAS_CHROM			Chromosome (Canvas, if applicable)
CANVAS_START			Start (Canvas, if applicable)
CANVAS_END			End (Canvas, if applicable)
CANVAS_DEL_SIZE			Size (Canvas, if applicable)
ILMN_ID				Illumina sample ID
BRIDGE_ID			BRIDGE sample ID
PROJECT				BRIDGE project
CHEM				Sequencing chemistry (read length)
SEX_X				Likely sex karyotype of individual (inferred from WGS data) in terms of the number of X chromosomes. Thus XXY sample is labelled as “F” and XO as “M”.
FLAG_SAMPLE_QUAL		Quality issue with sample. TRUE if sample had excessive number of Canvas calls OR excessive amount of bps deleted along with more than 50 Canvas deletions.
ETHNICITY			Likely ethnicity of individual (inferred from WGS data)
UNRELATED			Unrelated
AFFECTED			Affected
GENES				Genes within or overlapping predicted deletion (defined as any exon)
GENES_STRICT			Genes within or overlapping predicted deletion (defined as any CDS of canonical transcript of protein-coding, HGNC gene.)
MANTA_GT			Genotype (Manta, if applicable)
MANTA_GQ			Genotype quality (Manta, if applicable)
MANTA_FILTER			Filter field (Manta, if applicable) (MaxDepth=Sample site depth is > 3× the mean chromosome depth near 1 or both variant break-ends; MaxMQ0Frac=For a small variant (< 1000 bases), the fraction of reads with MAPQ=0 around either break-end exceeds 0.4; MGE10kb=Manta DEL or DUP call with length ≥ 10 kb; MinGQ=GQ score is < 20; NoPairSupport=For variants significantly larger than the paired read fragment size, no paired reads support the alternate allele; Ploidy=For DEL and DUP variants, the genotypes of overlapping variants with smaller size are inconsistent with diploid expectation)
MANTA_QUAL			Quality score (Manta, if applicable)
MANTA_PR			Spanning paired-read support for the REF and ALT alleles in the order listed (Manta, if applicable)
MANTA_SR			Split reads for the REF and ALT alleles in the order listed, for reads where P(allele|read) > 0.999 (Manta, if applicable)
MANTA_REF			Reference allele (Manta, if applicable)
MANTA_ALT			Alternate allele (Manta, if applicable)
MANTA_IMPRECISE			Imprecise structural variation. I.e. deletion scored on paired-end read evidence alone and not split reads. (Manta, if applicable)
MANTA_CIPOS			Confidence interval around POS (Manta, if applicable)
MANTA_CIEND			Confidence interval around END (Manta, if applicable)
MANTA_CIGAR			CIGAR alignment for each alternate indel allele (Manta, if applicable)
MANTA_MATEID			ID of mate break-end (Manta, if applicable)
MANTA_HOMLEN			Length of base pair identical microhomology at event breakpoints (Manta, if applicable)
MANTA_HOMSEQ			Sequence of base pair identical microhomology at event breakpoints (Manta, if applicable)
MANTA_SVINSLEN			Length of insertion (Manta, if applicable)
MANTA_SVINSSEQ			Sequence of insertion (Manta, if applicable)
MANTA_PAIR_COUNT		Read pairs supporting this variant where both reads are confidently mapped (Manta, if applicable)
MANTA_UPSTREAM_PAIR_COUNT	Confidently mapped reads supporting this variant at the upstream breakend (mapping may not be confident at downstream breakend) (Manta, if applicable)
MANTA_DOWNSTREAM_PAIR_COUNT	Confidently mapped reads supporting this variant at this downstream breakend (mapping may not be confident at upstream breakend) (Manta, if applicable)
MANTA_JUNC_QUAL			Provides the QUAL value for only the adjacency in question (Manta, if applicable)
CANVAS_GT			Genotype
CANVAS_RC			Mean counts per bin in the region
CANVAS_BC			Number of bins in the region
CANVAS_CN			Copy number
CANVAS_FILTER			Filter field (q10=Quality below 10; CLT10kb=Canvas call with length < 10 kb)
CANVAS_QUAL			Quality score
CANVAS_CALL_DUP			The Canvas call has been duplicated because it overlaps with more than one Manta call
MANTA_PROP_COVERED		Proportion of the Manta call that is covered by the Canvas call
CANVAS_PROP_COVERED		Proportion of the Canvas call that is covered by the Manta call
FLAGREG_N			Number of flagged regions that the deletion overlaps with (see previous section)
FLAGREG_MAX_TYPE		The type of the flagged region with greatest overlap with the deletion
FLAGREG_MAX_PROP_CNV		Overlap of the flagged region with greatest overlap with the deletion
DUKE0_TOTAL_PROP_CNV		Total proportion of the deletion that overlaps Duke non-unique regions
DUKE0_START			Start coordinate of deletion is within Duke non-unique region
DUKE0_STOP			End coordinate of deletion is within Duke non-unique region
READS_N				Number of reads within deletion boundaries
READ_DENSITY			READS_N / size(bp)
READS_MEAN_MAPQ			Mean MAPQ of reads within deletion boundaries
SNPS_N_HET_PASS			Number of PASS (OPR0.99) heterozygous SNPs within deletion boundaries
SNPS_DENSITY_HET_PASS		SNPS_N_HET_PASS / size(kbp)
AC_prefilt			Allele count (chrom, start, end, AND MANTA_SVINSSEQ identical) prior to filtering
AN				Allele number (calculated from sample number)
AF_prefilt			AC_prefilt/AN
AC_prefilt_*			Population-specific/Chemistry-specific AC prefilt
AN_*				Population-specific/Chemistry-specific AN prefilt
AF_prefilt_*			Population-specific/Chemistry-specific AF prefilt
overlap_internal_prefilt	Number of deletions within this dataset that overlap variant (reciprocal overlap of >=80%) prior to filtering. This gives a broader indication of the variability of the region than the AFs.
overlap_internal_prefilt_*	Population-specific/Chemistry-specific overlap

——————————————————————————————————————————
7) WARNINGS

Calling structural variants is challenging, these structural variant files do contain false positives. Estimated overall FDR of strict dataset=3%. Visual inspection of a variant of interest using IGV is recommended in the first instance to help assess whether it is real. To be confident that a structural variant is real, independent experimental validation is usually required.

——————————————————————————————————————————
