#####################################################################################################################
# Written by Keren Carss, 20170816
# Updated by Olga Shamardina
#
# Takes raw (melted) manta and canvas files, and bedintersect output (wao option) showing overlap of both,
# and makes tidy, merged, super table, minimum x% reciprocal overlap, no other filtering.
# Includes deletions called by canvas only, manta only, and both.
#
#####################################################################################################################

# Prepare
args <- commandArgs(TRUE)
SAMPLEID <- args[1]

source("config.sh")

library(plyr)

options(scipen=999)

# Set filenames
canvas.file <- paste0(OUTPUT, "/working_files/canvas_melted_dels/", SAMPLEID, ".CNV.vcf_melt.dels.txt")
manta.file <- paste0(OUTPUT, "/working_files/manta_melted_dels/", SAMPLEID, ".SV.vcf_melt.dels.txt")
overlap.file <- paste0(OUTPUT, "/working_files/manta_canvas_overlap_dels/", SAMPLEID, ".overlap.txt")

output.file <- paste0(OUTPUT, "/working_files/dels_all/", SAMPLEID, ".dels.txt")
output.bed <- paste0(OUTPUT, "/working_files/dels_all/", SAMPLEID, ".dels.bed")

# Read tables
canvas <- read.table(canvas.file, header=T, sep="\t", quote="")
manta <- read.table(manta.file, header=T, sep="\t", quote="")
overlap <- read.table(overlap.file)


# Calculate proportions covered
names(overlap) <- c("CHROM_MANTA", "START_MANTA", "END_MANTA", "CHROM_CANVAS", "START_CANVAS", "END_CANVAS", "OVERLAP_BPS")
overlap["MANTA_DEL_SIZE"] <- overlap$END_MANTA-overlap$START_MANTA
overlap["CANVAS_DEL_SIZE"] <- overlap$END_CANVAS-overlap$START_CANVAS
overlap["MANTA_PROP_COVERED"] <- overlap$OVERLAP_BPS/overlap$MANTA_DEL_SIZE
overlap["CANVAS_PROP_COVERED"] <- overlap$OVERLAP_BPS/overlap$CANVAS_DEL_SIZE
overlap["MANTA_ID"] <- paste(overlap$CHROM_MANTA, overlap$START_MANTA, overlap$END_MANTA)
overlap["CANVAS_ID"] <- paste(overlap$CHROM_CANVAS, overlap$START_CANVAS, overlap$END_CANVAS)

# Identify those with reciprocal overlap >70%
# These will be treated as one
# These should not need any further filering
overlap["CALLED_BY_BOTH_0.7"] <- ifelse(overlap$MANTA_PROP_COVERED>=0.7 & overlap$CANVAS_PROP_COVERED>=0.7, T, F)
overlap <- subset(overlap, CALLED_BY_BOTH_0.7==T)

# Bedtools duplicates manta call if more than one canvas call overlaps with it.
# Here remove dups, keeping one with highest canvas overlap
# AND NOT VICE VERSA (for this purpose, if a canvas call is duplicated coz it overlaps two manta calls, I want to keep them both, but add flag)
overlap <- overlap[order(overlap$MANTA_ID, -overlap$MANTA_PROP_COVERED),]
overlap["MANTA_DUP"] <- duplicated(overlap$MANTA_ID)
overlap <- overlap[overlap$MANTA_DUP==F, ]
overlap$MANTA_DUP <- NULL

overlap["CANVAS_CALL_DUP"] <- overlap$CANVAS_ID %in% unique(overlap$CANVAS_ID[duplicated(overlap$CANVAS_ID)])
overlap$CANVAS_CALL_DUP <- ifelse(overlap$CHROM_CANVAS==".", NA, overlap$CANVAS_CALL_DUP)
overlap <- overlap[,c(10:13,15)]


# Tidy canvas and manta dataframes
canvas$SAMPLE <- NULL
names(canvas) <- c("CANVAS_GT", "CANVAS_RC", "CANVAS_BC", "CANVAS_CN", "CANVAS_FILTER", "CANVAS_CHROM", "CANVAS_START", "CANVAS_QUAL", "CANVAS_END")
canvas["CANVAS_ID"] <- paste(canvas$CANVAS_CHROM, canvas$CANVAS_START-1, canvas$CANVAS_END)

manta$SAMPLE <- NULL
names(manta) <- c("MANTA_GT", "MANTA_GQ", "MANTA_PR", "MANTA_SR", "MANTA_FILTER", "MANTA_CHROM", "MANTA_START", "MANTA_REF", "MANTA_ALT", "MANTA_QUAL", "MANTA_IMPRECISE", "MANTA_END", "MANTA_CIPOS", "MANTA_CIEND", "MANTA_CIGAR", "MANTA_MATEID", "MANTA_HOMLEN", "MANTA_HOMSEQ", "MANTA_SVINSLEN", "MANTA_SVINSSEQ", "MANTA_PAIR_COUNT", "MANTA_UPSTREAM_PAIR_COUNT", "MANTA_DOWNSTREAM_PAIR_COUNT", "MANTA_JUNC_QUAL", "MANTA_CANVAS_ILMN")
manta["MANTA_ID"] <- paste(manta$MANTA_CHROM, manta$MANTA_START-1, manta$MANTA_END)


# Join three data frames together
final <- join(overlap, canvas, by="CANVAS_ID", type="full")
final <- join(final, manta, by="MANTA_ID", type="full")
final["variant_caller"] <- ifelse(is.na(final$MANTA_CHROM), "Canvas", ifelse(is.na(final$CANVAS_CHROM), "Manta", "Manta and Canvas"))
final["ILMN_ID"] <- SAMPLEID
final["MANTA_DEL_SIZE"] <- final$MANTA_END-final$MANTA_START
final["CANVAS_DEL_SIZE"] <- final$CANVAS_END-final$CANVAS_START
final$CHROM <- ifelse(final$variant_caller=="Canvas", as.character(final$CANVAS_CHROM), as.character(final$MANTA_CHROM))
final$START <- ifelse(final$variant_caller=="Canvas", final$CANVAS_START, final$MANTA_START)
final$END <- ifelse(final$variant_caller=="Canvas", final$CANVAS_END, final$MANTA_END)


# Tidy
final <- final[,c("MANTA_CHROM", "MANTA_START", "MANTA_END", "MANTA_DEL_SIZE", "CANVAS_CHROM", "CANVAS_START", "CANVAS_END", "CANVAS_DEL_SIZE", "CHROM", "START", "END", "variant_caller", "ILMN_ID", "MANTA_GT", "MANTA_GQ", "MANTA_FILTER", "MANTA_QUAL", "MANTA_PR", "MANTA_SR", "MANTA_REF", "MANTA_ALT", "MANTA_IMPRECISE", "MANTA_CIPOS", "MANTA_CIEND", "MANTA_CIGAR", "MANTA_MATEID", "MANTA_HOMLEN", "MANTA_HOMSEQ", "MANTA_SVINSLEN", "MANTA_SVINSSEQ", "MANTA_PAIR_COUNT", "MANTA_UPSTREAM_PAIR_COUNT", "MANTA_DOWNSTREAM_PAIR_COUNT", "MANTA_JUNC_QUAL", "CANVAS_GT", "CANVAS_RC", "CANVAS_BC", "CANVAS_CN", "CANVAS_FILTER", "CANVAS_QUAL", "CANVAS_CALL_DUP", "MANTA_PROP_COVERED", "CANVAS_PROP_COVERED", "MANTA_CANVAS_ILMN")]
final <- final[order(final$CHROM, final$START, final$END),]

# replace "." with NA for numerical/empty columns:
for (column in c("MANTA_MATEID", "MANTA_SVINSLEN", "MANTA_JUNC_QUAL", "MANTA_HOMLEN")) {
    final[, column][final[, column]=="."] <- NA
}

# replace "." with "" for charater columns:
final <- data.frame(lapply(final, function(x) if (is.factor(x)) as.character(x) else {x}), stringsAsFactors=FALSE)
final[final == "."] <- ""

# replace "1" with "True" for flag columns:
final$MANTA_IMPRECISE <- ifelse(final$MANTA_IMPRECISE==1, "True", final$MANTA_IMPRECISE)
final$MANTA_CANVAS_ILMN <- ifelse(final$MANTA_CANVAS_ILMN==1, "True", final$MANTA_CANVAS_ILMN)


# Make bed file
bed.file <- final[,c(9:11)]
bed.file$START <- bed.file$START-1


# Write
write.table(final, output.file, col.names=T, row.names=F, sep="\t", quote=F)
write.table(bed.file, output.bed, col.names=F, row.names=F, sep="\t", quote=F)
