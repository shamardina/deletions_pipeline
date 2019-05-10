#####################################################################################################################
# Written by Keren Carss, 20180306
# Updated by Alba Sanchis, Olga Shamardina
#
# Script to do annotation, filtering (filt2), and chr splitting of deletions.
#
# Step 1: Annotate dels_filt1 data with flag region info
# Step 2: Annotate dels_filt1 data with read info (number and average MAPQ)
# Step 3: Annotate dels_filt1 data with n het SNPs
# Step 4: Annotate with gene info
# Step 5: Annotate with affected info and relatedness info
# Step 6: Tidy
#
#####################################################################################################################

# Prepare

args <- (commandArgs(TRUE))
SAMPLE <- args[1]

source("config.sh")

suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))

options(scipen=999)

# File names
dels.file <- paste0(OUTPUT, "/working_files/dels_filt1/", SAMPLE, ".dels.filt1.txt")
flags.file <- paste0(OUTPUT, "/working_files/dels_filt1_flag_overlap/", SAMPLE, ".dels.filt1.flag.txt")
readinfo.file <- paste0(OUTPUT, "/working_files/dels_filt1_ann/", SAMPLE, ".dels.filt1.readinfo.txt")
snps.file <- paste0(OUTPUT, "/working_files/dels_filt1_snps/", SAMPLE, ".dels.filt1.snps.txt")
genes.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chrAll_intersect_genes_merged_VAR.bed")
exons.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chrAll_intersect_genes_strict_merged_VAR.bed")
flagship.data.file <- DATATABLE
output.file <- paste0(OUTPUT, "/working_files/dels_filt1_ann/", SAMPLE, ".dels.filt1.ann1.txt")


###################################################################
# Read files

dels <- data.frame(fread(dels.file))
flags <- read.table(flags.file, header=F, sep="\t", quote="", fill=T)
readinfo <- read.table(readinfo.file, header=F, sep=" ", quote="", fill=T)
snps  <- data.frame(fread(snps.file, header=F, sep="\t", quote=""))
flagship.data <- read.table(flagship.data.file, header=T, sep="\t")
genes <- data.frame(fread(genes.file, header=T, sep="\t", quote=""))
exons <- data.frame(fread(exons.file, header=T, sep="\t", quote=""))


###################################################################
# Do some checks

if (nrow(readinfo)!=nrow(dels)) {
  stop(paste(SAMPLE, ": Check that the number of rows in readinfo file equals the number of deletions.", sep=""))
}


###################################################################
# Order dels

dels$tmp <- dels$CHROM
dels$tmp <- as.character(dels$tmp)
dels[which(dels$tmp=="X"), "tmp"] <- 23
dels <- dels[order(as.integer(as.character(dels$tmp)), dels$START, dels$END),]
dels$tmp <- NULL


###################################################################
# Remove duplicates
# rare bug

dels[is.na(dels$MANTA_SVINSSEQ),"MANTA_SVINSSEQ"] <- ""
dels$VAR <- paste(dels$CHROM, dels$START, dels$END, dels$MANTA_SVINSSEQ)

dels$VAR2 <- paste(dels$BRIDGE_ID, dels$VAR)
dels <- dels[order(dels$VAR, dels$VAR2, -dels$MANTA_QUAL),]
dels$duplicated <- duplicated(dels$VAR2)
dels <- subset(dels, duplicated==F)
dels$VAR <- NULL
dels$VAR2 <- NULL
dels$duplicated <- NULL


#####################################################################################################################
# Step 1: Annotates dels_filt1 data with flag region info

# Prepare flags dataframe

flags$V9 <- NULL
names(flags) <- c("CHROM_CNV", "START_CNV", "STOP_CNV", "CHROM_FLAGREG", "START_FLAGREG", "STOP_FLAGREG", "TYPE_FLAGREG", "OVERLAP")
flags <- unique(flags)

flags["SIZE_CNV"] <- flags$STOP_CNV-flags$START_CNV
flags["SIZE_FLAGREG"] <- flags$STOP_FLAGREG-flags$START_FLAGREG
flags["PROP_CNV"] <- as.numeric(as.character(flags$OVERLAP))/as.numeric(as.character(flags$SIZE_CNV))
flags[which(is.na(flags$PROP_CNV)), "PROP_CNV"] <- 0
flags["PROP_FLAGREG"] <- as.numeric(as.character(flags$OVERLAP))/as.numeric(as.character(flags$SIZE_FLAGREG))
flags[which(is.na(flags$PROP_FLAGREG)), "PROP_FLAGREG"] <- 0
flags["VAR"] <- paste(flags$CHROM_CNV, flags$START_CNV, flags$STOP_CNV)
flags <- subset(flags, CHROM_CNV!=0) # Dodgy call in test

# Annotate flags dataframe with the number of FLAGREGs that overlap each deletion

temp1 <- aggregate(flags$PROP_FLAGREG, by=list(flags$VAR), FUN=length)
names(temp1) <- c("VAR", "N_FLAGREG")
flags <- join(flags, temp1, by="VAR")
flags[which(flags$SIZE_FLAGREG==0), "N_FLAGREG"] <- 0


###################################################################
# Duke 0 specific annotation

# Calculate total proportion of deletion that is non unique

flags_duke0 <- subset(flags, TYPE_FLAGREG=="duke_uniq_0")
temp1b <- aggregate(flags_duke0$PROP_CNV, by=list(flags_duke0$VAR), FUN=sum)
names(temp1b) <- c("VAR", "PROP_CNV_DUKE0")
flags_duke0 <- join(flags_duke0, temp1b, by="VAR")

# Calculate whether the 5' breakpoint (START) of each deletion is within a non-unique (DUKE0) sequence
# In order for a breakpoint to be considered to be within non-unique sequence it fulfils the following criteria:
# Predicted breakpoint is within boundaries of non unique region AND
# The deletion breakpoint is fairly deep within the non-unique region (>=10bps within). This also confers a minimum limit on the size of the non unique region that is used for filtering

flags_duke0$START_IN_DUKE0_ind <- ifelse(flags_duke0$START_CNV>=flags_duke0$START_FLAGREG+10 & flags_duke0$START_CNV<=flags_duke0$STOP_FLAGREG-10, 1, 0)
temp2 <- aggregate(flags_duke0$START_IN_DUKE0_ind, by=list(flags_duke0$VAR), FUN=sum)
names(temp2) <- c("VAR", "START_IN_DUKE0")
temp2[,2] <- ifelse(temp2[,2]==0, F, T)
flags_duke0 <- join(flags_duke0, temp2, by="VAR")

# Calculate whether the 3' breakpoint (STOP) of each deletion is within a non-unique (DUKE0) sequence
# Criteria as above

flags_duke0$STOP_IN_DUKE0_ind <- ifelse(flags_duke0$STOP_CNV>=flags_duke0$START_FLAGREG+10 & flags_duke0$STOP_CNV<=flags_duke0$STOP_FLAGREG-10, 1, 0)
temp3 <- aggregate(flags_duke0$STOP_IN_DUKE0_ind, by=list(flags_duke0$VAR), FUN=sum)
names(temp3) <- c("VAR", "STOP_IN_DUKE0")
temp3[,2] <- ifelse(temp3[,2]==0, F, T)
flags_duke0 <- join(flags_duke0, temp3, by="VAR")

flags_duke0 <- unique(flags_duke0[,c(13,15,17,19)])


###################################################################
# If there is more than one FLAGREG that overlaps with a CNV,
# select the one that has the highest coverage of the CNV, or if the same,
# the highest coverage of the FLAGREG.

flags <- flags[order(flags$VAR, -flags$PROP_CNV, -flags$PROP_FLAGREG),]
flags$duplicated <- duplicated(flags$VAR)
flags <- subset(flags, duplicated==F)
flags$duplicated <- NULL

row.names(flags) <- NULL

if (nrow(flags)!=nrow(dels)) {
  stop(paste(SAMPLE, ": Check that the number of dels in flag file equals the number of deletions.", sep=""))
}

# Add duke 0 info

flags <- join(flags, flags_duke0, by="VAR")


###################################################################
# Tidy flags dataframe

flags <- flags[,c(13, 14, 7, 11, 15:17)]
names(flags) <- c("VAR", "FLAGREG_N", "FLAGREG_MAX_TYPE", "FLAGREG_MAX_PROP_CNV", "DUKE0_TOTAL_PROP_CNV", "DUKE0_START", "DUKE0_STOP")

# Prepare dels dataframe

dels["startbed"] <- as.numeric(dels$START)-1
dels$startbed <- ifelse(dels$startbed==-1,0, dels$startbed)
dels["VAR"] <- paste(dels$CHROM, dels$startbed, dels$END)


###################################################################
# Annotate and tidy

dels <- join(dels, flags, by="VAR")
dels[which(is.na(dels$DUKE0_START)), "DUKE0_START"] <- FALSE
dels[which(is.na(dels$DUKE0_STOP)), "DUKE0_STOP"] <- FALSE
dels[which(is.na(dels$DUKE0_TOTAL_PROP_CNV)), "DUKE0_TOTAL_PROP_CNV"] <- 0


#####################################################################################################################
# Step 2: Annotates dels_filt1 data with read info (number and average MAPQ)

names(readinfo) <- c("CHROM", "START", "END", "READS_N", "READS_MEAN_MAPQ")
readinfo$VAR <- paste(readinfo$CHROM, readinfo$START, readinfo$END)
readinfo <- readinfo[,4:6]
readinfo <- unique(readinfo)
dels <- join(dels, readinfo, by="VAR")
dels$size <- dels$END-dels$START
dels$READ_DENSITY <- round((dels$READS_N / dels$size), 2)


#####################################################################################################################
# Step 3: Annotate dels_filt1 data with n het SNPs

names(snps) <- c("CHROM", "START", "END", "SAMPLE", "SNPS_N_HET_PASS")
snps$VAR <- paste(snps$CHROM, snps$START, snps$END)
snps <- snps[,5:6]
snps <- unique(snps)
dels <- join(dels, snps, by="VAR")
dels$SNPS_N_HET_PASS[is.na(dels$SNPS_N_HET_PASS)] <- 0  # "NA" produced by join should be 0
dels$SNPS_DENSITY_HET_PASS <- round((dels$SNPS_N_HET_PASS / (dels$size/1000)), 2)


#####################################################################################################################
# Step 4: Annotate with gene info

genes$in_dels <- genes$VAR %in% dels$VAR
genes <- unique(subset(genes, in_dels==T))
if (nrow(genes)!=nrow(dels)) {
  stop(paste(SAMPLE, ": Check that the number of dels in gene file equals the number of deletions.", sep=""))
}
dels <- join(dels, genes[, c("GENES", "VAR")], by="VAR")

colnames(exons)[which(colnames(exons)=="GENES")] <- "GENES_STRICT"
exons$in_dels <- exons$VAR %in% dels$VAR
exons <- unique(subset(exons, in_dels==T))
if (nrow(exons)!=nrow(dels)) {
  stop(paste(SAMPLE, ": Check that the number of dels in gene strict file equals the number of deletions.", sep=""))
}
dels <- join(dels, exons[, c("GENES_STRICT", "VAR")], by="VAR")


#####################################################################################################################
# Step 5: Annotate with ethnicity, affected info and relatedness info

flagship.data <- flagship.data[,c(1,6,7,9)]
names(flagship.data) <- c("BRIDGE_ID", names(flagship.data)[2:4])
dels <- join(dels, flagship.data, by="BRIDGE_ID")


#####################################################################################################################
# Step 6: Tidy

dels <- dels[, c("CHROM", "START", "END", "variant_caller", "size", "MANTA_CHROM", "MANTA_START", "MANTA_END", "MANTA_DEL_SIZE", "CANVAS_CHROM", "CANVAS_START", "CANVAS_END", "CANVAS_DEL_SIZE", "ILMN_ID", "BRIDGE_ID", "PROJECT", "CHEM", "SEX_X", "ETHNICITY", "UNRELATED", "AFFECTED", "GENES", "GENES_STRICT", "MANTA_GT", "MANTA_GQ", "MANTA_FILTER", "MANTA_QUAL", "MANTA_PR", "MANTA_SR", "MANTA_REF", "MANTA_ALT", "MANTA_IMPRECISE", "MANTA_CIPOS", "MANTA_CIEND", "MANTA_CIGAR", "MANTA_MATEID", "MANTA_HOMLEN", "MANTA_HOMSEQ", "MANTA_SVINSLEN", "MANTA_SVINSSEQ", "MANTA_PAIR_COUNT", "MANTA_UPSTREAM_PAIR_COUNT", "MANTA_DOWNSTREAM_PAIR_COUNT", "MANTA_JUNC_QUAL", "CANVAS_GT", "CANVAS_RC", "CANVAS_BC", "CANVAS_CN", "CANVAS_FILTER", "CANVAS_QUAL", "CANVAS_CALL_DUP", "MANTA_PROP_COVERED", "CANVAS_PROP_COVERED", "FLAGREG_N", "FLAGREG_MAX_TYPE", "FLAGREG_MAX_PROP_CNV", "DUKE0_TOTAL_PROP_CNV", "DUKE0_START", "DUKE0_STOP", "READS_N", "READ_DENSITY", "READS_MEAN_MAPQ", "SNPS_N_HET_PASS", "SNPS_DENSITY_HET_PASS")]
write.table(dels, output.file, col.names=T, row.names=F, sep="\t", quote=F)


#####################################################################################################################
