#####################################################################################################################
# Written by Keren Carss, 20171121
# Updated by Olga Shamardina
#
# Annotate with PARs and sample info
# Do some filtering (remove huge deletions and het M X)
#
#####################################################################################################################

# Prepare
args <- commandArgs(TRUE)
SAMPLEID <- args[1]

source("config.sh")

options(scipen=999)

library(plyr)

# Set input filenames
dels.file <- paste0(OUTPUT, "/working_files/dels_all/", SAMPLEID, ".dels.txt")
PAR.file <- paste0(OUTPUT, "/working_files/dels_all_PAR_overlap/", SAMPLEID, ".dels.PAR.txt")

# Read tables
dels <- read.table(dels.file, header=T, sep="\t", quote="")
manifest <- read.table(MANIFEST, header=T)
PAR <- read.table(PAR.file, sep="\t", quote="")

# Annotate with PAR info
dels$bedstart <- dels$START-1
dels$VAR <- paste(dels$CHROM, dels$bedstart, dels$END)

PAR$VAR <- paste(PAR$V1, PAR$V2, PAR$V3)
PAR <- PAR[,c(4:5)]
names(PAR) <- c("PAR", "VAR")

dels <- join(dels, PAR, by="VAR")

# Annotate with sample info
manifest$SVVCF_PATH <- NULL
dels <- join(dels, manifest, by="ILMN_ID")

# Tidy
dels <- dels[,c(1:11,47,12:13,48:53,14:46)]
dels$bedstart <- NULL
dels$VAR <- NULL

# Remove huge deletions
dels <- subset(dels, MANTA_DEL_SIZE<=50000000 | is.na(MANTA_DEL_SIZE))
dels <- subset(dels, CANVAS_DEL_SIZE<=50000000 | is.na(CANVAS_DEL_SIZE))

# Remove het M X ("M" in terms of number of the X chromosomes here)
dels$HET <- ifelse( ((dels$CANVAS_GT=="0/1" & is.na(dels$MANTA_GT)) | (dels$CANVAS_GT=="0/1" & dels$MANTA_GT=="0/1") | (is.na(dels$CANVAS_GT) & dels$MANTA_GT=="0/1")), T, F)
dels$HET_M_X <- ifelse(dels$SEX_X=="M" & dels$CHROM=="X" & dels$PAR==0 & dels$HET==T, T, F)
dels <- subset(dels, HET_M_X==F)

dels$HET <- NULL
dels$HET_M_X <- NULL

# Make beds
bed.file <- dels[,c(9:11)]
bed.file$START <- bed.file$START-1

# Get BRIDGE_ID
SAMPLEID2 <- as.character(dels[1,16])

# Set output filenames
output.file <- paste0(OUTPUT, "/working_files/dels_filt1/", SAMPLEID2, ".dels.filt1.txt")
output.bed <- paste0(OUTPUT, "/working_files/dels_filt1/", SAMPLEID2, ".dels.filt1.bed")

# Write
write.table(dels, output.file, col.names=T, row.names=F, sep="\t", quote=F)
write.table(bed.file, output.bed, col.names=F, row.names=F, sep="\t", quote=F)
