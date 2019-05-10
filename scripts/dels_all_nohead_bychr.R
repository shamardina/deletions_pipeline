#####################################################################################################################
# Written by Keren Carss, 20170829
# Updated by Olga Shamardina
#
# Reformat dels_all files for pre-filt AF calculations
#
#####################################################################################################################

# Prepare
args <- commandArgs(TRUE)
SAMPLEID <- args[1]

source("config.sh")

options(scipen=999)

# Set filenames
dels.file <- paste0(OUTPUT, "/working_files/dels_filt1/", SAMPLEID, ".dels.filt1.txt")

# Read tables
dels <- read.table(dels.file, header=T, sep="\t", quote="")

# Fix sex=F
if (dels$SEX_X[1]==FALSE) dels$SEX_X <- "F"

# Split data frame by chr and write files with no headers
for (chrom in c(1:22, "X")) write.table(subset(dels, CHROM == chrom), paste0(OUTPUT, "/working_files/dels_all_nohead_bychr/", SAMPLEID, ".dels.chr", chrom, ".txt"), col.names=F, row.names=F, sep="\t", quote=F)
