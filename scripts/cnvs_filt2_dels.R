#####################################################################################################################
# Written by Keren Carss, 20180731
# Updated by Alba Sanchis, Olga Shamardina
#
# Script to do filtering (filt2 lenient and strict), and chr splitting of deletions
#
#####################################################################################################################

# Prepare

source("config.sh")

args <- (commandArgs(TRUE))
SAMPLE <- args[1]

library(data.table)

options(scipen=999)


###################################################################
# File names

dels.file <- paste0(OUTPUT, "/working_files/dels_filt1_ann/", SAMPLE, ".dels.filt1.ann1.txt")
output.dir.str <- paste0(OUTPUT, "/working_files/dels_filt2_str/")
output.dir.lnt <- paste0(OUTPUT, "/working_files/dels_filt2_lnt/")
flagged.file <- paste0(OUTPUT, "/working_files/flagging/all_flagged.txt")
manta.qt.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/manta_quantiles.txt")

# Read files

dels <- data.frame(fread(dels.file, header=T, sep="\t", quote=""))
chromosomes <- c(1:22, "X")
flagged_samples <- read.table(flagged.file)
manta_qt <- read.delim(manta.qt.file)


###################################################################
# Lenient filtering

dels <- subset(dels, variant_caller!="Canvas" | CANVAS_QUAL>=10)
dels <- subset(dels, variant_caller!="Manta" | MANTA_FILTER=="PASS")

# Flagged samples

flagged <- dels$BRIDGE_ID[1] %in% flagged_samples$V1
dels$FLAG_SAMPLE_QUAL <- flagged

# Split data frame by chromosome and write files with no headers and sample flagging

dummy <- lapply(chromosomes, FUN=function(x){write.table(subset(dels, CHROM==x), paste(output.dir.lnt, SAMPLE, ".dels.filt2.lnt.chr", x, ".txt", sep=""), col.names=F, row.names=F, sep="\t", quote=F)} )


###################################################################
# Do strict filtering

# All calls

# Remove variants overlapping flag regions
dels <- subset(dels, FLAGREG_MAX_TYPE=="." | FLAGREG_MAX_PROP_CNV<0.7)

# Canvas calls

dels <- subset(dels, variant_caller!="Canvas" | FLAG_SAMPLE_QUAL==F)
dels <- subset(dels, variant_caller!="Canvas" | SNPS_DENSITY_HET_PASS <= 0.5)

# Canvas and manta calls
# No further filtering

# Manta calls

# Remove variants where evidence from reads suggests they are not real (lenient)
dels <- subset(dels, variant_caller!="Manta" | READS_MEAN_MAPQ>=45 | MANTA_GT=="1/1")
#dels <- subset(dels, variant_caller!="Manta" | READS_MEAN_MAPQ>40 | READ_DENSITY<0.5)

dels <- subset(dels, variant_caller!="Manta" | MANTA_IMPRECISE=="" | is.na(MANTA_IMPRECISE) )
dels <- subset(dels, variant_caller!="Manta" | MANTA_DEL_SIZE<=1000000)
dels <- subset(dels, variant_caller!="Manta" | DUKE0_START==F)
dels <- subset(dels, variant_caller!="Manta" | DUKE0_STOP==F)
dels <- subset(dels, variant_caller!="Manta" | SNPS_N_HET_PASS==0)

# Below thresholds are batch specific removing bottom 20%ile of variants, calculated on filt1 dataset
dels <- subset(dels, variant_caller!="Manta" | CHEM!=100 | MANTA_GQ>manta_qt[manta_qt$CHEM==100, "MANTA_GQ"])  # 53
dels <- subset(dels, variant_caller!="Manta" | CHEM!=125 | MANTA_GQ>manta_qt[manta_qt$CHEM==125, "MANTA_GQ"])  # 48
dels <- subset(dels, variant_caller!="Manta" | CHEM!=150 | MANTA_GQ>manta_qt[manta_qt$CHEM==150, "MANTA_GQ"])  # 34

dels <- subset(dels, variant_caller!="Manta" | CHEM!=100 | MANTA_QUAL>manta_qt[manta_qt$CHEM==100,"MANTA_QUAL"])  # 199
dels <- subset(dels, variant_caller!="Manta" | CHEM!=125 | MANTA_QUAL>manta_qt[manta_qt$CHEM==125,"MANTA_QUAL"])  # 184
dels <- subset(dels, variant_caller!="Manta" | CHEM!=150 | MANTA_QUAL>manta_qt[manta_qt$CHEM==150,"MANTA_QUAL"])  # 136


###################################################################
# Split data frame by chr and write files with no headers

dummy <- lapply(chromosomes, FUN=function(x){write.table(subset(dels, CHROM==x), paste(output.dir.str, SAMPLE, ".dels.filt2.str.chr", x, ".txt", sep=""), col.names=F, row.names=F, sep="\t", quote=F)} )

# DONE
#####################################################################################################################
