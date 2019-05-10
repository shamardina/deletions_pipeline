########################################################################
# Written by Olga Shamardina
#
# Identify which samples need to be flagged
#
########################################################################

source("config.sh")

options(scipen=999)

### Category 1 ###
# Samples with excessive number of canvas calls (due to non-uniform coverage) based on batch_report.csv
# Excessive is defined as CNV count >5SDs>mean calculated on included samples

CNVs.count <- read.delim(MANIFEST)[, c("BRIDGE_ID", "CNVCount")]
M <- mean(CNVs.count$CNV)
S <- sd(CNVs.count$CNV)

flagged1 <- as.character(CNVs.count$BRIDGE_ID[CNVs.count$CNV > M + 5*S])
write(flagged1, paste0(OUTPUT, "/working_files/flagging/flagged1.txt"))

### Category 2 ###
# Based on lenient dataset:
# samples with excessive number of bps deleted (number of deleted base pairs is > 99.7 %ile)
# AND having more than 50 Canvas calls

excessive.counts <- read.delim(paste0(OUTPUT, "/working_files/flagging/deletions_counts.txt"), header=FALSE, col.names=c("ID", "deleted_bp", "dels"))

threshold <- quantile(excessive.counts$deleted_bp, prob=0.997)
print(threshold)
flagged2 <- as.character(excessive.counts[excessive.counts$deleted_bp > threshold & excessive.counts$dels > 50, "ID"])
write(flagged2, paste0(OUTPUT, "/working_files/flagging/flagged2.txt"))

### Write both lists merged ###
write(unique(c(flagged1, flagged2)), paste0(OUTPUT, "/working_files/flagging/all_flagged.txt"))
