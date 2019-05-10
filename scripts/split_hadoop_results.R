########################################################################
# Written by Olga Shamardina
#
# Merge Hadoop's results and then split by sample
#
########################################################################

source("config.sh")

options(scipen=999)

library(data.table)

data = data.frame()
for (chrom in c(1:22, "X")) {
    data.chrom <- fread(paste0(HADOOPOUTPUT, "/", chrom, ".csv"), header=TRUE, data.table=FALSE)
    data.chrom$chrom <- chrom
    data <- rbind(data, data.chrom[,c("chrom", "start", "end", "sample", "count")])
}

for (sample in unique(data$sample)) {
    write.table(data[data$sample == sample,], paste0(OUTPUT, "/working_files/dels_filt1_snps/", sample, ".dels.filt1.snps.txt"), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
}
