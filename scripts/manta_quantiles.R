########################################################################
# Written by Olga Shamardina
#
# Batch specific MANTA_GQ and MANTA_QUAL bottom 20%ile, calculated on filt1 dataset
#
########################################################################

library(data.table)

options(scipen=999)

source("config.sh")

dels <- data.frame()

for (chrom in c(1:22, "X")) {
    new <- fread(paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", chrom, ".txt"), data.table=FALSE, select=c("CHROM", "variant_caller", "CHEM", "MANTA_QUAL", "MANTA_GQ"))
    dels <- rbind(dels, new[new$variant_caller!="Canvas",])
}

qual <- aggregate(MANTA_QUAL ~ CHEM, data=dels, FUN=function(x) {quantile(x, probs=0.2)})
gq <- aggregate(MANTA_GQ ~ CHEM, data=dels, FUN=function(x) {quantile(x, probs=0.2)})

write.table(merge(gq, qual), paste0(OUTPUT, "/working_files/dels_filt1_by_chr/manta_quantiles.txt"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
