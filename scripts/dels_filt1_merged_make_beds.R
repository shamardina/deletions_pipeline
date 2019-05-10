#####################################################################################################################
# Written by Keren Carss, 20171121
# Updated by Olga Shamardina
#
# Reformat dels_merged files to make beds
# Includes ethnicity-specific and chem-specific capability
#
#########################################################################

# Prepare
args <- commandArgs(TRUE)
CHROM <- args[1]

source("config.sh")

options(scipen=999)

library(data.table)

# Filenames
dels.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", CHROM, ".txt")
out.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", CHROM, ".bed")
out.file.type <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", CHROM, ".")

# Populations and chemistry lists:
population <- c("African" = "AFR", "East-Asian" = "EAS", "European" = "EUR", "Finnish-European" = "FIN", "Other" = "OTH", "South-Asian" = "SAS")
chem <- c(100, 125, 150)

# Read
dels <- fread(dels.file, header=T, sep="\t", quote="")
dels <- data.frame(dels)
samples <- read.table(MANIFEST, header=T, sep="\t", quote="")

# Make beds
dels$bed_start <- dels$START-1

bed <- dels[,c("CHROM", "bed_start", "END", "BRIDGE_ID")]
bed <- bed[order(bed$bed_start, bed$END),]
write.table(bed[, colnames(bed) != "BRIDGE_ID"], out.file, col.names=F, row.names=F, quote=F, sep="\t")
write.table(bed, paste0(out.file.type, "samp.bed"), col.names=F, row.names=F, quote=F, sep="\t")

for (pop in names(population)) {
    samp_sub <- subset(samples, ethnicity==pop)$BRIDGE_ID
    dels_pop <- dels[which(dels$BRIDGE_ID %in% samp_sub),]
    bed_pop <- dels_pop[,c("CHROM", "bed_start", "END", "BRIDGE_ID")]
    bed_pop <- bed_pop[order(bed_pop$bed_start, bed_pop$END),]
    write.table(bed_pop[, colnames(bed_pop) != "BRIDGE_ID"], paste0(out.file.type, population[pop], ".bed"), col.names=F, row.names=F, quote=F, sep="\t")
    write.table(bed_pop, paste0(out.file.type, "samp.", population[pop], ".bed"), col.names=F, row.names=F, quote=F, sep="\t")
}

for (RL in chem) {
    dels_bp <- subset(dels, CHEM==RL)
    bed_bp <- dels_bp[,c("CHROM", "bed_start", "END", "BRIDGE_ID")]
    bed_bp <- bed_bp[order(bed_bp$bed_start, bed_bp$END),]
    write.table(bed_bp[, colnames(bed_bp) != "BRIDGE_ID"], paste0(out.file.type, RL, "bp.bed"), col.names=F, row.names=F, quote=F, sep="\t")
    write.table(bed_bp, paste0(out.file.type, "samp.", RL, "bp.bed"), col.names=F, row.names=F, quote=F, sep="\t")
}
