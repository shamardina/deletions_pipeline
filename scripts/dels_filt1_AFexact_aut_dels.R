#####################################################################################################################
# Written by Keren Carss, 20171121
# Updated by Olga Shamardina
#
# From merged deletion table, calculate exact AFs and write them to a table
#
#####################################################################################################################

# Prepare
args <- commandArgs(TRUE)
chromosome <- args[1]

source("config.sh")

options(scipen=999)
library(plyr)
library(data.table)

# File names
input.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", chromosome, ".txt")
output.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", chromosome, ".AFexact.txt")

# Populations and chemistry lists:
population <- c("African" = "AFR", "East-Asian" = "EAS", "European" = "EUR", "Finnish-European" = "FIN", "Other" = "OTH", "South-Asian" = "SAS")
chem <- c(100, 125, 150)

# Read tables
dels <- fread(input.file, header=T, sep="\t", quote="")
dels <- data.frame(dels)
samples <- read.table(MANIFEST, header=T, sep="\t", quote="")

# Prepare data frame
dels$hom <- ifelse(dels$MANTA_GT=="1/1" | dels$CANVAS_GT=="1/1", T,F)
dels[is.na(dels$hom),"hom"] <- F
dels[is.na(dels$MANTA_SVINSSEQ),"MANTA_SVINSSEQ"] <- ""
dels$VAR <- paste(dels$CHROM, dels$START, dels$END, dels$MANTA_SVINSSEQ)

# Remove dups
# rare bug
dels <- unique(dels)
dels$VAR2 <- paste(dels$BRIDGE_ID, dels$VAR)
dels <- dels[order(dels$VAR, dels$VAR2, -dels$MANTA_QUAL),]
dels$duplicated <- duplicated(dels$VAR2)
dels <- subset(dels, duplicated==F)
dels$VAR2 <- NULL
dels$duplicated <- NULL

# Count vars, correcting for homozygosity, add info to varcount
dels_hom_double <- rbind(subset(dels, hom==T), dels)
varcount <- count(dels_hom_double$VAR)
names(varcount) <- c("VAR", "AC_prefilt_all")
varcount$AN_all <- 2*(length(unique(dels$BRIDGE_ID))) # 13037*2
varcount$AF_prefilt_all <- round(varcount$AC_prefilt_all/varcount$AN_all, 7)

# Repeat for all populations
for (pop in names(population)) {
    AC <- paste0("AC_prefilt_", population[pop])
    AN <- paste0("AN_", population[pop])
    AF <- paste0("AF_prefilt_", population[pop])
    samp_sub <- subset(samples, ethnicity==pop)$BRIDGE_ID
    dels_sub <- dels[which(dels$BRIDGE_ID %in% samp_sub),]
    dels_sub_hom_double <- rbind(subset(dels_sub, hom==T), dels_sub)
    varcount_tmp <- count(dels_sub_hom_double$VAR)
    names(varcount_tmp) <- c("VAR", AC)
    varcount <- join(varcount, varcount_tmp, by="VAR")
    varcount[is.na(varcount[[AC]]), AC] <- 0
    varcount[[AN]] <- length(samp_sub)*2
    varcount[[AF]] <- round(varcount[[AC]]/varcount[[AN]], 7)
}

# Repeat for all chemistries
for (RL in chem) {
    AC <- paste0("AC_prefilt_", RL, "bp")
    AN <- paste0("AN_", RL, "bp")
    AF <- paste0("AF_prefilt_", RL, "bp")
    samp_sub <- subset(samples, CHEM==RL)$BRIDGE_ID
    dels_sub <- dels[which(dels$BRIDGE_ID %in% samp_sub),]
    dels_sub_hom_double <- rbind(subset(dels_sub, hom==T), dels_sub)
    varcount_tmp <- count(dels_sub_hom_double$VAR)
    names(varcount_tmp) <- c("VAR", AC)
    varcount <- join(varcount, varcount_tmp, by="VAR")
    varcount[is.na(varcount[[AC]]), AC] <- 0
    varcount[[AN]] <- length(samp_sub)*2
    varcount[[AF]] <- round(varcount[[AC]]/varcount[[AN]], 7)
}

# Write
write.table(varcount, output.file, col.names=T, row.names=F, quote=F, sep="\t")
