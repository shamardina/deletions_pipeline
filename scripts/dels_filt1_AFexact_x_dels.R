#####################################################################################################################
# Written by Keren Carss, 20171121
# Updated by Olga Shamardina
#
# From merged deletion table, calculate exact AFs and write them to a table
#
#####################################################################################################################

# Prepare
source("config.sh")

options(scipen=999)
library(data.table)
library(plyr)

# File names
input.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chrX.txt")
output.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chrX.AFexact.txt")

# Populations and chemistry lists:
population <- c("African" = "AFR", "East-Asian" = "EAS", "European" = "EUR", "Finnish-European" = "FIN", "Other" = "OTH", "South-Asian" = "SAS")
chem <- c(100, 125, 150)

# Read tables
dels <- fread(input.file, header=T, sep="\t", quote="")
dels <- data.frame(dels)
samples <- read.table(MANIFEST, header=T, sep="\t", quote="")

# Prepare data frame
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

# Separate males and females ("M" in terms of number of the X chromosomes here)
delsm <- subset(dels, SEX_X=="M")
delsf <- subset(dels, SEX_X=="F")

# Calculate AFs for Females (correcting for homozygosity, just like autosomes) and males (no correction), and put together
delsf$hom <- ifelse(delsf$MANTA_GT=="1/1" | delsf$CANVAS_GT=="1/1", T,F)
delsf[is.na(delsf$hom),"hom"] <- F
delsf_hom_double <- rbind(subset(delsf, hom==T), delsf)
varcountf <- count(delsf_hom_double$VAR)
names(varcountf) <- c("VAR", "AC_prefilt_all_f")

varcountm <- count(delsm$VAR)
names(varcountm) <- c("VAR", "AC_prefilt_all_m")

varcount <- data.frame(unique(dels$VAR))
names(varcount) <- c("VAR")
varcount <- join(varcount, varcountf, by="VAR")
varcount <- join(varcount, varcountm, by="VAR")

varcount[is.na(varcount$AC_prefilt_all_f),"AC_prefilt_all_f"] <- 0
varcount[is.na(varcount$AC_prefilt_all_m),"AC_prefilt_all_m"] <- 0

varcount$AC_prefilt_all <- varcount$AC_prefilt_all_f + varcount$AC_prefilt_all_m
varcount$AC_prefilt_all_f <- NULL
varcount$AC_prefilt_all_m <- NULL

# Add AN and AF
varcount$AN_all <- (2*(length(unique(delsf$BRIDGE_ID)))) + ( length(unique(delsm$BRIDGE_ID)) )
varcount$AF_prefilt_all <- varcount$AC_prefilt_all/varcount$AN_all

# Repeat for all populations
for (pop in names(population)) {
    AC <- paste0("AC_prefilt_", population[pop])
    ACf <- paste0("AC_prefilt_", population[pop], "_f")
    ACm <- paste0("AC_prefilt_", population[pop], "_m")
    AN <- paste0("AN_", population[pop])
    AF <- paste0("AF_prefilt_", population[pop])

    samp_sub <- subset(samples, ethnicity==pop)$BRIDGE_ID
    dels_sub <- dels[which(dels$BRIDGE_ID %in% samp_sub),]

    delsm_sub <- subset(dels_sub, SEX_X=="M")
    delsf_sub <- subset(dels_sub, SEX_X=="F")

    delsf_sub$hom <- ifelse(delsf_sub$MANTA_GT=="1/1" | delsf_sub$CANVAS_GT=="1/1", T,F)
    delsf_sub[is.na(delsf_sub$hom),"hom"] <- ""
    delsf_sub_hom_double <- rbind(subset(delsf_sub, hom==T), delsf_sub)
    varcountf_tmp <- count(delsf_sub_hom_double$VAR)
    names(varcountf_tmp) <- c("VAR", ACf)

    varcountm_tmp <- count(delsm_sub$VAR)
    names(varcountm_tmp) <- c("VAR", ACm)

    varcount_tmp <- data.frame(unique(dels_sub$VAR))
    names(varcount_tmp) <- c("VAR")
    varcount_tmp <- join(varcount_tmp, varcountf_tmp, by="VAR")
    varcount_tmp <- join(varcount_tmp, varcountm_tmp, by="VAR")

    varcount_tmp[is.na(varcount_tmp[[ACf]]),ACf] <- 0
    varcount_tmp[is.na(varcount_tmp[[ACm]]),ACm] <- 0

    varcount_tmp[[AC]] <- varcount_tmp[[ACf]] + varcount_tmp[[ACm]]

    varcount_tmp[[ACf]] <- NULL
    varcount_tmp[[ACm]] <- NULL

    varcount <- join(varcount, varcount_tmp, by="VAR")
    varcount[is.na(varcount[[AC]]),AC] <- 0

    varcount[[AN]] <- (2*(length(unique(delsf_sub$BRIDGE_ID)))) + ( length(unique(delsm_sub$BRIDGE_ID)) )
    varcount[[AF]] <- round(varcount[[AC]]/varcount[[AN]], 7)
}

# Repeat for all chemistries
for (RL in chem) {
    AC <- paste0("AC_prefilt_", RL, "bp")
    ACf <- paste0("AC_prefilt_", RL, "bp_f")
    ACm <- paste0("AC_prefilt_", RL, "bp_m")
    AN <- paste0("AN_", RL, "bp")
    AF <- paste0("AF_prefilt_", RL, "bp")

    samp_sub <- subset(samples, CHEM==RL)$BRIDGE_ID
    dels_sub <- dels[which(dels$BRIDGE_ID %in% samp_sub),]

    delsm_sub <- subset(dels_sub, SEX_X=="M")
    delsf_sub <- subset(dels_sub, SEX_X=="F")

    delsf_sub$hom <- ifelse(delsf_sub$MANTA_GT=="1/1" | delsf_sub$CANVAS_GT=="1/1", T,F)
    delsf_sub[is.na(delsf_sub$hom),"hom"] <- ""
    delsf_sub_hom_double <- rbind(subset(delsf_sub, hom==T), delsf_sub)
    varcountf_tmp <- count(delsf_sub_hom_double$VAR)
    names(varcountf_tmp) <- c("VAR", ACf)

    varcountm_tmp <- count(delsm_sub$VAR)
    names(varcountm_tmp) <- c("VAR", ACm)

    varcount_tmp <- data.frame(unique(dels_sub$VAR))
    names(varcount_tmp) <- c("VAR")
    varcount_tmp <- join(varcount_tmp, varcountf_tmp, by="VAR")
    varcount_tmp <- join(varcount_tmp, varcountm_tmp, by="VAR")

    varcount_tmp[is.na(varcount_tmp[[ACf]]),ACf] <- 0
    varcount_tmp[is.na(varcount_tmp[[ACm]]),ACm] <- 0

    varcount_tmp[[AC]] <- varcount_tmp[[ACf]] + varcount_tmp[[ACm]]
    varcount_tmp[[ACf]] <- NULL
    varcount_tmp[[ACm]] <- NULL

    varcount <- join(varcount, varcount_tmp, by="VAR")
    varcount[is.na(varcount[[AC]]),AC] <- 0

    varcount[[AN]] <- (2*(length(unique(delsf_sub$BRIDGE_ID)))) + ( length(unique(delsm_sub$BRIDGE_ID)) )
    varcount[[AF]] <- round(varcount[[AC]]/varcount[[AN]], 7)
}

# Write
write.table(varcount, output.file, col.names=T, row.names=F, quote=F, sep="\t")
