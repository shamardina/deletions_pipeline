#####################################################################################################################
# Written by Keren Carss, 20180319
# Updated by Alba Sanchis, Olga Shamardina
#
# Adds prefilt AFs and overlaps.
#
############################################################################


# Prepare
args <- commandArgs(TRUE)
chromosome <- args[1]
dataset.type <- args[2]
if (! dataset.type %in% c("strict", "lenient")) {
    cat("ERROR: unknown dataset type\n")
    quit("no", 1)
}
dataset.type.short <- c("strict" = "str", "lenient" = "lnt")

options(scipen=999)

source("config.sh")

library(plyr)
library(data.table)


# File names
dels.file <- ifelse(dataset.type == "strict",
                    paste0(OUTPUT, "/working_files/dels_filt2_str_by_chr/dels.filt2.str.chr", chromosome, ".txt"),
                    paste0(OUTPUT, "/working_files/dels_filt2_lnt_by_chr/dels.filt2.lnt.chr", chromosome, ".txt"))

AFexpre.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", chromosome, ".AFexact.txt")
ovpre.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", chromosome, ".samp.overlap.internal.txt")
ovpre.AFR.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", chromosome, ".AFR.overlap.internal.txt")
ovpre.EAS.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", chromosome, ".EAS.overlap.internal.txt")
ovpre.EUR.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", chromosome, ".EUR.overlap.internal.txt")
ovpre.FIN.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", chromosome, ".FIN.overlap.internal.txt")
ovpre.OTH.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", chromosome, ".OTH.overlap.internal.txt")
ovpre.SAS.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", chromosome, ".SAS.overlap.internal.txt")
ovpre.100bp.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", chromosome, ".100bp.overlap.internal.txt")
ovpre.125bp.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", chromosome, ".125bp.overlap.internal.txt")
ovpre.150bp.file <- paste0(OUTPUT, "/working_files/dels_filt1_by_chr/dels.filt1.chr", chromosome, ".150bp.overlap.internal.txt")

output.file <- paste0(OUTPUT, "/dataset/", dataset.type, "/dels.filt2.", dataset.type.short[dataset.type], ".chr", chromosome, ".final.txt")
output.gel <- paste0(OUTPUT, "/dataset/", dataset.type, "_nonGEL/dels.filt2.", dataset.type.short[dataset.type], ".chr", chromosome, ".final.nonGEL.txt")


# Read tables
dels <- read.table(dels.file, header=T, sep="\t", quote="")
AFexpre <- read.table(AFexpre.file, header=T, sep="\t", quote="")
ovpre <- read.table(ovpre.file, header=F, sep="\t", quote="")
ovpre.AFR <- read.table(ovpre.AFR.file, header=F, sep="\t", quote="")
ovpre.EAS <- read.table(ovpre.EAS.file, header=F, sep="\t", quote="")
ovpre.EUR <- read.table(ovpre.EUR.file, header=F, sep="\t", quote="")
ovpre.FIN <- read.table(ovpre.FIN.file, header=F, sep="\t", quote="")
ovpre.OTH <- read.table(ovpre.OTH.file, header=F, sep="\t", quote="")
ovpre.SAS <- read.table(ovpre.SAS.file, header=F, sep="\t", quote="")
ovpre.100bp <- read.table(ovpre.100bp.file, header=F, sep="\t", quote="")
ovpre.125bp <- read.table(ovpre.125bp.file, header=F, sep="\t", quote="")
ovpre.150bp <- read.table(ovpre.150bp.file, header=F, sep="\t", quote="")


# Prepare dels table
dels$start_bed <- dels$START-1
dels[is.na(dels$MANTA_SVINSSEQ),"MANTA_SVINSSEQ"] <- ""
dels$VAR <- paste(dels$CHROM, dels$START, dels$END, dels$MANTA_SVINSSEQ)
dels$VARbed <- paste(dels$CHROM, dels$start_bed, dels$END)


# Prepare overlap samp tables
ovpre$V4 <- NULL
ovpre <- unique(ovpre)


# Annotate with pre-filtering exact allele frequency data
dels <- join(dels, AFexpre, by="VAR")


# Annotate with pre-filtering overlap data
ovpre$VARbed <- paste(ovpre$V1, ovpre$V2, ovpre$V3)
ovpre <- ovpre[,4:5]
names(ovpre) <- c("overlap_internal_prefilt", "VARbed")
dels <- join(dels, ovpre, by="VARbed")

ovpre.AFR$VAR <- paste(ovpre.AFR$V1, ovpre.AFR$V2, ovpre.AFR$V3)
ovpre.AFR <- unique(ovpre.AFR[,4:5])
names(ovpre.AFR) <- c("overlap_internal_prefilt_AFR", "VARbed")
dels <- join(dels, ovpre.AFR, by="VARbed")

ovpre.EAS$VAR <- paste(ovpre.EAS$V1, ovpre.EAS$V2, ovpre.EAS$V3)
ovpre.EAS <- unique(ovpre.EAS[,4:5])
names(ovpre.EAS) <- c("overlap_internal_prefilt_EAS", "VARbed")
dels <- join(dels, ovpre.EAS, by="VARbed")

ovpre.EUR$VAR <- paste(ovpre.EUR$V1, ovpre.EUR$V2, ovpre.EUR$V3)
ovpre.EUR <- unique(ovpre.EUR[,4:5])
names(ovpre.EUR) <- c("overlap_internal_prefilt_EUR", "VARbed")
dels <- join(dels, ovpre.EUR, by="VARbed")

ovpre.FIN$VAR <- paste(ovpre.FIN$V1, ovpre.FIN$V2, ovpre.FIN$V3)
ovpre.FIN <- unique(ovpre.FIN[,4:5])
names(ovpre.FIN) <- c("overlap_internal_prefilt_FIN", "VARbed")
dels <- join(dels, ovpre.FIN, by="VARbed")

ovpre.OTH$VAR <- paste(ovpre.OTH$V1, ovpre.OTH$V2, ovpre.OTH$V3)
ovpre.OTH <- unique(ovpre.OTH[,4:5])
names(ovpre.OTH) <- c("overlap_internal_prefilt_OTH", "VARbed")
dels <- join(dels, ovpre.OTH, by="VARbed")

ovpre.SAS$VAR <- paste(ovpre.SAS$V1, ovpre.SAS$V2, ovpre.SAS$V3)
ovpre.SAS <- unique(ovpre.SAS[,4:5])
names(ovpre.SAS) <- c("overlap_internal_prefilt_SAS", "VARbed")
dels <- join(dels, ovpre.SAS, by="VARbed")

ovpre.100bp$VAR <- paste(ovpre.100bp$V1, ovpre.100bp$V2, ovpre.100bp$V3)
ovpre.100bp <- unique(ovpre.100bp[,4:5])
names(ovpre.100bp) <- c("overlap_internal_prefilt_100bp", "VARbed")
dels <- join(dels, ovpre.100bp, by="VARbed")

ovpre.125bp$VAR <- paste(ovpre.125bp$V1, ovpre.125bp$V2, ovpre.125bp$V3)
ovpre.125bp <- unique(ovpre.125bp[,4:5])
names(ovpre.125bp) <- c("overlap_internal_prefilt_125bp", "VARbed")
dels <- join(dels, ovpre.125bp, by="VARbed")

ovpre.150bp$VAR <- paste(ovpre.150bp$V1, ovpre.150bp$V2, ovpre.150bp$V3)
ovpre.150bp <- unique(ovpre.150bp[,4:5])
names(ovpre.150bp) <- c("overlap_internal_prefilt_150bp", "VARbed")
dels <- join(dels, ovpre.150bp, by="VARbed")


# Tidy
dels$VARbed <- NULL
dels$VAR <- NULL
dels$start_bed <- NULL

dels <- dels[order(dels$START, dels$END, dels$BRIDGE_ID),]


# Write
write.table(dels, output.file, col.names=T, row.names=F, sep="\t", quote=F)


# GEL subset
dels_gel <- subset(dels, PROJECT!="GEL")
write.table(dels_gel, output.gel, col.names=T, row.names=F, sep="\t", quote=F)
