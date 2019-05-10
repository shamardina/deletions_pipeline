#!/bin/bash

set -e
set -o pipefail
set -u


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
cd ${DIR}
source config.sh
source functions.sh
module load parallel


echo
echo "Step 1: Convert VCFs to table (get only deletions), $(date)"
module load samtools/1.8
mkdir -p ${OUTPUT}/working_files/{canvas_melted_dels,manta_melted_dels}
chromosome_regions=$(echo $(seq 1 22) X | tr ' ' ,)

# Canvas deletions:
Canvas_VCF_header='[%SAMPLE\t%GT\t%RC\t%BC\t%CN]\t%FILTER\t%CHROM\t%POS\t%QUAL\t%END\n'
Canvas_out_header="$(echo $Canvas_VCF_header | sed 's/\\t\|\\n\|\[\|\]//g; s/%//; s/%/\\t/g')"
parallel --halt 1 -j $NJOBS --colsep \\t --verbose "echo -e '$Canvas_out_header' > ${OUTPUT}/working_files/canvas_melted_dels/{2}.CNV.vcf_melt.dels.txt; bcftools query -r $chromosome_regions -i 'SVTYPE==\"CNV\" && CN<2' -f '$Canvas_VCF_header' {1} >> ${OUTPUT}/working_files/canvas_melted_dels/{2}.CNV.vcf_melt.dels.txt" :::: <(paste <(get_samples 1 vcf) <(get_samples 1 iid))

# Manta deletions:
Manta_VCF_header='[%SAMPLE\t%GT\t%GQ\t%PR\t%SR]\t%FILTER\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%IMPRECISE\t%END\t%CIPOS\t%CIEND\t%CIGAR\t%MATEID\t%HOMLEN\t%HOMSEQ\t%SVINSLEN\t%SVINSSEQ\t%PAIR_COUNT\t%UPSTREAM_PAIR_COUNT\t%DOWNSTREAM_PAIR_COUNT\t%JUNCTION_QUAL\t%ColocalizedCanvas\n'
Manta_out_header="$(echo $Manta_VCF_header | sed 's/\\t\|\\n\|\[\|\]//g; s/%//; s/%/\\t/g')"
parallel --halt 1 -j $NJOBS --colsep \\t --verbose "echo -e '$Manta_out_header' > ${OUTPUT}/working_files/manta_melted_dels/{2}.SV.vcf_melt.dels.txt; bcftools query -r $chromosome_regions -i 'SVTYPE==\"DEL\"' -f '$Manta_VCF_header' {1} >> ${OUTPUT}/working_files/manta_melted_dels/{2}.SV.vcf_melt.dels.txt" :::: <(paste <(get_samples 1 vcf) <(get_samples 1 iid))

module unload samtools
update_stage 1
echo "Step 1: finished, sample processing step in ${LIST} updated, $(date)"


echo
echo "Step 2: Calculate manta-canvas overlap, $(date)"
mkdir -p ${OUTPUT}/working_files/manta_canvas_overlap_dels
module load bedtools/2.22.0
parallel --halt 1 -j $NJOBS --verbose "intersectBed -wao -a <(cut -f 7,8,13 ${OUTPUT}/working_files/manta_melted_dels/{}.SV.vcf_melt.dels.txt | sed 1d | vcf2bed) -b <(cut -f 7,8,10 ${OUTPUT}/working_files/canvas_melted_dels/{}.CNV.vcf_melt.dels.txt | sed 1d | vcf2bed) > ${OUTPUT}/working_files/manta_canvas_overlap_dels/{}.overlap.txt" :::: <(get_samples 2 iid)
module unload bedtools
update_stage 2
echo "Step 2: finished, sample processing step in ${LIST} updated, $(date)"


echo
echo "Step 3: Make single dels table (ie re-combine canvas and manta calls), $(date)"
mkdir -p ${OUTPUT}/working_files/dels_all
parallel --halt 1 -j $NJOBS --verbose "Rscript make_single_dels_tab.R {}" :::: <(get_samples 3 iid)
update_stage 3
echo "Step 3: finished, sample processing step in ${LIST} updated, $(date)"


echo
echo "Step 4: Overlap with PARs, $(date)"
mkdir -p ${OUTPUT}/working_files/dels_all_PAR_overlap
module load bedtools/2.22.0
parallel --halt 1 -j $NJOBS --verbose "intersectBed -c -f 0.5 -a ${OUTPUT}/working_files/dels_all/{}.dels.bed -b ${RESOURCES}/PAR_coordinates.bed > ${OUTPUT}/working_files/dels_all_PAR_overlap/{}.dels.PAR.txt" :::: <(get_samples 4 iid)
module unload bedtools
update_stage 4
echo "Step 4: finished, sample processing step in ${LIST} updated, $(date)"


echo
echo "Step 5: Annotation and filtering stage 1 (removes het calls on male X and anything >50Mb only), $(date)"
mkdir -p ${OUTPUT}/working_files/dels_filt1
parallel --halt 1 -j $NJOBS --verbose "Rscript dels_ann_filt1.R {}" :::: <(get_samples 5 iid)
update_stage 5
echo "Step 5: finished, sample processing step in ${LIST} updated, $(date)"


echo
echo "Step 6: Calculate overlap with flagged regions of genome (eg repetitive), $(date)"
mkdir -p ${OUTPUT}/working_files/dels_filt1_flag_overlap
module load bedtools/2.22.0
parallel --halt 1 -j $NJOBS --verbose "intersectBed -wao -a ${OUTPUT}/working_files/dels_filt1/{}.dels.filt1.bed -b ${RESOURCES}/CNV_regions_to_flag.bed > ${OUTPUT}/working_files/dels_filt1_flag_overlap/{}.dels.filt1.flag.txt" :::: <(get_samples 6 short)
module unload bedtools
update_stage 6
echo "Step 6: finished, sample processing step in ${LIST} updated, $(date)"


echo
echo "Step 7: Do merging and pre-filtering frequency calculations, $(date)"
# First split dels by chromosome
mkdir -p ${OUTPUT}/working_files/dels_all_nohead_bychr
parallel --halt 1 -j $NJOBS --verbose "Rscript dels_all_nohead_bychr.R {}" :::: <(get_samples 7 short)

# Then aggregate samples
dels_all_header="MANTA_CHROM\tMANTA_START\tMANTA_END\tMANTA_DEL_SIZE\tCANVAS_CHROM\tCANVAS_START\tCANVAS_END\tCANVAS_DEL_SIZE\tCHROM\tSTART\tEND\tPAR\tvariant_caller\tILMN_ID\tBRIDGE_ID_LIB\tBRIDGE_ID\tPROJECT\tBATCH\tCHEM\tSEX_X\tMANTA_GT\tMANTA_GQ\tMANTA_FILTER\tMANTA_QUAL\tMANTA_PR\tMANTA_SR\tMANTA_REF\tMANTA_ALT\tMANTA_IMPRECISE\tMANTA_CIPOS\tMANTA_CIEND\tMANTA_CIGAR\tMANTA_MATEID\tMANTA_HOMLEN\tMANTA_HOMSEQ\tMANTA_SVINSLEN\tMANTA_SVINSSEQ\tMANTA_PAIR_COUNT\tMANTA_UPSTREAM_PAIR_COUNT\tMANTA_DOWNSTREAM_PAIR_COUNT\tMANTA_JUNC_QUAL\tCANVAS_GT\tCANVAS_RC\tCANVAS_BC\tCANVAS_CN\tCANVAS_FILTER\tCANVAS_QUAL\tCANVAS_CALL_DUP\tMANTA_PROP_COVERED\tCANVAS_PROP_COVERED\tMANTA_CANVAS_ILMN"
mkdir -p ${OUTPUT}/working_files/dels_filt1_by_chr
for c in $(echo $(seq 1 22) X); do
    echo -e $dels_all_header > ${OUTPUT}/working_files/dels_filt1_by_chr/dels.filt1.chr${c}.txt
done
parallel --halt 1 -j $NJOBS --verbose "cat ${OUTPUT}/working_files/dels_all_nohead_bychr/{1}.dels.chr{2}.txt >> ${OUTPUT}/working_files/dels_filt1_by_chr/dels.filt1.chr{2}.txt" ::: $(get_samples 1000 short) ::: $(echo $(seq 1 22) X)

# Then AF exact prefilt calculations (autosomes and X have different scripts)
parallel --halt 1 -j $NJOBS --verbose "Rscript dels_filt1_AFexact_aut_dels.R {}" ::: $(seq 1 22)
Rscript dels_filt1_AFexact_x_dels.R

# Make beds
parallel --halt 1 -j $NJOBS --verbose "Rscript dels_filt1_merged_make_beds.R {}" ::: $(echo $(seq 1 22) X)

# Do prefilt internal overlap (0.8) including read length and ethnicity specific
module load bedtools/2.22.0
for t in "samp" "AFR" "EAS" "EUR" "FIN" "OTH" "SAS" "100bp" "125bp" "150bp"; do
    echo "Internal overlap for $t"
    parallel --halt 1 -j $NJOBS --verbose "intersectBed -c -f 0.8 -r -a ${OUTPUT}/working_files/dels_filt1_by_chr/dels.filt1.chr{}.${t}.bed -b ${OUTPUT}/working_files/dels_filt1_by_chr/dels.filt1.chr{}.${t}.bed > ${OUTPUT}/working_files/dels_filt1_by_chr/dels.filt1.chr{}.${t}.overlap.internal.txt" ::: $(echo $(seq 1 22) X)
done
module unload bedtools
update_stage 7
echo "Step 7: finished, sample processing step in ${LIST} updated, $(date)"


echo
echo "Steps 1-7 of the pipeline have finished."
echo "Now we need to calculate the mean mapping quality and number of high-quality heterozygous SNPs within each deletion."
echo "Mapping quality: check, update if needed and run MAPQ_slurm_submit_sep.sh."
mkdir -p ${OUTPUT}/working_files/dels_filt1_ann
echo "Number of SNPs: check, update if needed and run hadoop_query.py at hdp-master2."


exit 0
