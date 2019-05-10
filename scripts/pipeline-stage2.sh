#!/bin/bash

### Stage 2
# To be run after Stage 1, MAPQ_slurm_submit_sep.sh and hadoop_query.py have finished successfully.
########

set -e
set -o pipefail
set -u


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
cd ${DIR}
source config.sh
source functions.sh
module load parallel


echo "Step 8 (\"Hadoop\"): after the results are ready, split them by sample"
mkdir -p ${OUTPUT}/working_files/dels_filt1_snps
Rscript split_hadoop_results.R
update_stage 8
echo "Step 8: finished, sample processing step in ${LIST} updated, $(date)"


echo "Step 9: Annotate deletions, $(date)"
# Prepare file with genes for annotation:
# Merging:
for c in $(seq 1 22) X; do sort -k2,2g -k3,3g ${OUTPUT}/working_files/dels_filt1_by_chr/dels.filt1.chr${c}.bed | uniq; done > ${OUTPUT}/working_files/dels_filt1_by_chr/dels.filt1.chrAll.bed
# Making BED files for gene annotation:
awk -F \\t -v OFS=\\t 'NR>1 && $3!="Y" && $3!="MT" {print $3,$4,$5,$1}' ${RESOURCES}/hsapiens_gene_ensembl.txt | sort -k1,1V -k2,2g -k3,3g | uniq > ${OUTPUT}/working_files/dels_filt1_by_chr/genes.bed
awk -F \\t -v OFS=\\t '$1!="Y" && $1!="MT" {print $1,$2,$3,$5}' ${RESOURCES}/GRCh37_exons_prot_coding_HGNC_canonical_CDS.bed | sort -k1,1V -k2,2g -k3,3g | uniq > ${OUTPUT}/working_files/dels_filt1_by_chr/genes_strict.bed
module load bedtools
for annot in genes genes_strict; do
    # Annotating:
    if [[ $annot == "genes" ]]; then
	shift1=0
    else
	shift1=1
    fi
    bedtools intersect -sorted -loj -a <(awk -F \\t -v OFS=\\t -v shift1=$shift1 '{$2+=shift1; $3+=1; print}' ${OUTPUT}/working_files/dels_filt1_by_chr/dels.filt1.chrAll.bed) -b <(awk -F \\t -v OFS=\\t '{$3+=1; print}' ${OUTPUT}/working_files/dels_filt1_by_chr/${annot}.bed) | cut -f 1-3,7 | sort | uniq | sort -k1,1V -k2,2g -k3,3g | awk -F \\t -v OFS=\\t -v shift1=$shift1 '{$2-=shift1; print}' > ${OUTPUT}/working_files/dels_filt1_by_chr/dels.filt1.chrAll_intersect_${annot}.bed
    # Merging annotations:
    python merge_del_gene_info.py ${OUTPUT}/working_files/dels_filt1_by_chr/dels.filt1.chrAll_intersect_${annot}.bed
    awk -F \\t -v OFS=\\t '{if (NR == 1) print $0,"VAR"; else {$3-=1; print $0,$1" "$2" "$3}}' ${OUTPUT}/working_files/dels_filt1_by_chr/dels.filt1.chrAll_intersect_${annot}_merged.bed > ${OUTPUT}/working_files/dels_filt1_by_chr/dels.filt1.chrAll_intersect_${annot}_merged_VAR.bed
done
module unload bedtools

# Uses output of previous steps to annotate main deletions file
parallel --halt 1 -j $NJOBS --verbose "Rscript cnvs_ann_for_filt2.R {}" :::: <(get_samples 9 short)
update_stage 9
echo "Step 9: finished, sample processing step in ${LIST} updated, $(date)"


echo "Step 10: Make list of flagged samples $(date)"
mkdir -p ${OUTPUT}/working_files/flagging
any_sample=$(get_samples 1000 short | sed -n 1p)
CANVAS_QUAL=$(get_column_number ${OUTPUT}/working_files/dels_filt1_ann/${any_sample}.dels.filt1.ann1.txt CANVAS_QUAL)
MANTA_FILTER=$(get_column_number ${OUTPUT}/working_files/dels_filt1_ann/${any_sample}.dels.filt1.ann1.txt MANTA_FILTER)
module load bedtools
parallel --halt 1 -j $NJOBS --verbose "awk -F \\\\t -v OFS=\\\\t -v CQ=${CANVAS_QUAL} -v MF=${MANTA_FILTER} 'NR>1 && ((\$4!=\"Canvas\" || \$CQ>=10) && (\$4!=\"Manta\" || \$MF==\"PASS\")) {print \$1,\$2-1,\$3,\$4}' ${OUTPUT}/working_files/dels_filt1_ann/{}.dels.filt1.ann1.txt | sort -V | bedtools merge -c 4 -o collapse -i - | awk -F \\\\t -v OFS=\\\\t '{len=len+\$3-\$2; split(\$4,a,\",\"); for (i in a) {if (a[i]==\"Canvas\") count++}} END {print \"'{}'\", len, count}'" :::: <(get_samples 1000 short) > ${OUTPUT}/working_files/flagging/deletions_counts.txt
module unload bedtools
Rscript flag_samples.R
update_stage 10
echo "Step 10: finished, sample processing step in ${LIST} updated, $(date)"


echo "Step 11: Annotate deletions $(date)"
# Calculate the 20%ile of MANTA_GQ and MANTA_QUAL on dels_filt1_by_chr dataset:
Rscript manta_quantiles.R
# Make strict dataset and lenient datasets:
mkdir -p ${OUTPUT}/working_files/dels_filt2_str
mkdir -p ${OUTPUT}/working_files/dels_filt2_lnt
parallel --halt 1 -j $NJOBS --verbose "Rscript cnvs_filt2_dels.R {}" :::: <(get_samples 11.1 short)
update_stage 11.1

# Merge results and split by chr
header="CHROM\tSTART\tEND\tvariant_caller\tsize\tMANTA_CHROM\tMANTA_START\tMANTA_END\tMANTA_DEL_SIZE\tCANVAS_CHROM\tCANVAS_START\tCANVAS_END\tCANVAS_DEL_SIZE\tILMN_ID\tBRIDGE_ID\tPROJECT\tCHEM\tSEX_X\tETHNICITY\tUNRELATED\tAFFECTED\tGENES\tGENES_STRICT\tMANTA_GT\tMANTA_GQ\tMANTA_FILTER\tMANTA_QUAL\tMANTA_PR\tMANTA_SR\tMANTA_REF\tMANTA_ALT\tMANTA_IMPRECISE\tMANTA_CIPOS\tMANTA_CIEND\tMANTA_CIGAR\tMANTA_MATEID\tMANTA_HOMLEN\tMANTA_HOMSEQ\tMANTA_SVINSLEN\tMANTA_SVINSSEQ\tMANTA_PAIR_COUNT\tMANTA_UPSTREAM_PAIR_COUNT\tMANTA_DOWNSTREAM_PAIR_COUNT\tMANTA_JUNC_QUAL\tCANVAS_GT\tCANVAS_RC\tCANVAS_BC\tCANVAS_CN\tCANVAS_FILTER\tCANVAS_QUAL\tCANVAS_CALL_DUP\tMANTA_PROP_COVERED\tCANVAS_PROP_COVERED\tFLAGREG_N\tFLAGREG_MAX_TYPE\tFLAGREG_MAX_PROP_CNV\tDUKE0_TOTAL_PROP_CNV\tDUKE0_START\tDUKE0_STOP\tREADS_N\tREAD_DENSITY\tREADS_MEAN_MAPQ\tSNPS_N_HET_PASS\tSNPS_DENSITY_HET_PASS\tFLAG_SAMPLE_QUAL"

# Lenient
mkdir -p ${OUTPUT}/working_files/dels_filt2_lnt_by_chr
parallel --halt 1 -j $NJOBS --verbose "echo -e $header > ${OUTPUT}/working_files/dels_filt2_lnt_by_chr/dels.filt2.lnt.chr{}.txt" ::: $(echo $(seq 1 22) X)
parallel --halt 1 -j $NJOBS --verbose "cat ${OUTPUT}/working_files/dels_filt2_lnt/{1}.dels.filt2.lnt.chr{2}.txt >> ${OUTPUT}/working_files/dels_filt2_lnt_by_chr/dels.filt2.lnt.chr{2}.txt" ::: $(get_samples 11.2 short) ::: $(echo $(seq 1 22) X)
update_stage 11.2

# Strict
mkdir -p ${OUTPUT}/working_files/dels_filt2_str_by_chr
parallel --halt 1 -j $NJOBS --verbose "echo -e $header > ${OUTPUT}/working_files/dels_filt2_str_by_chr/dels.filt2.str.chr{}.txt" ::: $(echo $(seq 1 22) X)
parallel --halt 1 -j $NJOBS --verbose "cat ${OUTPUT}/working_files/dels_filt2_str/{1}.dels.filt2.str.chr{2}.txt >> ${OUTPUT}/working_files/dels_filt2_str_by_chr/dels.filt2.str.chr{2}.txt" ::: $(get_samples 11.3 short) ::: $(echo $(seq 1 22) X)
update_stage 11.3
echo "Step 11: finished, sample processing step in ${LIST} updated, $(date)"


echo "Step 12: $(date)"
mkdir -p ${OUTPUT}/dataset/strict
mkdir -p ${OUTPUT}/dataset/strict_nonGEL
parallel --halt 1 -j $NJOBS --verbose "Rscript dels_filt2_final_annotation.R {} strict" ::: $(echo $(seq 1 22) X)

mkdir -p ${OUTPUT}/dataset/lenient
mkdir -p ${OUTPUT}/dataset/lenient_nonGEL
parallel --halt 1 -j $NJOBS --verbose "Rscript dels_filt2_final_annotation.R {} lenient" ::: $(echo $(seq 1 22) X)
update_stage 12
echo "Step 12: finished, sample processing step in ${LIST} updated, $(date)"


exit 0
