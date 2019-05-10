#!/bin/bash

########################################################################
# Written by Olga Shamardina
#
# Auxiliary functions for the pipeline
#
########################################################################

set -e
set -o pipefail
set -u


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
cd ${DIR}
source config.sh


function get_samples {
    stage=$1
    column=$2
    if [[ $column == "vcf" ]]; then
	F=$((1 + $(sed 1q ${MANIFEST} | tr \\t \\n | grep -nx SVVCF_PATH | cut -f 1 -d :)))
    elif [[ $column == "bam" ]]; then
	F=$((1 + $(sed 1q ${MANIFEST} | tr \\t \\n | grep -nx BAM_PATH | cut -f 1 -d :)))
    elif [[ $column == "id" ]]; then
	F=1
    elif [[ $column == "short" ]]; then
	F=3
    elif [[ $column == "iid" ]]; then
	F=4
    else
	F='1-'
    fi
    join -t $'\t' <(awk -F "\t" '$2<'${stage}'' ${LIST} | sort) <(sed 1d ${MANIFEST} | sort) | cut -f $F
}
export -f get_samples


function update_stage {
    stage=$1
    awk -F "\t" -v OFS="\t" '{if ($2<'${stage}') $2='${stage}'; print $0}' ${LIST} > ${LIST}-new
    mv ${LIST}-new ${LIST}
}
export -f update_stage


function vcf2bed {
    awk -F \\t -v OFS=\\t '{print $1,$2-1,$3}' $1
}
export -f vcf2bed

function get_column_number {
    input=$1
    column=$2
    sed 1q $input | tr \\t \\n | grep -w -n -m 1 $column | cut -f 1 -d :
}
export -f get_column_number
