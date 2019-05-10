#!/bin/bash

########################################################################
# Written by Olga Shamardina
#
# Submits Get_MAPQ_perDEL.sh as slurm job for all samples in the manifest
#
########################################################################

### Estimation of run time based on linear fit for 6 points + 50%: #####
# int((a + b*NoD)/60*1.5 + 1)
LFA=233.1
LFB=0.09229
########################################################################

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
cd ${DIR}
source config.sh

sed 1d ${MANIFEST} | cut -f 2,10 | while read -a B; do
    BED=${OUTPUT}/working_files/dels_filt1/${B[0]}.dels.filt1.bed
    BAM=${B[1]}
    OUT=${OUTPUT}/working_files/dels_filt1_ann/${B[0]}.dels.filt1.readinfo.txt
    T=$(echo "print int(($LFA + $LFB*$(cat $BED | wc -l))/60*1.5 + 1)" | python)
    sbatch -J ${B[0]} --nodes=1 --ntasks=1 --time=${T}:00 --mail-type=FAIL -o "slurm/slurm-%j.out" Get_MAPQ_perDEL.sh $BED $BAM $OUT
done
