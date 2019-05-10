#!/bin/bash

########################################################################
# Written by Olga Shamardina
#
# In a BAM, finds number of reads and their mean mapping quality for regions in BED
#
########################################################################

set -e
set -o pipefail
set -u

BED=$1
BAM=$2
OUT=$3

module load samtools

echo "Start time: $(date --rfc-3339='seconds')"
echo "Sample: $SLURM_JOB_NAME"
echo "Input: $BED $BAM $OUT"
sort -V $BED | while read -a b; do
    samtools view $BAM ${b[0]}:${b[1]}-${b[2]} | awk -F "\t" '{a=a+$5} END {if (NR==0) M="NA"; else M=a/NR; print "'${b[0]}'","'${b[1]}'","'${b[2]}'",NR,M}'
done > $OUT
echo "Finish time: $(date --rfc-3339='seconds')"

exit 0
