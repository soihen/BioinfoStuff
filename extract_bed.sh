#!/usr/bin/bash

# --------------------------- Description --------------------------- #
# Output genomic regions with coverage depth >= 20x from a BAM file   #
# ------------------------------------------------------------------- #
#                                                                     #
# Usage:                                                              #
#   [admin@kai]$bash extract_bed.sh [input BAM] [output BED filename] #
#                                                                     #
# ------------------------------------------------------------------- #

input_bam=$1
output_bed=$2

# -------------------- set parameters ------------------- #
bedtools="/public/software/bedtools2/bin/bedtools"
# ------------------------------------------------------- #


$bedtools genomecov -ibam ${input_bam} -bg | awk 'OFS="\t" {if ($4 > 20) print $1,$2,$3}' > tmp.bed
sort -k1,1 -k2,2n tmp.bed > tmp.sort.bed
bedtools merge -i tmp.sort.bed > ${output_bed}

rm tmp.bed tmp.sort.bed
