#! /usr/bin/bash

# ----------------------------------- Description ------------------------------------- #
#      Noted: You need to purchase sentieon's license for use of this pipeline.         #
# ------------------------------------------------------------------------------------- #
# Perform:                                                                              #
#   Trimming (fastp) +                                                                  #
#   Alignment (bwa) +                                                                   #
#   statistic analysis +                                                                #
#   Deduplication (optional) +                                                          #
#   Indel-rearrangement +                                                               #
#   BQSR +                                                                              #
#   SNP/INDEL calling (Haplotyper)  +                                                   #
#   Variant Position Normalisation (bcftools) +                                         #
#   Filter based on depth coverage (bcftools)                                           #
# ------------------------------------------------------------------------------------- #
# Usage:                                                                                #
#   [admin@kai]$bash snp_sentieon.sh [input_folder] [output_folder]                     #
#                                                                                       #
# input_folder should contain pair-end sequenced Fastq files.                           #
# Each sample is expected to have two Fastq files, with naming convention of:           #
# ${sampleID}_R[1|2].fastq.gz                                                           #
# ------------------------------------------------------------------------------------- #
# 1. Depends on DNA capture methods (i.e. targeted amplicon based or hybrid capture     #
#    based), user can choose to either perform or not perform deduplication step.       #
# 2. BED file is also NOT used by default,                                              #
#    as BED file will only specify regions of variant calling.                          #
#    This aim can be achieve equally by filtering obtained VCF file based on            #
#    genomic coordinate.                                                                #
# ------------------------------------------------------------------------------------- #

input_folder=$1
output_folder=$2

# -------------------- set parameters ------------------- #
sentieon_license="192.168.1.186:8990"
thread=16
# whether perform deduplicate step (true || false)
dedup=true
# path to software 
fastp="/data/ngs/softs/fastp/fastp"
sentieon="/data/ngs/softs/sentieon/sentieon-genomics-201808.08/bin/sentieon"
bcftools="/public/software/bcftools-1.9/bcftools"
# additional files
# reference fasta file
ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
# sentieon provided three additional vcf file that can be used in Base Recalibration
k1="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf.gz"
k2="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
k3="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
# dbSNP database
dbsnp="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf.gz"
# ------------------------------------------------------- #


if [[ ! -d $output_folder ]]; 
then
    mkdir $output_folder
fi


# ------------------------------------------ #

export SENTIEON_LICENSE=${sentieon_license};

for ifile in $input_folder/*_R1.fastq.gz;
do 
    sampleID=`basename ${ifile%%"_R1"*}`

    # step1 - trim reads
    trim_dir=$output_folder/trim/;
    if [[ ! -d $trim_dir ]]; then
        mkdir $trim_dir
    fi

    $fastp --in1 $input_folder/${sampleID}_R1.fastq.gz \
    --in2 $input_folder/${sampleID}_R2.fastq.gz \
    --out1 $trim_dir/${sampleID}_trim_R1.fastq.gz \
    --out2 $trim_dir/${sampleID}_trim_R2.fastq.gz \
    -c --length_required 3 --detect_adapter_for_pe -p \
    --thread ${thread} \
    --html $trim_dir/${sampleID}.trim.html \
    --json $trim_dir/${sampleID}.trim.json;

    # step2 - align & sort
    align_dir=$output_folder/align/;
    if [[ ! -d $align_dir ]]; then
        mkdir $align_dir
    fi

    ($sentieon bwa mem -M -R "@RG\tID:${sampleID}\tSM:${sampleID}\tPL:illumina" \
    -t ${thread} -K 10000000 ${ref} \
    $trim_dir/${sampleID}_trim_R1.fastq.gz \
    $trim_dir/${sampleID}_trim_R2.fastq.gz \
    || echo -n 'error' ) \
    | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}.sorted.bam \
    -t ${thread} --sam2bam -i -;

    # step3 - statistic analysis
    ${sentieon} driver -t ${thread} -r ${ref} -i ${align_dir}/${sampleID}.sorted.bam \
    --algo GCBias --summary ${align_dir}/${sampleID}.gc_summary.txt \
    ${align_dir}/${sampleID}.gc_metric.txt \
    --algo MeanQualityByCycle ${align_dir}/${sampleID}.mq_metric.txt \
    --algo QualDistribution ${align_dir}/${sampleID}.qd_metric.txt \
    --algo InsertSizeMetricAlgo ${align_dir}/${sampleID}.is_metric.txt \
    --algo AlignmentStat ${align_dir}/${sampleID}.aln_metric.txt;

    # step4 (optional) - remove duplicates
    if [ $dedup == true ]; then
        ${sentieon} driver -t ${thread} \
        -i ${align_dir}/${sampleID}.sorted.bam \
        --algo LocusCollector \
        --fun score_info ${align_dir}/${sampleID}.score.txt;

        ${sentieon} driver -t ${thread} \
        -i ${align_dir}/${sampleID}.sorted.bam \
        --algo Dedup --score_info ${align_dir}/${sampleID}.score.txt \
        --metrics ${align_dir}/${sampleID}.dedup_metrics.txt \
        ${align_dir}/${sampleID}.sorted.dedup.bam;
    fi

    # step5 - indel re-alignment
    if [ $dedup == true ]; then
        ${sentieon} driver -t ${thread} -r ${ref} \
        -i ${align_dir}/${sampleID}.sorted.dedup.bam \
        --algo Realigner -k ${k1} -k ${k2} -k ${k3} \
        ${align_dir}/${sampleID}.realigned.bam
    else
        ${sentieon} driver -t ${thread} -r ${ref} \
        -i ${align_dir}/${sampleID}.sorted.bam \
        --algo Realigner -k ${k1} -k ${k2} -k ${k3} \
        ${align_dir}/${sampleID}.realigned.bam
    fi

    # step6 - BQSR
    ${sentieon} driver -t ${thread} -r ${ref} \
    -i ${align_dir}/${sampleID}.realigned.bam \
    --algo QualCal -k ${k1} -k ${k2} -k ${k3} ${align_dir}/${sampleID}.recal.table;

    ${sentieon} driver -t ${thread} -r ${ref} -i ${align_dir}/${sampleID}.realigned.bam \
    -q ${align_dir}/${sampleID}.recal.table --algo QualCal -k ${k1} -k ${k2} -k ${k3} \
    ${align_dir}/${sampleID}.recal.post.table \
    --algo ReadWriter ${align_dir}/${sampleID}.recal.bam;

    # step7 - Haplotyper
    snp_dir=$output_folder/snp/;
    if [[ ! -d $snp_dir ]]; then
        mkdir $snp_dir
    fi

    ${sentieon} driver -r ${ref} -t ${thread} \
    -i ${align_dir}/${sampleID}.recal.bam \
    --algo Haplotyper --emit_conf=10 --call_conf=10 \
    -d ${dbsnp} \
    ${snp_dir}/${sampleID}.germline.raw.vcf;

    # step8 - normalise variants positions
    $bcftools norm -m -both -f ${ref} \
    ${snp_dir}/${sampleID}.germline.raw.vcf \
    -o ${snp_dir}/${sampleID}.germline.norm.vcf

    # step9 - filter variants without enough coverage
    $bcftools filter -i "FORMAT/DP>20" \
    ${snp_dir}/${sampleID}.germline.norm.vcf > \
    ${snp_dir}/${sampleID}.germline.filtered.vcf
done

