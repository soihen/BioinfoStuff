#! /usr/bin/bash

# ----------------------------------- Description ------------------------------------- #
# Perform:                                                                              #
#   Trimming (fastp) +                                                                  #
#   Alignment (Sentieon-bwa) +                                                          #
#   statistic analysis +                                                                #
#   Deduplication (optional) +                                                          #
#   Allelic specific CNV calling (Sequencza)                                            #
# ------------------------------------------------------------------------------------- #
# Usage:                                                                                #
#   [admin@kai]$bash snv_calling_TNscope.sh [input_folder] [output_folder]              #
#                                                                                       #
# 1. You MUST specifiy whether use tumor-only mode ('single') or tumor/normal mode      #
#    ('match') of somatic calling of TNscope.                                           #
# 2. Depends on DNA capture methods (i.e. targeted amplicon based or hybrid capture     #
#    based), user can choose to either perform or not perform deduplication step.       #
# ------------------------------------------------------------------------------------- #


# -------------------- set parameters ------------------- #
# mode of somatic calling: (matched || single)
mode="matched"

# whether perform deduplicate step (true || false)
dedup=true

# path to software 
fastp="/data/ngs/softs/fastp/fastp"
sentieon="/data/ngs/softs/sentieon/sentieon-genomics-201808.08/bin/sentieon"
bcftools="/public/software/bcftools-1.9/bcftools"

sentieon_license="192.168.1.186:8990"
thread=16

# additional files
# reference fasta file
ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
dbsnp="/public/database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf"
k1="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf.gz"
k2="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
k3="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"

# ------------------------------------------------------- #


if [[  $1 == '-h'  ]]; then
    echo "Usage: ./snv_calling_TNscope.sh [input_folder] [output_folder]"
    echo "-------------------------------------------------------------------------"
    echo "[input_folder] should contain fastq files with following naming system:"
    echo "  single -- \${sampleID}_R[1|2].fastq.gz"
    echo "  matched -- \${sampleID}_R[1|2].tumor.fastq.gz"
    exit 0
fi



input_folder=$1
output_folder=$2

if [[ ! -d $output_folder ]]; 
then
    mkdir $output_folder
fi

export SENTIEON_LICENSE=${sentieon_license};


# tumor-normal matched mode
if [[  $mode == 'matched' ]]; then
    for ifile in $input_folder/*_R1.tumor.fastq.gz;
    do 
        sampleID=`basename ${ifile%%"_R1"*}`

        # step1 - trim reads
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- trimming reads";

        trim_dir=$output_folder/trim/;
        if [[ ! -d $trim_dir ]]; then
            mkdir $trim_dir
        fi

        $fastp --in1 $input_folder/${sampleID}_R1.tumor.fastq.gz \
        --in2 $input_folder/${sampleID}_R2.tumor.fastq.gz \
        --out1 $trim_dir/${sampleID}_trim_R1.tumor.fastq.gz \
        --out2 $trim_dir/${sampleID}_trim_R2.tumor.fastq.gz \
        -c --length_required 3 --detect_adapter_for_pe -p \
        --thread ${thread} \
        --html $trim_dir/${sampleID}.tumor.trim.html \
        --json $trim_dir/${sampleID}.tumor.trim.json;

        $fastp --in1 $input_folder/${sampleID}_R1.normal.fastq.gz \
        --in2 $input_folder/${sampleID}_R2.normal.fastq.gz \
        --out1 $trim_dir/${sampleID}_trim_R1.normal.fastq.gz \
        --out2 $trim_dir/${sampleID}_trim_R2.normal.fastq.gz \
        -c --length_required 3 --detect_adapter_for_pe -p \
        --thread ${thread} \
        --html $trim_dir/${sampleID}.normal.trim.html \
        --json $trim_dir/${sampleID}.normal.trim.json;


        # step2 - align & sort
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- alignment & sorting";

        align_dir=$output_folder/align/;
        if [[ ! -d $align_dir ]]; then
            mkdir $align_dir
        fi

        ($sentieon bwa mem -M -R "@RG\tID:${sampleID}.tumor\tSM:${sampleID}.tumor\tPL:illumina" \
        -t ${thread} -K 10000000 ${ref} \
        $trim_dir/${sampleID}_trim_R1.tumor.fastq.gz \
        $trim_dir/${sampleID}_trim_R2.tumor.fastq.gz \
        || echo -n 'error' ) \
        | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}.tumor.sorted.bam \
        -t ${thread} --sam2bam -i -;

        ($sentieon bwa mem -M -R "@RG\tID:${sampleID}.normal\tSM:${sampleID}.normal\tPL:illumina" \
        -t ${thread} -K 10000000 ${ref} \
        $trim_dir/${sampleID}_trim_R1.normal.fastq.gz \
        $trim_dir/${sampleID}_trim_R2.normal.fastq.gz \
        || echo -n 'error' ) \
        | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}.normal.sorted.bam \
        -t ${thread} --sam2bam -i -;


        # step3 (optional) - remove duplicates
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- deduplication";
        
        if [[ $dedup == true ]]; then
            ${sentieon} driver -t ${thread} \
            -i ${align_dir}/${sampleID}.tumor.sorted.bam \
            --algo LocusCollector \
            --fun score_info ${align_dir}/${sampleID}.tumor.score.txt;

            ${sentieon} driver -t ${thread} \
            -i ${align_dir}/${sampleID}.tumor.sorted.bam \
            --algo Dedup --rmdup \
            --score_info ${align_dir}/${sampleID}.tumor.score.txt \
            --metrics ${align_dir}/${sampleID}.tumor.dedup_metrics.txt \
            ${align_dir}/${sampleID}.tumor.sorted.dedup.bam;

            ${sentieon} driver -t ${thread} \
            -i ${align_dir}/${sampleID}.normal.sorted.bam \
            --algo LocusCollector \
            --fun score_info ${align_dir}/${sampleID}.normal.score.txt;

            ${sentieon} driver -t ${thread} \
            -i ${align_dir}/${sampleID}.normal.sorted.bam \
            --algo Dedup --rmdup \
            --score_info ${align_dir}/${sampleID}.normal.score.txt \
            --metrics ${align_dir}/${sampleID}.normal.dedup_metrics.txt \
            ${align_dir}/${sampleID}.normal.sorted.dedup.bam;
        fi


        # step4 - Base quality score recalibration 
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- BQSR";

        # set path of the bam file
        if [[ $dedup == true ]]; then
            normal_bam=${align_dir}/${sampleID}.normal.sorted.dedup.bam
            tumor_bam=${align_dir}/${sampleID}.tumor.sorted.dedup.bam
        else
            normal_bam=${align_dir}/${sampleID}.normal.sorted.bam
            tumor_bam=${align_dir}/${sampleID}.tumor.sorted.bam
        fi

        ${sentieon} driver -t ${thread} -r ${ref} \
        -i ${tumor_bam} --algo QualCal -k ${k1} -k ${k2} -k ${k3} \
        ${align_dir}/${sampleID}.tumor.recal.table;

        ${sentieon} driver -t ${thread} -r ${ref} \
        -i ${normal_bam} --algo QualCal -k ${k1} -k ${k2} -k ${k3} \
        ${align_dir}/${sampleID}.normal.recal.table;


        # step5 - variant calling TNscope
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- variant calling";

        snv_dir=$output_folder/snv/;
        if [[ ! -d $snv_dir ]]; then
            mkdir $snv_dir
        fi
        
        $sentieon driver -t ${thread} -r ${ref} \
        -i $tumor_bam -q ${align_dir}/${sampleID}.tumor.recal.table \
        -i $normal_bam -q ${align_dir}/${sampleID}.normal.recal.table \
        --algo TNscope \
        --tumor_sample ${sampleID}.tumor --normal_sample ${sampleID}.normal \
        --dbsnp ${dbsnp} ${snv_dir}/${sampleID}.raw.vcf
        
    done

# tumor-only mode 
elif [[  $mode == 'single'  ]]; then
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
        $input_folder/${sampleID}_R1.trimmed.fastq.gz \
        $input_folder/${sampleID}_R2.trimmed.fastq.gz \
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

        # set path of the bam file
        if [ $dedup == true ]; then
            bam=${align_dir}/${sampleID}.sorted.dedup.bam
        else
            bam=${align_dir}/${sampleID}.sorted.bam
        fi

        # step5 - BQSR
        ${sentieon} driver -t ${thread} -r ${ref} \
        -i ${bam} --algo QualCal -k ${k1} -k ${k2} -k ${k3} \
        ${align_dir}/${sampleID}.recal.table;

        # step6 - variant calling (TNscope)
        snv_dir=$output_folder/snv/;
        if [[ ! -d $snv_dir ]]; then
            mkdir $snv_dir
        fi
        ${sentieon} driver -t ${thread} -r ${ref} \
        -i ${bam} -q ${align_dir}/${sampleID}.recal.table \
        --algo TNscope --tumor_sample ${sampleID} \
        --dbsnp ${dbsnp} ${snv_dir}/${sampleID}.tumor.raw.vcf;

    done

fi
