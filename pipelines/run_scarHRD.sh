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
#   [admin@kai]$bash run_sequenza.sh [input_folder] [output_folder]                     #
#                                                                                       #
# input_folder should contain pairs of matched tumor/normal pair-end sequenced Fastqs   #
#                                                                                       #
# Naming convention:                                                                    #
# tumor:   ${sampleID}_R[1|2].tumor.fastq.gz                                            #
# normal:  ${sampleID}_R[1|2].normal.fastq.gz                                           #
#                                                                                       #
# 1. Depends on DNA capture methods (i.e. targeted amplicon based or hybrid capture     #
#    based), user can choose to either perform or not perform deduplication step.       #
# ------------------------------------------------------------------------------------- #


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
gcwig="/public/database/GATK_Resource_Bundle/hg19/hg19.gc50.wig.gz"

# ------------------------------------------------------- #


input_folder=$1
output_folder=$2

if [[ ! -d $output_folder ]]; 
then
    mkdir $output_folder
fi

export SENTIEON_LICENSE=${sentieon_license};

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

    ($sentieon bwa mem -M -R "@RG\tID:${sampleID}\tSM:${sampleID}\tPL:illumina" \
    -t ${thread} -K 10000000 ${ref} \
    $trim_dir/${sampleID}_trim_R1.tumor.fastq.gz \
    $trim_dir/${sampleID}_trim_R2.tumor.fastq.gz \
    || echo -n 'error' ) \
    | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}.tumor.sorted.bam \
    -t ${thread} --sam2bam -i -;

    ($sentieon bwa mem -M -R "@RG\tID:${sampleID}\tSM:${sampleID}\tPL:illumina" \
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


    # step4 - generate seqz file
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- generating seqz file";

    cnv_dir=$output_folder/sequenza_cnv/;
    if [[ ! -d $cnv_dir ]]; then
        mkdir $cnv_dir
    fi

    if [[ $dedup == true ]]; then
        normal_bam=${align_dir}/${sampleID}.normal.sorted.dedup.bam
        tumor_bam=${align_dir}/${sampleID}.tumor.sorted.dedup.bam
    else
        normal_bam=${align_dir}/${sampleID}.normal.sorted.bam
        tumor_bam=${align_dir}/${sampleID}.tumor.sorted.bam
    fi

    sequenza-utils bam2seqz \
    --chromosome chr1 chr2 chr3 chr4 chr5 chr6 chr7 \
    chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 \
    chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
    --parallel ${thread} \
    -n $normal_bam \
    -t $tumor_bam \
    --fasta ${ref} -gc ${gcwig} -f illumina \
    -o $cnv_dir/${sampleID}.seqz.gz;

    ls $cnv_dir/${sampleID}*.seqz.gz | \
    xargs -i -P ${thread} sh -c \
    "sequenza-utils seqz_binning \
    -s {} -w 50 -o {}.bin.gz"

    zcat $cnv_dir/${sampleID}_chr{1..22}.seqz.gz.bin.gz \
    $cnv_dir/${sampleID}_chrX.seqz.gz.bin.gz \
    $cnv_dir/${sampleID}_chrY.seqz.gz.bin.gz \
    > $cnv_dir/${sampleID}.bin.seqz;

    sed -i -e 1b -e '/^chromosome/d' $cnv_dir/${sampleID}.bin.seqz;
    bgzip $cnv_dir/${sampleID}.bin.seqz;
    tabix -f -s 1 -b 2 -e 2 -S 1 $cnv_dir/${sampleID}.bin.seqz.gz;


    # step5 - process seqz data, normalization, segmentation and estimate cellularity and ploidy
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- CNV calling ";

    Rscript -e "library(sequenza); \
    extract <- sequenza.extract('$cnv_dir/${sampleID}.bin.seqz.gz'); \
    cp <- sequenza.fit(extract); \
    sequenza.results(sequenza.extract = extract, \
                     cp.table = cp, \
                     sample.id = '${sampleID}', \
                     out.dir = '$cnv_dir/${sampleID}')"

    # step6 - scarHRD
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- scarHRD ";
    
    Rscript -e "library(scarHRD); \
    hrd <- scar_score('$cnv_dir/${sampleID}.bin.seqz.gz', \
                reference = 'grch37', \
                seqz = 'TURE'); \
    write.csv(hrd, row.names = 'TURE', file = '$cnv_dir/${sampleID}.hrd_score.csv')"
done

 