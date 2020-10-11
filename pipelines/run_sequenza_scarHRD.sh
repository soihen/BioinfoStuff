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
# 2. As this assay designed to sequence on backbone region and gene region at different #
#    depths, two BED files must be prepared to indicate those two regions separately.   #
# ------------------------------------------------------------------------------------- #


# -------------------- set parameters ------------------- #
sentieon_license="192.168.1.186:8990"
thread=8

# whether perform deduplicate step (true || false)
dedup=true

# path to software 
fastp="/data/ngs/softs/fastp/fastp"
sentieon="/data/ngs/softs/sentieon/sentieon-genomics-201808.08/bin/sentieon"
bcftools="/public/software/bcftools-1.9/bcftools"
samtools="/public/software/samtools-1.9/samtools"
bedtools="/public/software/bedtools2/bin/bedtools"

# additional files
# reference fasta file
ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
gcwig="/public/database/GATK_Resource_Bundle/hg19/hg19.gc50.wig.gz"

backbone_bed="/public/home/kai/hrd_project/source_data/twist/backbone_100k_Liftover_success_hg19.bed"
gene_bed="/public/home/kai/hrd_project/source_data/twist/150genes.probe_merged.p.rev1_hg19.bed"

# ------------------------------------------------------- #



# -------------------------------
# create output directory
input_folder=$1
output_folder=$2


if [[ ! -d $output_folder ]]; 
then
    mkdir $output_folder
fi

trim_dir=$output_folder/trim/;
if [[ ! -d $trim_dir ]]; then
    mkdir $trim_dir
fi

align_dir=$output_folder/align/;
if [[ ! -d $align_dir ]]; then
    mkdir $align_dir
fi

cnv_dir=$output_folder/sequenza_cnv/;
if [[ ! -d $cnv_dir ]]; then
    mkdir $cnv_dir
fi

hrd_dir=${output_folder}/HRD/;
if [[  ! -d $hrd_dir  ]];
then 
    mkdir $hrd_dir
fi

qc_dir=$output_folder/qc/;
if [[ ! -d $qc_dir ]]; then
    mkdir $qc_dir
fi

# -------------------------------


# merge two BED files
cat ${backbone_bed} ${gene_bed} > $output_folder/concat.bed
sort -k1,1 -k2,2n $output_folder/concat.bed > $output_folder/concat_sorted.bed
$bedtools merge -i $output_folder/concat_sorted.bed > $output_folder/merged.bed
bed="$output_folder/merged.bed"

rm $output_folder/concat.bed $output_folder/concat_sorted.bed


echo "sampleID,fastq_size,raw_reads,raw_bases,clean_reads,clean_bases,\
qc30_rate,mapping_rate(%),on-target_percent(%),\
mean_depth,mean_dedup_depth,dup_rate(%),\
average_insert_size,std_insert_size,\
backbone_50x_depth_percent(%),backbone_100x_depth_percent(%),\
backbone_150x_depth_percent(%),gene_200x_depth_percent(%),\
gene_300x_depth_percent(%),gene_400x_depth_percent(%),\
gene_500x_depth_percent(%)" > $qc_dir/QC_summary.csv


export SENTIEON_LICENSE=${sentieon_license};

for ifile in $input_folder/*_R1.tumor.fastq.gz;
do 
    sampleID=`basename ${ifile%%"_R1"*}`

    # step1 - trim reads
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- trimming reads";

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


    # step4 - quality control
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- qc";

    $samtools depth -d 0 -b ${bed} ${align_dir}/${sampleID}.normal.sorted.bam > $qc_dir/${sampleID}.normal.depth.txt;
    $samtools depth -d 0 -b ${bed} ${align_dir}/${sampleID}.normal.sorted.dedup.bam > $qc_dir/${sampleID}.normal.dedup.depth.txt;
    $samtools depth -d 0 ${align_dir}/${sampleID}.normal.sorted.dedup.bam > $qc_dir/${sampleID}.normal.dedup.total_depth.txt;

    $samtools depth -d 0 -b ${bed} ${align_dir}/${sampleID}.tumor.sorted.bam > $qc_dir/${sampleID}.tumor.depth.txt;
    $samtools depth -d 0 -b ${bed} ${align_dir}/${sampleID}.tumor.sorted.dedup.bam > $qc_dir/${sampleID}.tumor.dedup.depth.txt;
    $samtools depth -d 0 ${align_dir}/${sampleID}.tumor.sorted.dedup.bam > $qc_dir/${sampleID}.tumor.dedup.total_depth.txt;

    $samtools stats -@ ${thread} ${align_dir}/${sampleID}.normal.sorted.dedup.bam > ${qc_dir}/${sampleID}.normal.stats.txt;
    $samtools stats -@ ${thread} ${align_dir}/${sampleID}.tumor.sorted.dedup.bam > ${qc_dir}/${sampleID}.tumor.stats.txt;

    tumor_r1=$(du $input_folder/${sampleID}_R1.tumor.fastq.gz -sh |awk '{print $1}');
    tumor_r2=$(du $input_folder/${sampleID}_R2.tumor.fastq.gz -sh |awk '{print $1}');
    normal_r1=$(du $input_folder/${sampleID}_R1.normal.fastq.gz -sh |awk '{print $1}');
    normal_r2=$(du $input_folder/${sampleID}_R2.normal.fastq.gz -sh |awk '{print $1}');

    normal_raw_reads=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
    print(fh['summary']['before_filtering']['total_reads'])"`

    normal_clean_reads=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
    print(fh['summary']['after_filtering']['total_reads'])"`

    normal_raw_bases=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
    print(fh['summary']['before_filtering']['total_bases'])"`

    normal_clean_bases=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
    print(fh['summary']['after_filtering']['total_bases'])"`

    tumor_raw_reads=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
    print(fh['summary']['before_filtering']['total_reads'])"`

    tumor_clean_reads=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
    print(fh['summary']['after_filtering']['total_reads'])"`

    tumor_raw_bases=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
    print(fh['summary']['before_filtering']['total_bases'])"`

    tumor_clean_bases=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
    print(fh['summary']['after_filtering']['total_bases'])"`

    normal_qc_rate=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
    print(fh['summary']['after_filtering']['q30_rate'])"`

    tumor_qc_rate=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
    print(fh['summary']['after_filtering']['q30_rate'])"`

    normal_mapped_reads=$(awk -F"\t" '$2 == "reads mapped:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt);
    normal_all_reads=$(awk -F"\t" '$2 == "sequences:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt);
    normal_mapping_rate=`bc <<< "scale=4; ${normal_mapped_reads} / ${normal_all_reads} * 100"`

    tumor_mapped_reads=$(awk -F"\t" '$2 == "reads mapped:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt);
    tumor_all_reads=$(awk -F"\t" '$2 == "sequences:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt);
    tumor_mapping_rate=`bc <<< "scale=4; ${tumor_mapped_reads} / ${tumor_all_reads} * 100"`

    normal_mean_depth=$(awk '{sum+=$3}END{print sum/NR}' $qc_dir/${sampleID}.normal.depth.txt)
    normal_mean_dedup_depth=$(awk '{sum+=$3}END{print sum/NR}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
    normal_dup_rate=`bc <<< "scale=4; (1 - ${normal_mean_dedup_depth} / ${normal_mean_depth}) * 100"`
    
    tumor_mean_depth=$(awk '{sum+=$3}END{print sum/NR}' $qc_dir/${sampleID}.tumor.depth.txt)
    tumor_mean_dedup_depth=$(awk '{sum+=$3}END{print sum/NR}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
    tumor_dup_rate=`bc <<< "scale=4; (1 - ${tumor_mean_dedup_depth} / ${tumor_mean_depth}) * 100"`

    normal_total_depth=$(awk '{sum+=$3}END{print sum}' $qc_dir/${sampleID}.normal.dedup.total_depth.txt)
    normal_on_target=$(awk -v all_bases=${normal_total_depth} '{sum+=$3}END{print sum/all_bases*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)

    tumor_total_depth=$(awk '{sum+=$3}END{print sum}' $qc_dir/${sampleID}.tumor.dedup.total_depth.txt)
    tumor_on_target=$(awk -v all_bases=${tumor_total_depth} '{sum+=$3}END{print sum/all_bases*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)

    normal_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt);
    normal_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt);

    tumor_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt);
    tumor_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt);


    normal_50x=$(awk 'BEGIN {count=0} {if ($3 > 50) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
    normal_100x=$(awk 'BEGIN {count=0} {if ($3 > 100) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
    normal_150x=$(awk 'BEGIN {count=0} {if ($3 > 150) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
    normal_200x=$(awk 'BEGIN {count=0} {if ($3 > 200) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
    normal_300x=$(awk 'BEGIN {count=0} {if ($3 > 300) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
    normal_400x=$(awk 'BEGIN {count=0} {if ($3 > 400) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
    normal_500x=$(awk 'BEGIN {count=0} {if ($3 > 500) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)

    tumor_50x=$(awk 'BEGIN {count=0} {if ($3 > 50) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
    tumor_100x=$(awk 'BEGIN {count=0} {if ($3 > 100) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
    tumor_150x=$(awk 'BEGIN {count=0} {if ($3 > 150) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
    tumor_200x=$(awk 'BEGIN {count=0} {if ($3 > 200) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
    tumor_300x=$(awk 'BEGIN {count=0} {if ($3 > 300) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
    tumor_400x=$(awk 'BEGIN {count=0} {if ($3 > 400) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
    tumor_500x=$(awk 'BEGIN {count=0} {if ($3 > 500) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)

    echo "${sampleID}.normal,${normal_r1}/${normal_r2},${normal_raw_reads},${normal_raw_bases},${normal_clean_reads},${normal_clean_bases},${normal_qc_rate},${normal_mapping_rate},${normal_on_target},${normal_mean_depth},${normal_mean_dedup_depth},${normal_dup_rate::-2},${normal_insert_size}, ${normal_insert_std},${normal_50x},${normal_100x},${normal_150x},${normal_200x},${normal_300x},${normal_400x},${normal_500x}" \
    >> ${qc_dir}/QC_summary.csv

    echo "${sampleID}.tumor,${tumor_r1}/${tumor_r2},${tumor_raw_reads},${tumor_raw_bases},${tumor_clean_reads},${tumor_clean_bases},${tumor_qc_rate},${tumor_mapping_rate},${tumor_on_target},${tumor_mean_depth},${tumor_mean_dedup_depth},${tumor_dup_rate::-2},${tumor_insert_size}, ${tumor_insert_std},${tumor_50x},${tumor_100x},${tumor_150x},${tumor_200x},${tumor_300x},${tumor_400x},${tumor_500x}" \
    >> ${qc_dir}/QC_summary.csv


    # step5 - generate seqz file
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- generating seqz file";

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


    # step6 - process seqz data, normalization, segmentation and estimate cellularity and ploidy
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- CNV calling ";

    Rscript -e "library(sequenza); \
    extract <- sequenza.extract('$cnv_dir/${sampleID}.bin.seqz.gz'); \
    cp <- sequenza.fit(extract); \
    sequenza.results(sequenza.extract = extract, \
                     cp.table = cp, \
                     sample.id = '${sampleID}', \
                     out.dir = '$cnv_dir/${sampleID}')"

    # step7 - scarHRD
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- scarHRD ";

    ploidy=$(awk 'NR ==2 {print $2}' $cnv_dir/${sampleID}/${sampleID}_alternative_solutions.txt);

    awk -F"\t" -v ploidy=${ploidy} -v sample=${sampleID} \
    'BEGIN {OFS="\t"; print "SampleID","Chromosome","Start_position","End_position","total_cn","A_cn","B_cn","ploidy"};
     NR!=1 {OFS="\t"; print sample,$1,$2,$3,$10,$11,$12,ploidy}' \
    $cnv_dir/${sampleID}/${sampleID}_segments.txt > $hrd_dir/${sampleID}.pre_hrd.tsv;

    Rscript -e "library(scarHRD); \
    hrd <- scar_score('$hrd_dir/${sampleID}.pre_hrd.tsv', \
                reference = 'grch37', \
                seqz = 'FALSE'); \
    write.csv(hrd, row.names = '${sampleID}', file = '$hrd_dir/${sampleID}.hrd.csv')"
done

echo "LOGGING: `date --rfc-3339=seconds` -- Analysis finished";