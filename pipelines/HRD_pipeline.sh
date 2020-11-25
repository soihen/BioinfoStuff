#! /usr/bin/bash

# ----------------------------------- Description ------------------------------------- #
# Perform:                                                                              #
#   Trimming (fastp) +                                                                  #
#   Alignment (Sentieon-bwa) +                                                          #
#   statistic analysis +                                                                #
#   Deduplication (optional) +                                                          #
#   Allelic specific CNV calling (Sequencza) +                                          #
#   Calculate HRD score (scarHRD) +                                                     #
#   SNV calling on HR-related genes                                                     #
# ------------------------------------------------------------------------------------- #
# Usage:                                                                                #
#   [admin@kai]$bash HRD_pipeline.sh [input_folder] [output_folder] \                   #
#               [BED_loci] [BED_gene]                                                   # 
#                                                                                       #
# input_folder should contain pairs of matched tumor/normal pair-end sequenced Fastqs   #
#                                                                                       #
# Naming convention:                                                                    #
# tumor:   ${sampleID}_tumor_R[1|2].fastq.gz                                            #
# normal:  ${sampleID}_normal_R[1|2].fastq.gz                                           #
#                                                                                       #
# 1. Depends on DNA capture methods (i.e. targeted amplicon based or hybrid capture     #
#    based), user can choose to either perform or not perform deduplication step.       #
# 2. As this assay designed to sequence on backbone region and gene region at different #
#    depths, two BED files must be prepared to indicate those two regions separately.   #
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
samtools="/public/software/samtools-1.9/samtools"
bamdst="/public/software/bamdst/bamdst"
bgzip="/public/software/htslib-1.9/bin/bgzip"
tabix="/public/software/htslib-1.9/bin/tabix"
vep="/public/software/98vep/ensembl-vep/vep"

merge_mnv="/public/home/kai/BioinfoStuff/tertiary_analysis/merge_mnv.py"
anno_hgvs="/public/home/kai/BioinfoStuff/tertiary_analysis/anno_hgvs.py"

# path to databases
ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
dbsnp="/public/database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf"
k1="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf.gz"
k2="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
k3="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
refflat="/public/database/GATK_Resource_Bundle/refFlat.txt"
vep_dir="/public/software/vep_98/"
cache_version="98"
clinic_transcripts="/public/home/kai/database/LRG/parsed_LRG.tsv"
gcwig="/public/database/GATK_Resource_Bundle/hg19/hg19.gc50.wig.gz"

# ------------------------------------------------------- #

if [[  $1 == '-h'  ]]; then
    echo "Usage: ./HRD_pipeline.sh [input_folder] [output_folder] [BED1] [BED2] [BED3]"
    echo "-------------------------------------------------------------------------"
    echo "[input_folder] should contain fastq files with following naming system:"
    echo "               single -- \${sampleID}_R[1|2].fastq.gz"
    echo "               matched -- \${sampleID}_normal/tumor_R[1|2].fastq.gz"
    echo "[BED1] -- SNP loci"
    echo "[BED2] -- Genes"
    echo "[BED3] -- merged BED file"
    exit 0
fi


# -------------------------------
# create output directory
input_folder=$1
output_folder=$2
bed1=$3
bed2=$4
bed3=$5


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

snv_dir=${output_folder}/SNV/;
if [[  ! -d $snv_dir  ]];
then 
    mkdir $snv_dir
fi

# -------------------------------
echo "LOGGING: `date --rfc-3339=seconds` -- Analysis started"
echo "========================================================"
echo "LOGGING: -- settings -- input folder -- ${input_folder}"
echo "LOGGING: -- settings -- output folder -- ${output_folder}"
echo "LOGGING: -- settings -- BED file for SNP loci -- ${bed1}"
echo "LOGGING: -- settings -- BED file for HR genes -- ${bed2}"
echo "LOGGING: -- settings -- merged BED file -- ${bed3}"
echo "========================================================"


echo "sampleID,fastq_size,raw_reads,raw_bases,clean_reads,clean_bases,\
qc30_rate,mapping_rate(%),on-target_percent(%),dup_rate(%),\
mean_depth(SNP loci),mean_dedup_depth(SNP loci),\
mean_depth(HR genes),mean_dedup_depth(HR genes),\
average_insert_size,std_insert_size,\
Uniformity_0.1X(SNP loci)(%),Uniformity_0.2X(SNP loci)(%),\
Uniformity_0.5X(SNP loci)(%),Uniformity_1X(SNP loci)(%),\
Uniformity_0.1X(HR genes)(%),Uniformity_0.2X(HR genes)(%),\
Uniformity_0.5X(HR genes)(%),Uniformity_1X(HR genes)(%),\
100x_depth_percent(SNP loci)(%),150x_depth_percent(SNP loci)(%),200x_depth_percent(SNP loci)(%),\
300x_depth_percent(HR genes)(%),400x_depth_percent(HR genes)(%),500x_depth_percent(HR genes)(%)" \
> $qc_dir/QC_summary.csv


export SENTIEON_LICENSE=${sentieon_license};

for ifile in $input_folder/*_tumor_R1.fastq.gz;
do 
    sampleID=`basename ${ifile%%"_tumor"*}`

    # step1 - trim reads
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- trimming reads";

    $fastp --in1 $input_folder/${sampleID}_tumor_R1.fastq.gz \
    --in2 $input_folder/${sampleID}_tumor_R2.fastq.gz \
    --out1 $trim_dir/${sampleID}_trim_R1.tumor.fastq.gz \
    --out2 $trim_dir/${sampleID}_trim_R2.tumor.fastq.gz \
    -c --length_required 3 --detect_adapter_for_pe -p \
    --thread ${thread} \
    --html $trim_dir/${sampleID}.tumor.trim.html \
    --json $trim_dir/${sampleID}.tumor.trim.json;

    $fastp --in1 $input_folder/${sampleID}_normal_R1.fastq.gz \
    --in2 $input_folder/${sampleID}_normal_R2.fastq.gz \
    --out1 $trim_dir/${sampleID}_trim_R1.normal.fastq.gz \
    --out2 $trim_dir/${sampleID}_trim_R2.normal.fastq.gz \
    -c --length_required 3 --detect_adapter_for_pe -p \
    --thread ${thread} \
    --html $trim_dir/${sampleID}.normal.trim.html \
    --json $trim_dir/${sampleID}.normal.trim.json;


    # step2 - align & sort
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- alignment & sorting";

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
        --algo Dedup \
        --score_info ${align_dir}/${sampleID}.tumor.score.txt \
        --metrics ${align_dir}/${sampleID}.tumor.dedup_metrics.txt \
        ${align_dir}/${sampleID}.tumor.sorted.dedup.bam;

        ${sentieon} driver -t ${thread} \
        -i ${align_dir}/${sampleID}.normal.sorted.bam \
        --algo LocusCollector \
        --fun score_info ${align_dir}/${sampleID}.normal.score.txt;

        ${sentieon} driver -t ${thread} \
        -i ${align_dir}/${sampleID}.normal.sorted.bam \
        --algo Dedup \
        --score_info ${align_dir}/${sampleID}.normal.score.txt \
        --metrics ${align_dir}/${sampleID}.normal.dedup_metrics.txt \
        ${align_dir}/${sampleID}.normal.sorted.dedup.bam;
    fi


    # step4 - quality control
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- qc";

    if [[ ! -d $qc_dir/${sampleID}.normal ]]; then
        mkdir $qc_dir/${sampleID}.normal
    fi

    if [[ ! -d $qc_dir/${sampleID}.tumor ]]; then
        mkdir $qc_dir/${sampleID}.tumor
    fi

    if [[ ! -d $qc_dir/${sampleID}.normal.loci ]]; then
        mkdir $qc_dir/${sampleID}.normal.loci
    fi

    if [[ ! -d $qc_dir/${sampleID}.tumor.loci ]]; then
        mkdir $qc_dir/${sampleID}.tumor.loci
    fi

    if [[ ! -d $qc_dir/${sampleID}.normal.genes ]]; then
        mkdir $qc_dir/${sampleID}.normal.genes
    fi

    if [[ ! -d $qc_dir/${sampleID}.tumor.genes ]]; then
        mkdir $qc_dir/${sampleID}.tumor.genes
    fi

    $bamdst -p $bed3 -o $qc_dir/${sampleID}.normal \
    ${align_dir}/${sampleID}.normal.sorted.dedup.bam;

    $bamdst -p $bed3 -o $qc_dir/${sampleID}.tumor \
    ${align_dir}/${sampleID}.tumor.sorted.dedup.bam;

    $bamdst -p $bed1 -o $qc_dir/${sampleID}.normal.loci \
    ${align_dir}/${sampleID}.normal.sorted.dedup.bam;

    $bamdst -p $bed1 -o $qc_dir/${sampleID}.tumor.loci \
    ${align_dir}/${sampleID}.tumor.sorted.dedup.bam;

    $bamdst -p $bed2 -o $qc_dir/${sampleID}.normal.genes \
    ${align_dir}/${sampleID}.normal.sorted.dedup.bam;

    $bamdst -p $bed2 -o $qc_dir/${sampleID}.tumor.genes \
    ${align_dir}/${sampleID}.tumor.sorted.dedup.bam;

    $samtools stats -@ ${thread} ${align_dir}/${sampleID}.normal.sorted.dedup.bam > ${qc_dir}/${sampleID}.normal.stats.txt;
    $samtools stats -@ ${thread} ${align_dir}/${sampleID}.tumor.sorted.dedup.bam > ${qc_dir}/${sampleID}.tumor.stats.txt;

    tumor_r1=$(du $input_folder/${sampleID}_tumor_R1.fastq.gz -shL |awk '{print $1}');
    tumor_r2=$(du $input_folder/${sampleID}_tumor_R2.fastq.gz -shL |awk '{print $1}');
    normal_r1=$(du $input_folder/${sampleID}_normal_R1.fastq.gz -shL |awk '{print $1}');
    normal_r2=$(du $input_folder/${sampleID}_normal_R2.fastq.gz -shL |awk '{print $1}');

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

    normal_mapping_rate=$(grep "Fraction of Mapped Reads" $qc_dir/${sampleID}.normal/coverage.report | awk -F"\t" '{print $2}');
    tumor_mapping_rate=$(grep "Fraction of Mapped Reads" $qc_dir/${sampleID}.tumor/coverage.report | awk -F"\t" '{print $2}');

    normal_dup_rate=$(grep "Fraction of PCR duplicate reads" $qc_dir/${sampleID}.normal/coverage.report |awk -F"\t" '{print $2}');
    tumor_dup_rate=$(grep "Fraction of PCR duplicate reads" $qc_dir/${sampleID}.tumor/coverage.report |awk -F"\t" '{print $2}');

    normal_on_target=$(grep "Fraction of Target Reads in all reads" $qc_dir/${sampleID}.normal/coverage.report |awk -F"\t" '{print $2}');
    tumor_on_target=$(grep "Fraction of Target Reads in all reads" $qc_dir/${sampleID}.tumor/coverage.report |awk -F"\t" '{print $2}');

    normal_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt);
    normal_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt);

    tumor_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt);
    tumor_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt);
    
    normal_100x=$(less -S $qc_dir/${sampleID}.normal.loci/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 100) count+=1} END {print count/NR*100}');
    normal_150x=$(less -S $qc_dir/${sampleID}.normal.loci/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 150) count+=1} END {print count/NR*100}');
    normal_200x=$(less -S $qc_dir/${sampleID}.normal.loci/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 200) count+=1} END {print count/NR*100}');

    normal_300x=$(less -S $qc_dir/${sampleID}.normal.genes/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 300) count+=1} END {print count/NR*100}');
    normal_400x=$(less -S $qc_dir/${sampleID}.normal.genes/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 400) count+=1} END {print count/NR*100}');
    normal_500x=$(less -S $qc_dir/${sampleID}.normal.genes/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 500) count+=1} END {print count/NR*100}');

    tumor_100x=$(less -S $qc_dir/${sampleID}.tumor.loci/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 100) count+=1} END {print count/NR*100}');
    tumor_150x=$(less -S $qc_dir/${sampleID}.tumor.loci/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 150) count+=1} END {print count/NR*100}');
    tumor_200x=$(less -S $qc_dir/${sampleID}.tumor.loci/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 200) count+=1} END {print count/NR*100}');

    tumor_300x=$(less -S $qc_dir/${sampleID}.tumor.genes/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 300) count+=1} END {print count/NR*100}');
    tumor_400x=$(less -S $qc_dir/${sampleID}.tumor.genes/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 400) count+=1} END {print count/NR*100}');
    tumor_500x=$(less -S $qc_dir/${sampleID}.tumor.genes/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 500) count+=1} END {print count/NR*100}');
    
    normal_mean_loci=$(grep "Average depth" $qc_dir/${sampleID}.normal.loci/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
    normal_mean_dedup_loci=$(grep "Average depth(rmdup)" $qc_dir/${sampleID}.normal.loci/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
    normal_01x_loci=$(less -S $qc_dir/${sampleID}.normal.loci/depth.tsv.gz | awk -v depth=${normal_mean_loci} 'BEGIN {count=0} {if ($4 > depth*0.1) count+=1} END {print count/NR*100}');
    normal_02x_loci=$(less -S $qc_dir/${sampleID}.normal.loci/depth.tsv.gz | awk -v depth=${normal_mean_loci} 'BEGIN {count=0} {if ($4 > depth*0.2) count+=1} END {print count/NR*100}');
    normal_05x_loci=$(less -S $qc_dir/${sampleID}.normal.loci/depth.tsv.gz | awk -v depth=${normal_mean_loci} 'BEGIN {count=0} {if ($4 > depth*0.5) count+=1} END {print count/NR*100}');
    normal_1x_loci=$(less -S $qc_dir/${sampleID}.normal.loci/depth.tsv.gz | awk -v depth=${normal_mean_loci} 'BEGIN {count=0} {if ($4 > depth) count+=1} END {print count/NR*100}');
    
    tumor_mean_loci=$(grep "Average depth" $qc_dir/${sampleID}.tumor/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
    tumor_mean_dedup_loci=$(grep "Average depth(rmdup)" $qc_dir/${sampleID}.tumor/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
    tumor_01x_loci=$(less -S $qc_dir/${sampleID}.tumor.loci/depth.tsv.gz | awk -v depth=${tumor_mean_loci} 'BEGIN {count=0} {if ($4 > depth*0.1) count+=1} END {print count/NR*100}');
    tumor_02x_loci=$(less -S $qc_dir/${sampleID}.tumor.loci/depth.tsv.gz | awk -v depth=${tumor_mean_loci} 'BEGIN {count=0} {if ($4 > depth*0.2) count+=1} END {print count/NR*100}');
    tumor_05x_loci=$(less -S $qc_dir/${sampleID}.tumor.loci/depth.tsv.gz | awk -v depth=${tumor_mean_loci} 'BEGIN {count=0} {if ($4 > depth*0.5) count+=1} END {print count/NR*100}');
    tumor_1x_loci=$(less -S $qc_dir/${sampleID}.tumor.loci/depth.tsv.gz | awk -v depth=${tumor_mean_loci} 'BEGIN {count=0} {if ($4 > depth) count+=1} END {print count/NR*100}');
    
    normal_mean_genes=$(grep "Average depth" $qc_dir/${sampleID}.normal.genes/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
    normal_mean_dedup_genes=$(grep "Average depth(rmdup)" $qc_dir/${sampleID}.normal.genes/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
    normal_01x_genes=$(less -S $qc_dir/${sampleID}.normal.genes/depth.tsv.gz | awk -v depth=${normal_mean_genes} 'BEGIN {count=0} {if ($4 > depth*0.1) count+=1} END {print count/NR*100}');
    normal_02x_genes=$(less -S $qc_dir/${sampleID}.normal.genes/depth.tsv.gz | awk -v depth=${normal_mean_genes} 'BEGIN {count=0} {if ($4 > depth*0.2) count+=1} END {print count/NR*100}');
    normal_05x_genes=$(less -S $qc_dir/${sampleID}.normal.genes/depth.tsv.gz | awk -v depth=${normal_mean_genes} 'BEGIN {count=0} {if ($4 > depth*0.5) count+=1} END {print count/NR*100}');
    normal_1x_genes=$(less -S $qc_dir/${sampleID}.normal.genes/depth.tsv.gz | awk -v depth=${normal_mean_genes} 'BEGIN {count=0} {if ($4 > depth) count+=1} END {print count/NR*100}');

    tumor_mean_genes=$(grep "Average depth" $qc_dir/${sampleID}.tumor.genes/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
    tumor_mean_dedup_genes=$(grep "Average depth(rmdup)" $qc_dir/${sampleID}.tumor.genes/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
    tumor_01x_genes=$(less -S $qc_dir/${sampleID}.tumor.genes/depth.tsv.gz | awk -v depth=${tumor_mean_genes} 'BEGIN {count=0} {if ($4 > depth*0.1) count+=1} END {print count/NR*100}');
    tumor_02x_genes=$(less -S $qc_dir/${sampleID}.tumor.genes/depth.tsv.gz | awk -v depth=${tumor_mean_genes} 'BEGIN {count=0} {if ($4 > depth*0.2) count+=1} END {print count/NR*100}');
    tumor_05x_genes=$(less -S $qc_dir/${sampleID}.tumor.genes/depth.tsv.gz | awk -v depth=${tumor_mean_genes} 'BEGIN {count=0} {if ($4 > depth*0.5) count+=1} END {print count/NR*100}');
    tumor_1x_genes=$(less -S $qc_dir/${sampleID}.tumor.genes/depth.tsv.gz | awk -v depth=${tumor_mean_genes} 'BEGIN {count=0} {if ($4 > depth) count+=1} END {print count/NR*100}');

    echo "${sampleID}.normal,${normal_r1}/${normal_r2},${normal_raw_reads},${normal_raw_bases},${normal_clean_reads},${normal_clean_bases},\
    ${normal_qc_rate},${normal_mapping_rate},${normal_on_target},${normal_dup_rate},${normal_mean_loci},${normal_mean_dedup_loci},${normal_mean_genes},${normal_mean_dedup_genes},\
    ${normal_insert_size},${normal_insert_std},${normal_01x_loci},${normal_02x_loci},${normal_05x_loci},${normal_1x_loci},\
    ${normal_01x_genes},${normal_02x_genes},${normal_05x_genes},${normal_1x_genes},${normal_100x},${normal_150x},${normal_200x},${normal_300x},${normal_400x},${normal_500x}" \
    >> ${qc_dir}/QC_summary.csv

    echo "${sampleID}.tumor,${tumor_r1}/${tumor_r2},${tumor_raw_reads},${tumor_raw_bases},${tumor_clean_reads},${tumor_clean_bases},\
    ${tumor_qc_rate},${tumor_mapping_rate},${tumor_on_target},${tumor_dup_rate},${tumor_mean_loci},${tumor_mean_dedup_loci},${tumor_mean_genes},${tumor_mean_dedup_genes},\
    ${tumor_insert_size},${tumor_insert_std},${tumor_01x_loci},${tumor_02x_loci},${tumor_05x_loci},${tumor_1x_loci},\
    ${tumor_01x_genes},${tumor_02x_genes},${tumor_05x_genes},${tumor_1x_genes},${tumor_100x},${tumor_150x},${tumor_200x},${tumor_300x},${tumor_400x},${tumor_500x}" \
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

    # step8 - variant calling 
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

    # step6 - variant calling TNscope
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- variant calling";

    $sentieon driver -t ${thread} -r ${ref} \
    -i $tumor_bam -q ${align_dir}/${sampleID}.tumor.recal.table \
    -i $normal_bam -q ${align_dir}/${sampleID}.normal.recal.table \
    --algo TNscope \
    --tumor_sample ${sampleID}.tumor --normal_sample ${sampleID}.normal \
    --dbsnp ${dbsnp} ${snv_dir}/${sampleID}.raw.vcf
    

    # step7 - normalisation + remove low quality variants
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- normalise VCF + filter low-support";
    
    # remove SV
    grep -v "SVTYPE=BND" ${snv_dir}/${sampleID}.raw.vcf \
    >  ${snv_dir}/${sampleID}.step1_snv.vcf;

    # split multiallelic sites + left-alignment
    $bcftools norm -m -both -f ${ref} \
    ${snv_dir}/${sampleID}.step1_snv.vcf \
    -o ${snv_dir}/${sampleID}.step2_norm.vcf;

    # bgzip and make index file
    $bgzip ${snv_dir}/${sampleID}.step2_norm.vcf;
    $tabix -p vcf ${snv_dir}/${sampleID}.step2_norm.vcf.gz;

    # filter off-target variants
    $bcftools view -R ${bed2} \
    ${snv_dir}/${sampleID}.step2_norm.vcf.gz \
    > ${snv_dir}/${sampleID}.step3_on_target.vcf;

    # filter low-support variants
    grep "#\|PASS" ${snv_dir}/${sampleID}.step3_on_target.vcf > \
    ${snv_dir}/${sampleID}.step4_filter.vcf

    $bcftools filter -i "(FORMAT/AF[0]) >= 0.05" \
    ${snv_dir}/${sampleID}.step4_filter.vcf > \
    ${snv_dir}/${sampleID}.step5_filter.vcf;

    $bcftools filter -i "(FORMAT/AD[0:0]+AD[0:1]) >= 50" \
    ${snv_dir}/${sampleID}.step5_filter.vcf > \
    ${snv_dir}/${sampleID}.step6_filter.vcf; 

    # merge MNV
    python3 $merge_mnv ${snv_dir}/${sampleID}.step6_filter.vcf ${ref} \
    -o $snv_dir/${sampleID}.step7_MNV_merged.vcf;

    # step8 - annotation
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- VEP Annotation";

    $vep --dir ${vep_dir} --cache --offline --cache_version ${cache_version} \
    --assembly GRCh37 --format vcf --fa ${ref} --force_overwrite --vcf \
    --gene_phenotype --use_given_ref --refseq --check_existing \
    --hgvs --hgvsg --transcript_version --max_af \
    --vcf_info_field ANN -i ${snv_dir}/${sampleID}.step7_MNV_merged.vcf \
    -o ${snv_dir}/${sampleID}.step8_anno.vcf;
    
    python3 $anno_hgvs ${snv_dir}/${sampleID}.step8_anno.vcf \
    $clinic_transcripts $refflat -o ${snv_dir}/${sampleID}.step9_anno.vcf;
    
done

cat $hrd_dir/*csv | awk 'NR==1 || !/(HRD-sum)/' > $hrd_dir/HRD_results.csv

echo "LOGGING: `date --rfc-3339=seconds` -- Analysis finished";