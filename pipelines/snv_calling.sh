#! /usr/bin/bash

# ----------------------------------- Description ------------------------------------- #
# Perform:                                                                              #
#   Trimming (fastp) +                                                                  #
#   Alignment (Sentieon-bwa) +                                                          #
#   statistic analysis +                                                                #
#   Deduplication (optional) +                                                          #
#   Variant Calling +                                                                   #
#   Normalisation +                                                                     #
#   Remove low-quality variants +                                                       #
#   VEP annotation +                                                                    #
#   Remove potential germline mutation +                                                #
#   VIC annotation                                                                      #
# ------------------------------------------------------------------------------------- #
# Usage:                                                                                #
#   [admin@kai]$bash snv_calling.sh [input_folder] [output_folder]                      #
#                                                                                       #
# 1. You MUST specifiy whether use tumor-only mode ('single') or tumor/normal mode      #
#    ('match') of somatic calling of TNscope.                                           #
# 2. Depends on DNA capture methods (i.e. targeted amplicon based or hybrid capture     #
#    based), user can choose to either perform or not perform deduplication step.       #
# 3. Germline mutation removal is followed by:                                          #
#    Sukhai MA, Misyura M, Thomas M, Garg S, Zhang T, Stickle N, Virtanen C, Bedard PL, #
#    Siu LL, Smets T, Thijs G, Van Vooren S, Kamel-Reid S, Stockley TL.                 #
#    Somatic Tumor Variant Filtration Strategies to Optimize Tumor-Only Molecular       #
#    Profiling Using Targeted Next-Generation Sequencing Panels.                        #
#    J Mol Diagn. 2019 Mar;21(2):261-273. doi: 10.1016/j.jmoldx.2018.09.008.            #
#    Epub 2018 Dec 19. PMID: 30576869.                                                  #
# ------------------------------------------------------------------------------------- #


# -------------------- set parameters ------------------- #
# mode of somatic calling: (matched || single)
mode="single"

# calling method (TNseq || TNscope)
calling_method="TNscope"

# whether perform deduplicate step (true || false)
dedup=true

# path to software & scripts
fastp="/data/ngs/softs/fastp/fastp"
sentieon="/data/ngs/softs/sentieon/sentieon-genomics-201808.08/bin/sentieon"
bcftools="/public/software/bcftools-1.9/bcftools"
samtools="/public/software/samtools-1.9/samtools"
vep="/public/software/98vep/ensembl-vep/vep"
bgzip="/public/software/htslib-1.9/bin/bgzip"
tabix="/public/software/htslib-1.9/bin/tabix"
VIC="/public/software/VIC"
annovar="/public/software/annovar/"

merge_mnv="/public/home/kai/BioinfoStuff/tertiary_analysis/merge_mnv.py"
anno_hgvs="/public/home/kai/BioinfoStuff/tertiary_analysis/anno_hgvs.py"
rm_common_variant="/public/home/kai/BioinfoStuff/rm_common_variant.py"

sentieon_license="192.168.1.186:8990"
thread=8

# additional files
# reference fasta file
ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
dbsnp="/public/database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf"
k1="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf.gz"
k2="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
k3="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
refflat="/public/database/GATK_Resource_Bundle/refFlat.txt"
vep_dir="/public/software/vep_98/"
cache_version="98"
clinic_transcripts="/public/home/kai/database/LRG/parsed_LRG.tsv"


# ------------------------------------------------------- #


if [[  $1 == '-h'  ]]; then
    echo "Usage: ./snv_calling.sh [input_folder] [output_folder] [BED]"
    echo "-------------------------------------------------------------------------"
    echo "[input_folder] should contain fastq files with following naming system:"
    echo "  single -- \${sampleID}_R[1|2].fastq.gz"
    echo "  matched -- \${sampleID}_normal/tumor_R[1|2].fastq.gz"
    exit 0
fi


# -------------------------------
# create output directory

input_folder=$1
output_folder=$2
bed=$3
cancer_type=$4


if [[ ! -d $input_folder  ]];
then
    echo "Error: input_folder does not Found!"
    exit 1
fi


if [[ ! -f $bed  ]];
then
    echo "Error: BED file does not Found!"
    exit 1
fi


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

snv_dir=${output_folder}/SNV/;
if [[  ! -d $snv_dir  ]];
then 
    mkdir $snv_dir
fi

qc_dir=$output_folder/qc/;
if [[ ! -d $qc_dir ]]; then
    mkdir $qc_dir
fi

# -------------------------------

echo "LOGGING: `date --rfc-3339=seconds` -- Analysis started"
echo "========================================================"
echo "LOGGING: -- settings -- input folder -- ${input_folder}"
echo "LOGGING: -- settings -- output folder -- ${output_folder}"
echo "LOGGING: -- settings -- BED file -- ${bed}"
echo "LOGGING: -- settings -- cancer type -- ${cancer_type}"
echo "========================================================"


echo "sampleID,fastq_size,raw_reads,raw_bases,clean_reads,clean_bases,\
qc30_rate,mapping_rate(%),on-target_percent(%),\
mean_depth,mean_dedup_depth,dup_rate(%),\
average_insert_size,std_insert_size,\
Uniformity_0.1X(%),Uniformity_0.2X(%),\
Uniformity_0.5X(%),Uniformity_1X(%),\
50x_depth_percent(%),100x_depth_percent(%),\
150x_depth_percent(%),200x_depth_percent(%),\
300x_depth_percent(%),400x_depth_percent(%),\
500x_depth_percent(%)" \
> $qc_dir/QC_summary.csv


export SENTIEON_LICENSE=${sentieon_license};


# tumor-normal matched mode
if [[  $mode == 'matched' ]]; then
    for ifile in $input_folder/*_normal_R1.fastq.gz;
    do 
        sampleID=`basename ${ifile%%"_normal"*}`

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

        tumor_r1=$(du $input_folder/${sampleID}_tumor_R1.fastq.gz -sh |awk '{print $1}');
        tumor_r2=$(du $input_folder/${sampleID}_tumor_R2.fastq.gz -sh |awk '{print $1}');
        normal_r1=$(du $input_folder/${sampleID}_normal_R1.fastq.gz -sh |awk '{print $1}');
        normal_r2=$(du $input_folder/${sampleID}_normal_R2.fastq.gz -sh |awk '{print $1}');

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

        normal_01x=$(awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($3 > depth*0.1) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
        normal_02x=$(awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($3 > depth*0.2) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
        normal_05x=$(awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($3 > depth*0.5) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
        normal_1x=$(awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($3 > depth) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)

        normal_50x=$(awk 'BEGIN {count=0} {if ($3 > 50) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
        normal_100x=$(awk 'BEGIN {count=0} {if ($3 > 100) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
        normal_150x=$(awk 'BEGIN {count=0} {if ($3 > 150) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
        normal_200x=$(awk 'BEGIN {count=0} {if ($3 > 200) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
        normal_300x=$(awk 'BEGIN {count=0} {if ($3 > 300) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
        normal_400x=$(awk 'BEGIN {count=0} {if ($3 > 400) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)
        normal_500x=$(awk 'BEGIN {count=0} {if ($3 > 500) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.normal.dedup.depth.txt)

        tumor_01x=$(awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($3 > depth*0.1) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
        tumor_02x=$(awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($3 > depth*0.2) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
        tumor_05x=$(awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($3 > depth*0.5) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
        tumor_1x=$(awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($3 > depth) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)

        tumor_50x=$(awk 'BEGIN {count=0} {if ($3 > 50) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
        tumor_100x=$(awk 'BEGIN {count=0} {if ($3 > 100) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
        tumor_150x=$(awk 'BEGIN {count=0} {if ($3 > 150) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
        tumor_200x=$(awk 'BEGIN {count=0} {if ($3 > 200) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
        tumor_300x=$(awk 'BEGIN {count=0} {if ($3 > 300) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
        tumor_400x=$(awk 'BEGIN {count=0} {if ($3 > 400) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)
        tumor_500x=$(awk 'BEGIN {count=0} {if ($3 > 500) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.tumor.dedup.depth.txt)

        echo "${sampleID}.normal,${normal_r1}/${normal_r2},${normal_raw_reads},${normal_raw_bases},${normal_clean_reads},${normal_clean_bases},${normal_qc_rate},${normal_mapping_rate},${normal_on_target},${normal_mean_depth},${normal_mean_dedup_depth},${normal_dup_rate},${normal_insert_size},${normal_insert_std},${normal_01x},${normal_02x},${normal_05x},${normal_1x},${normal_50x},${normal_100x},${normal_150x},${normal_200x},${normal_300x},${normal_400x},${normal_500x}" \
        >> ${qc_dir}/QC_summary.csv

        echo "${sampleID}.tumor,${tumor_r1}/${tumor_r2},${tumor_raw_reads},${tumor_raw_bases},${tumor_clean_reads},${tumor_clean_bases},${tumor_qc_rate},${tumor_mapping_rate},${tumor_on_target},${tumor_mean_depth},${tumor_mean_dedup_depth},${tumor_dup_rate},${tumor_insert_size},${tumor_insert_std},${tumor_01x},${tumor_02x},${tumor_05x},${tumor_1x},${tumor_50x},${tumor_100x},${tumor_150x},${tumor_200x},${tumor_300x},${tumor_400x},${tumor_500x}" \
        >> ${qc_dir}/QC_summary.csv


        # step5 - Base quality score recalibration 
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
        $bcftools view -R $bed \
        ${snv_dir}/${sampleID}.step2_norm.vcf.gz \
        > ${snv_dir}/${sampleID}.step3_on_target.vcf;

        # filter low-support variants
        grep "#\|PASS" ${snv_dir}/${sampleID}.step3_on_target.vcf > \
        ${snv_dir}/${sampleID}.step4_filter.vcf

        $bcftools filter -i "(FORMAT/AF[0]) >= 0.05" \
        ${snv_dir}/${sampleID}.step4_filter.vcf > \
        ${snv_dir}/${sampleID}.step5_filter.vcf;

        $bcftools filter -i "(FORMAT/AD[0:0]+AD[0:1]) >= 250" \
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

    done

# tumor-only mode 
elif [[  $mode == 'single'  ]]; then
    for ifile in $input_folder/*_R1.fastq.gz;
    do
        sampleID=`basename ${ifile%%"_R1"*}`

        # step1 - trim reads
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- trimming reads";

        $fastp --in1 $input_folder/${sampleID}_R1.fastq.gz \
        --in2 $input_folder/${sampleID}_R2.fastq.gz \
        --out1 $trim_dir/${sampleID}_trim_R1.fastq.gz \
        --out2 $trim_dir/${sampleID}_trim_R2.fastq.gz \
        -c --length_required 3 --detect_adapter_for_pe -p \
        --thread ${thread} \
        --html $trim_dir/${sampleID}.trim.html \
        --json $trim_dir/${sampleID}.trim.json;

        # step2 - align & sort
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- alignment & sorting";

        ($sentieon bwa mem -M -R "@RG\tID:${sampleID}\tSM:${sampleID}\tPL:illumina" \
        -t ${thread} -K 10000000 ${ref} \
        ${trim_dir}/${sampleID}_trim_R1.fastq.gz \
        ${trim_dir}/${sampleID}_trim_R2.fastq.gz \
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
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- deduplication";

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

        # step5 - quality control
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- qc";

        $samtools depth -d 0 -b ${bed} ${align_dir}/${sampleID}.sorted.bam > $qc_dir/${sampleID}.depth.txt;
        $samtools depth -d 0 -b ${bed} ${align_dir}/${sampleID}.sorted.dedup.bam > $qc_dir/${sampleID}.dedup.depth.txt;
        $samtools depth -d 0 ${align_dir}/${sampleID}.sorted.dedup.bam > $qc_dir/${sampleID}.dedup.total_depth.txt;
        $samtools stats -@ ${thread} ${align_dir}/${sampleID}.sorted.dedup.bam > ${qc_dir}/${sampleID}.stats.txt;

        tumor_r1=$(du $input_folder/${sampleID}_R1.fastq.gz -sh |awk '{print $1}');
        tumor_r2=$(du $input_folder/${sampleID}_R2.fastq.gz -sh |awk '{print $1}');

        tumor_raw_reads=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['total_reads'])"`

        tumor_clean_reads=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
        print(fh['summary']['after_filtering']['total_reads'])"`

        tumor_raw_bases=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['total_bases'])"`

        tumor_clean_bases=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
        print(fh['summary']['after_filtering']['total_bases'])"`

        tumor_qc_rate=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
        print(fh['summary']['after_filtering']['q30_rate'])"`

        tumor_mapped_reads=$(awk -F"\t" '$2 == "reads mapped:" {print $3}' ${qc_dir}/${sampleID}.stats.txt);
        tumor_all_reads=$(awk -F"\t" '$2 == "sequences:" {print $3}' ${qc_dir}/${sampleID}.stats.txt);
        tumor_mapping_rate=`bc <<< "scale=4; ${tumor_mapped_reads} / ${tumor_all_reads} * 100"`
        
        tumor_mean_depth=$(awk '{sum+=$3}END{print sum/NR}' $qc_dir/${sampleID}.depth.txt)
        tumor_mean_dedup_depth=$(awk '{sum+=$3}END{print sum/NR}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_dup_rate=`bc <<< "scale=4; (1 - ${tumor_mean_dedup_depth} / ${tumor_mean_depth}) * 100"`

        tumor_total_depth=$(awk '{sum+=$3}END{print sum}' $qc_dir/${sampleID}.dedup.total_depth.txt)
        tumor_on_target=$(awk -v all_bases=${tumor_total_depth} '{sum+=$3}END{print sum/all_bases*100}' $qc_dir/${sampleID}.dedup.depth.txt)

        tumor_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.stats.txt);
        tumor_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.stats.txt);

        tumor_01x=$(awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($3 > depth*0.1) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_02x=$(awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($3 > depth*0.2) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_05x=$(awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($3 > depth*0.5) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_1x=$(awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($3 > depth) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)

        tumor_50x=$(awk 'BEGIN {count=0} {if ($3 > 50) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_100x=$(awk 'BEGIN {count=0} {if ($3 > 100) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_150x=$(awk 'BEGIN {count=0} {if ($3 > 150) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_200x=$(awk 'BEGIN {count=0} {if ($3 > 200) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_300x=$(awk 'BEGIN {count=0} {if ($3 > 300) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_400x=$(awk 'BEGIN {count=0} {if ($3 > 400) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_500x=$(awk 'BEGIN {count=0} {if ($3 > 500) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)

        echo "${sampleID},${tumor_r1}/${tumor_r2},${tumor_raw_reads},${tumor_raw_bases},${tumor_clean_reads},${tumor_clean_bases},${tumor_qc_rate},${tumor_mapping_rate},${tumor_on_target},${tumor_mean_depth},${tumor_mean_dedup_depth},${tumor_dup_rate},${tumor_insert_size},${tumor_insert_std},${tumor_01x},${tumor_02x},${tumor_05x},${tumor_1x},${tumor_50x},${tumor_100x},${tumor_150x},${tumor_200x},${tumor_300x},${tumor_400x},${tumor_500x}" \
        >> $qc_dir/QC_summary.csv

        # set path of the bam file
        if [ $dedup == true ]; then
            bam=${align_dir}/${sampleID}.sorted.dedup.bam
        else
            bam=${align_dir}/${sampleID}.sorted.bam
        fi

        
        # step6 - BQSR
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- BQSR";

        ${sentieon} driver -t ${thread} -r ${ref} \
        -i ${bam} --algo QualCal -k ${k1} -k ${k2} -k ${k3} \
        ${align_dir}/${sampleID}.recal.table;


        # step7 - variant calling
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- variant calling";

        if [[  $calling_method == 'TNscope' ]];
        then
            ${sentieon} driver -t ${thread} -r ${ref} \
            -i ${bam} -q ${align_dir}/${sampleID}.recal.table \
            --algo TNscope --tumor_sample ${sampleID} \
            --dbsnp ${dbsnp} ${snv_dir}/${sampleID}.tumor.raw.vcf;

        elif [[  $calling_method == 'TNseq' ]];
        then
            ${sentieon} driver -t ${thread} -r ${ref} \
            -i ${bam} \
            -q ${align_dir}/${sampleID}.recal.table \
            --algo TNhaplotyper2 --tumor_sample ${sampleID} \
            --default_af 0.00000001 \
            ${snv_dir}/${sampleID}.tumor.tmp.vcf;

            ${sentieon} tnhapfilter --tumor_sample ${sampleID} \
            -v ${snv_dir}/${sampleID}.tumor.tmp.vcf \
            ${snv_dir}/${sampleID}.tumor.raw.vcf;

            rm ${snv_dir}/${sampleID}.tumor.tmp.vcf;
        fi

        # step8 - normalise
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- normalise VCF + filter low-support";

        # remove 'germline_risk' annotation
        $bcftools annotate -x FILTER/germline_risk \
        ${snv_dir}/${sampleID}.tumor.raw.vcf > \
        ${snv_dir}/${sampleID}.step1_deanno.vcf;

        # remove SV
        grep -v "SVTYPE=BND" ${snv_dir}/${sampleID}.step1_deanno.vcf \
        >  ${snv_dir}/${sampleID}.step2_snv.vcf;

        # split multiallelic sites + left-alignment
        $bcftools norm -m -both -f ${ref} \
        ${snv_dir}/${sampleID}.step2_snv.vcf \
        -o ${snv_dir}/${sampleID}.step3_norm.vcf;

        # bgzip and make index file
        $bgzip ${snv_dir}/${sampleID}.step3_norm.vcf;
        $tabix -p vcf ${snv_dir}/${sampleID}.step3_norm.vcf.gz;

        # filter off-target variants
        $bcftools view -R $bed \
        ${snv_dir}/${sampleID}.step3_norm.vcf.gz \
        > ${snv_dir}/${sampleID}.step4_on_target.vcf;

        # filter low-support variants
        grep "#\|PASS" ${snv_dir}/${sampleID}.step4_on_target.vcf > \
        ${snv_dir}/${sampleID}.step5_filter.vcf

        $bcftools filter -i "(FORMAT/AF[0]) >= 0.05" \
        ${snv_dir}/${sampleID}.step5_filter.vcf > \
        ${snv_dir}/${sampleID}.step6_filter.vcf;

        $bcftools filter -i "(FORMAT/AD[0:0]+AD[0:1]) >= 250" \
        ${snv_dir}/${sampleID}.step6_filter.vcf > \
        ${snv_dir}/${sampleID}.step7_filter.vcf; 

        # merge MNVs
        python3 $merge_mnv ${snv_dir}/${sampleID}.step6_filter.vcf ${ref} \
        -o $snv_dir/${sampleID}.step7_MNV_merged.vcf;

        # step9 - annotation
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- VEP Annotation";

        $vep --dir ${vep_dir} --cache --offline --cache_version ${cache_version} \
        --assembly GRCh37 --format vcf --fa ${ref} --force_overwrite --vcf \
        --gene_phenotype --use_given_ref --refseq --check_existing \
        --hgvs --hgvsg --transcript_version --max_af \
        --vcf_info_field ANN -i ${snv_dir}/${sampleID}.step7_MNV_merged.vcf \
        -o ${snv_dir}/${sampleID}.step8_anno.vcf;

        # step 10 - remove potential germline variants
        # 1) MAF >= 0.2% in 1000g, ESP, ExAC
        # 2) benign or likely_benign in ClinVar
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- remove common SNPs";

        python3 $rm_common_variant ${snv_dir}/${sampleID}.step9_anno.vcf > \
        ${snv_dir}/${sampleID}.step10_somatic.vcf;

        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- VIC Annotation";
        
        java -jar ${VIC}/target/VIC-1.0-jar-with-dependencies.jar -cancer_type ${cancer_type} \
        -b hg19 -db ${VIC}/vicdb/ -cosmic cosmic70 -clinvar clinvar_20191202 \
        -gnomad gnomad211_exome -dbnsfp dbnsfp35a \
        -d ${annovar}/humandb/ -table_annovar ${annovar}/table_annovar.pl \
        -convert2annovar ${annovar}/convert2annovar.pl \
        -annotate_variation ${annovar}/annotate_variation.pl \
        -input_type VCF -i ${snv_dir}/${sampleID}.step10_somatic.vcf \
        -o ${snv_dir}/${sampleID}.step11_classified
        
    done

fi

echo "LOGGING: `date --rfc-3339=seconds` -- Analysis finished";
