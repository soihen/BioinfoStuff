#!/usr/bin/bash


if [[  $1 == '-h'  ]]; then
    echo "Usage: ./qc.sh [input_folder] [output_folder] [BED file]"
    echo "-------------------------------------------------------------------------"
    echo "suitable for other pipelines as well"
    echo "[input_folder] should contain fastq files with following naming system:"
    echo "  single -- \${sampleID}_R[1|2].fastq.gz"
    echo "  matched -- \${sampleID}_R[1|2].tumor.fastq.gz"
    exit 0
fi


# -------------------- set parameters ------------------- #
mode="matched" # matched || single
thread=16

samtools="/public/software/samtools-1.9/samtools"
bedtools="/public/software/bedtools2/bin/bedtools"


# -------------------------------
# create output directory
input_folder=$1
output_folder=$2
bed=$3


trim_dir=$output_folder/trim/;
align_dir=$output_folder/align/;

qc_dir=$output_folder/qc_gene/;
if [[ ! -d $qc_dir ]]; then
    mkdir $qc_dir
fi

# -------------------------------

echo "sampleID,fastq_size,raw_reads,raw_bases,clean_reads,clean_bases,\
qc30_rate,mapping_rate(%),on-target_percent(%),\
mean_depth,mean_dedup_depth,dup_rate(%),\
average_insert_size,std_insert_size,\
50x_depth_percent(%),100x_depth_percent(%),\
150x_depth_percent(%),200x_depth_percent(%),\
300x_depth_percent(%),400x_depth_percent(%),\
500x_depth_percent(%)" > $qc_dir/QC_summary.csv


if [[ $mode == 'matched'  ]];
then
    for ifile in $input_folder/*_R1.tumor.fastq.gz;
    do
        sampleID=`basename ${ifile%%"_R1"*}`
        
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

        echo "${sampleID}.normal,${normal_r1}/${normal_r2},${normal_raw_reads},${normal_raw_bases},${normal_clean_reads},${normal_clean_bases},${normal_qc_rate},${normal_mapping_rate},${normal_on_target},${normal_mean_depth},${normal_mean_dedup_depth},${normal_dup_rate},${normal_insert_size}, ${normal_insert_std},${normal_50x},${normal_100x},${normal_150x},${normal_200x},${normal_300x},${normal_400x},${normal_500x}" \
        >> $qc_dir/QC_summary.csv

        echo "${sampleID}.tumor,${tumor_r1}/${tumor_r2},${tumor_raw_reads},${tumor_raw_bases},${tumor_clean_reads},${tumor_clean_bases},${tumor_qc_rate},${tumor_mapping_rate},${tumor_on_target},${tumor_mean_depth},${tumor_mean_dedup_depth},${tumor_dup_rate},${tumor_insert_size}, ${tumor_insert_std},${tumor_50x},${tumor_100x},${tumor_150x},${tumor_200x},${tumor_300x},${tumor_400x},${tumor_500x}" \
        >> $qc_dir/QC_summary.csv

    done


elif [[  $mode == 'single'  ]];
then
    for ifile in $input_folder/*_R1.fastq.gz;
    do 
        sampleID=`basename ${ifile%%"_R1"*}`
        
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

        tumor_50x=$(awk 'BEGIN {count=0} {if ($3 > 50) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_100x=$(awk 'BEGIN {count=0} {if ($3 > 100) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_150x=$(awk 'BEGIN {count=0} {if ($3 > 150) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_200x=$(awk 'BEGIN {count=0} {if ($3 > 200) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_300x=$(awk 'BEGIN {count=0} {if ($3 > 300) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_400x=$(awk 'BEGIN {count=0} {if ($3 > 400) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)
        tumor_500x=$(awk 'BEGIN {count=0} {if ($3 > 500) count+=1} END {print count/NR*100}' $qc_dir/${sampleID}.dedup.depth.txt)

        echo "${sampleID},${tumor_r1}/${tumor_r2},${tumor_raw_reads},${tumor_raw_bases},${tumor_clean_reads},${tumor_clean_bases},${tumor_qc_rate},${tumor_mapping_rate},${tumor_on_target},${tumor_mean_depth},${tumor_mean_dedup_depth},${tumor_dup_rate},${tumor_insert_size}, ${tumor_insert_std},${tumor_50x},${tumor_100x},${tumor_150x},${tumor_200x},${tumor_300x},${tumor_400x},${tumor_500x}" \
        >> $qc_dir/QC_summary.csv
    done

fi