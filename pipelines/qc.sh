#!/usr/bin/bash

# ----------------------------------- Description ------------------------------------- #
# Calculate:                                                                            #
#    1) fastq file sizes                                                                #
#    2) Raw Reads/Bases                                                                 #
#    3) Clean Reads/Bases                                                               #
#    4) Q30 rate                                                                        #
#    5) Mapping rate                                                                    #
#    6) on-target rate (by reads)                                                       #
#    7) Mean depth & Mean deduplicated depth                                            #
#    8) Duplicate rate                                                                  #
#    9) Insert size / standard deviation                                                #
#    10) Uniformality rate (> 0.1x; 0.2x; 0.5x; 1x; > 100x; 200x; 500x)                 #
# ------------------------------------------------------------------------------------- #


# --------------------------- set parameters --------------------------- #
# ---------------------------------------------------------------------- #
thread=8

# software
samtools="/public/software/samtools-1.9/samtools"
bamdst="/public/software/bamdst/bamdst"

# ------------------------------ argparser ----------------------------- #
# ---------------------------------------------------------------------- #
input_folder=$1
output_folder=$2
qc_dir=$3
bed=$4

if [[  $1 == '-h'  ]]; 
then
    echo "Usage: ./qc.sh [input_folder] [results_folder] [qc_folder] [BED file]"
    echo "-------------------------------------------------------------------------"
    echo "suitable for other pipelines as well"
    echo "[input_folder] should contain fastq files with following naming system:"
    echo "\${sampleID}_R[1|2].fastq.gz"
    exit 0
fi


if [[ ! -d $input_folder  ]];
then
    echo "Error: input_folder does not Found!"
    exit 1
fi

if [[ ! -d $output_folder  ]];
then
    echo "Error: results_folder does not Found!"
    exit 1
fi

if [[ ! -f $bed  ]];
then
    echo "Error: BED file does not Found!"
    exit 1
fi

# ----------------------  orgnise output dir  -------------------------- #
# ---------------------------------------------------------------------- #

trim_dir=$output_folder/trim/;
align_dir=$output_folder/align/;

if [[ ! -d $qc_dir ]]; then
    mkdir $qc_dir;
fi

# ---------------------------  LOGGING  -------------------------------- #
# ---------------------------------------------------------------------- #
echo "LOGGING: `date --rfc-3339=seconds` -- Analysis started"
echo "========================================================"
echo "LOGGING: -- settings -- input folder -- ${input_folder}"
echo "LOGGING: -- settings -- results folder -- ${output_folder}"
echo "LOGGING: -- settings -- QC folder -- ${qc_dir}"
echo "LOGGING: -- settings -- BED file -- ${bed}"
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
> $qc_dir/QC_summary.csv;


# ---------------------------------------------------------------------- #
# ---------------------------  Pipeline  ------------------------------- #
# ---------------------------------------------------------------------- #
for ifile in $input_folder/*_R1.fastq.gz;
do 
    sampleID=`basename ${ifile%%"_R1"*}`;

    if [[ ! -d $qc_dir/${sampleID} ]]; then
        mkdir $qc_dir/${sampleID};
    fi;
    
    $bamdst -p $bed -o $qc_dir/${sampleID} \
    ${align_dir}/${sampleID}.sorted.dedup.bam;

    $samtools stats -@ ${thread} ${align_dir}/${sampleID}.sorted.dedup.bam > ${qc_dir}/${sampleID}.stats.txt;

    tumor_r1=$(du $input_folder/${sampleID}_R1.fastq.gz -shL |awk '{print $1}');
    tumor_r2=$(du $input_folder/${sampleID}_R2.fastq.gz -shL |awk '{print $1}');

    tumor_raw_reads=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
    print(fh['summary']['before_filtering']['total_reads'])"`;

    tumor_clean_reads=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
    print(fh['summary']['after_filtering']['total_reads'])"`;

    tumor_raw_bases=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
    print(fh['summary']['before_filtering']['total_bases'])"`;

    tumor_clean_bases=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
    print(fh['summary']['after_filtering']['total_bases'])"`;

    tumor_qc_rate=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
    print(fh['summary']['before_filtering']['q30_rate'])"`;

    tumor_mapping_rate=$(grep "Fraction of Mapped Reads" $qc_dir/${sampleID}/coverage.report | awk -F"\t" '{print $2}');
    
    tumor_mean_depth=$(grep "Average depth" $qc_dir/${sampleID}/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
    tumor_mean_dedup_depth=$(grep "Average depth(rmdup)" $qc_dir/${sampleID}/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
    tumor_dup_rate=$(grep "Fraction of PCR duplicate reads" $qc_dir/${sampleID}/coverage.report |awk -F"\t" '{print $2}');

    tumor_on_target=$(grep "Fraction of Target Reads in all reads" $qc_dir/${sampleID}/coverage.report |awk -F"\t" '{print $2}');

    tumor_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.stats.txt);
    tumor_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.stats.txt);

    tumor_50x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 50) count+=1} END {print count/NR*100}');
    tumor_100x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 100) count+=1} END {print count/NR*100}');
    tumor_150x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 150) count+=1} END {print count/NR*100}');
    tumor_200x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 200) count+=1} END {print count/NR*100}');
    tumor_300x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 300) count+=1} END {print count/NR*100}');
    tumor_400x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 400) count+=1} END {print count/NR*100}');
    tumor_500x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 500) count+=1} END {print count/NR*100}');

    tumor_01x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.1) count+=1} END {print count/NR*100}');
    tumor_02x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.2) count+=1} END {print count/NR*100}');
    tumor_05x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.5) count+=1} END {print count/NR*100}');
    tumor_1x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth) count+=1} END {print count/NR*100}');

    echo "${sampleID},${tumor_r1}/${tumor_r2},${tumor_raw_reads},${tumor_raw_bases},${tumor_clean_reads},${tumor_clean_bases},\
    ${tumor_qc_rate},${tumor_mapping_rate},${tumor_on_target},${tumor_mean_depth},${tumor_mean_dedup_depth},${tumor_dup_rate},\
    ${tumor_insert_size},${tumor_insert_std},${tumor_01x},${tumor_02x},${tumor_05x},${tumor_1x},\
    ${tumor_50x},${tumor_100x},${tumor_150x},${tumor_200x},${tumor_300x},${tumor_400x},${tumor_500x}" \
    >> ${qc_dir}/QC_summary.csv;
done
