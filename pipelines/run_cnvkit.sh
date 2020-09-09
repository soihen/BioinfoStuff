#!/usr/bin/bash

# ----------------------------------- Description ------------------------------------- #
# Perform:                                                                              #
#   CNV calling using cnvkit                                                            #
#                                                                                       #
# ------------------------------------------------------------------------------------- #
# Usage:                                                                                #
#   [admin@kai]$bash run_cnvkit.sh [input_folder] [output_folder] #
#                                                                                       #
# Depends on existance of normal samples, CNVkit can be run in 'flat' mode (no normal)  #
# or in 'normal' mode (has normals).                                                    #
# In flat mode, you only need to provide two parameters: [input_folder] where all tumor #
# bams are located; and a [output_folder];                                              #
# In normal mode, you need to provide an additional [input_normal_folder] where all     #
# normal bams are located.                                                              #
#                                                                                       #
# You also need to specify the capture method used, i.e. 'hybrid' capture & 'amplicon'  #
# capture. See CNVkit doc for details                                                   #
# ------------------------------------------------------------------------------------- #



# --------------------- set parameters ------------------ #

# tumor/normal matched mode or tumor-only mode (matched || single)
mode="matched"

# whether perform deduplicate step (true || false)
dedup=true

thread=16

# path to cnvkit package
cnvkit="/public/software/cnvkit/"
sentieon="/data/ngs/softs/sentieon/sentieon-genomics-201808.08/bin/sentieon"
bcftools="/public/software/bcftools-1.9/bcftools"
sentieon_license="192.168.1.186:8990"

# additional files
ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
refflat="/data/ngs/database/soft_database/GATK_Resource_Bundle/refFlat.txt"
bed="/public/shared_data/HRD测试数据/HRDinsert.bed"
pon=""
# ------------------------------------------------------- #


input_folder=$1
output_folder=$2

if [[ ! -d $output_folder ]]; 
then
    mkdir $output_folder
fi


export SENTIEON_LICENSE=${sentieon_license};


# step0 -- prepare BED file;
echo "LOGGING: `date --rfc-3339=seconds` -- prepare (anti-)target BED files"

cnv_dir=$output_folder/cnvkit/;
if [[  ! -d $cnv_dir  ]]; then
    mkdir $cnv_dir
fi

${cnvkit}/cnvkit.py target ${bed} --annotate ${refflat} \
--split -o $cnv_dir/target.bed;

${cnvkit}/cnvkit.py antitarget $cnv_dir/target.bed \
-g ${cnvkit}/data/access-5k-mappable.hg19.bed \
-o $cnv_dir/antitarget.bed;


if [[  $mode == "matched"  ]]; then
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


        # step4 - SNV/SNP calling using TNsnv
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- SNV/SNP calling";

        snp_dir=$output_folder/tnsnv/;
        if [[  ! -d $snp_dir  ]]; then
            mkdir $snp_dir
        fi

        if [[ $dedup == true ]]; then
            normal_bam=${align_dir}/${sampleID}.normal.sorted.dedup.bam
            tumor_bam=${align_dir}/${sampleID}.tumor.sorted.dedup.bam
        else
            normal_bam=${align_dir}/${sampleID}.normal.sorted.bam
            tumor_bam=${align_dir}/${sampleID}.tumor.sorted.bam
        fi

        ${sentieon} driver \
        -t ${thread} -r ${ref} \
        -i ${normal_bam} -i ${tumor_bam} \
        --algo TNsnv \
        --tumor_sample ${sampleID}.tumor \
        --normal_sample ${sampleID}.normal \
        ${snp_dir}/${sampleID}.raw.vcf;


        # step5 -- Normalise variants locations; filter low-depth variants
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- variant normalising & filtering";

        $bcftools norm -m -both -f ${ref} \
        ${snp_dir}/${sampleID}.raw.vcf \
        -o ${snp_dir}/${sampleID}.norm.vcf; 

        $bcftools filter -i "(FORMAT/DP[1] >= 30)" \
        ${snp_dir}/${sampleID}.norm.vcf > ${snp_dir}/${sampleID}.filter.vcf;


        if [[  ! -f $pon  ]];
        then
            # prepare baseline
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- prepare normal CNV baseline";

            ${cnvkit}/cnvkit.py coverage ${normal_bam} $cnv_dir/target.bed \
            -p ${thread} -o ${cnv_dir}/${sampleID}.normal.targetcoverage.cnn;
            
            ${cnvkit}/cnvkit.py coverage ${normal_bam} $cnv_dir/antitarget.bed \
            -p ${thread} -o ${cnv_dir}/${sampleID}.normal.antitargetcoverage.cnn;
        fi
    done


    if [[  ! -f $pon  ]];
    then
        # prepare baseline
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- generate CNV baseline";

        ${cnvkit}/cnvkit.py reference ${cnv_dir}/*normal*cnn \
        -f ${ref} -o ${cnv_dir}/baseline.cnn;

        pon=${cnv_dir}/baseline.cnn
    fi
    

    for ifile in $input_folder/*_R1.tumor.fastq.gz;
    do
        # step6 -- CNV calling
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- CNA calling";      

        ${cnvkit}/cnvkit.py coverage ${tumor_bam} $cnv_dir/target.bed \
        -p ${thread} -o ${cnv_dir}/${sampleID}.tumor.targetcoverage.cnn;

        ${cnvkit}/cnvkit.py coverage ${tumor_bam} $cnv_dir/antitarget.bed \
        -p ${thread} -o ${cnv_dir}/${sampleID}.tumor.antitargetcoverage.cnn;

        ${cnvkit}/cnvkit.py fix \
        ${cnv_dir}/${sampleID}.tumor.targetcoverage.cnn \
        ${cnv_dir}/${sampleID}.tumor.antitargetcoverage.cnn \
        $pon \
        -o ${cnv_dir}/${sampleID}.cnr;

        ${cnvkit}/cnvkit.py segment \
        ${cnv_dir}/${sampleID}.cnr \
        --drop-low-coverage \
        -p ${thread} --smooth-cbs \
        -v ${snp_dir}/${sampleID}.filter.vcf \
        -i ${sampleID}.tumor \
        -n ${sampleID}.normal \
        -o ${cnv_dir}/${sampleID}.cns;

        ${cnvkit}/cnvkit.py call \
        ${cnv_dir}/${sampleID}.cns \
        --center -o ${cnv_dir}/${sampleID}.call.cns;

        ${cnvkit}/cnvkit.py scatter \
        ${cnv_dir}/${sampleID}.cnr \
        -s ${cnv_dir}/${sampleID}.cns \
        -v ${snp_dir}/${sampleID}.filter.vcf \
        -i ${sampleID}.tumor \
        -n ${sampleID}.normal \
        -o ${cnv_dir}/${sampleID}.scatter.png;
    done

elif [[  $mode == 'single'  ]]; then
    echo
fi


