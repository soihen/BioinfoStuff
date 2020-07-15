#! /usr/bin/bash


# ----------------------------------- Description ------------------------------------- #
# Dependencies:                                                                         #
#       bcftools + bgzip + tabix                                                        #
#                                                                                       #
# Perform:                                                                              #
#     split multiallelic sites +                                                        #
#     left-alignment +                                                                  #
#     filter off-target variants +                                                      #
#     filter low-support variants                                                       #
# ------------------------------------------------------------------------------------- #

# --------------------------------- set parameters ------------------------------------ #
bcftools="/usr/local/bin/bcftools"
bgzip="/usr/local/bin/bgzip"
tabix="/usr/local/bin/tabix"
ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
target_bed="/data/ngs/database/bed/DF.bed"
minimal_depth="300"
# ------------------------------------------------------------------------------------- #


input_folder=$1
output_folder=$2


if [[ ! -d $output_folder ]]; 
then
    mkdir $output_folder
fi


for ifile in $input_folder/*vcf;
do 
    sampleID=`basename ${ifile%%"."*}`;

    # step1 - split multiallelic sites + left-alignment
    $bcftools norm -m -both -f ${ref} \
    ${ifile} \
    -o ${output_folder}/${sampleID}.step1_norm.vcf;

    # step2 - bgzip and make index file
    $bgzip ${output_folder}/${sampleID}.step1_norm.vcf;
    $tabix -p vcf ${output_folder}/${sampleID}.step1_norm.vcf.gz;

    # step3 - filter off-target variants
    $bcftools view -R $target_bed \
    ${output_folder}/${sampleID}.step1_norm.vcf.gz \
    > ${output_folder}/${sampleID}.step2_on_target.vcf;

    # step4 - filter low-support variants
    grep "#\|PASS" ${output_folder}/${sampleID}.step2_on_target.vcf > \
    ${output_folder}/${sampleID}.step3_filter.vcf

    $bcftools filter -i "(FORMAT/AD[0:0]+AD[0:1]) >= ${minimal_depth}" \
    ${output_folder}/${sampleID}.step3_filter.vcf > \
    ${output_folder}/${sampleID}.step4.vcf;
done









