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
bcftools="/public/software/bcftools-1.9/bcftools"
bgzip="/public/software/htslib-1.9/bin/bgzip"
tabix="/public/software/htslib-1.9/bin/tabix"
ref="/public/database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
target_bed="/public/home/kai/各类标准品/WES/Exome_RefSeq_merged_probes_hg19.bed"
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

    # step1 - remove BND variant
    grep -v "SVTYPE=BND" $ifile > ${output_folder}/${sampleID}.step1_snv.vcf;

    # step2 - split multiallelic sites + left-alignment
    $bcftools norm -m -both -f ${ref} \
    ${output_folder}/${sampleID}.step1_snv.vcf \
    -o ${output_folder}/${sampleID}.step2_norm.vcf; 

    # step3 - bgzip and make index file
    $bgzip ${output_folder}/${sampleID}.step2_norm.vcf;
    $tabix -p vcf ${output_folder}/${sampleID}.step2_norm.vcf.gz;

    # step4 - filter off-target variants
    $bcftools view -R $target_bed \
    ${output_folder}/${sampleID}.step2_norm.vcf.gz \
    > ${output_folder}/${sampleID}.step3_on_target.vcf;

    # step5 - filter low-support variants
    grep "#\|PASS" ${output_folder}/${sampleID}.step3_on_target.vcf > \
    ${output_folder}/${sampleID}.step4_filter.vcf

    $bcftools filter -i "(FORMAT/AD[0:0]+AD[0:1]) >= ${minimal_depth}" \
    ${output_folder}/${sampleID}.step4_filter.vcf > \
    ${output_folder}/${sampleID}.step5_filter.vcf;
done 









