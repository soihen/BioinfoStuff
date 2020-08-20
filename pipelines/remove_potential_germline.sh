#! /usr/bin/bash


# ------------------ Description --------------------- #
# Dependencies:                                        #
#       vep + rm_common_variant.py                     #
#                                                      #
# Perform:                                             #
#     remove common variant (i.e. MAF >= 0.01 )        #
# ---------------------------------------------------- #

# ------------------------ set parameters ------------------------- #
vep="/data/ngs/database/soft_database/98vep/ensembl-vep/vep"
vep_dir="/data/ngs/database/soft_database/vep_98"
cache_version="98"
ref="/public/database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"

rm_common_variant="/public/home/kai/test_kai/BioinfoStuff/rm_common_variant.py"
# ----------------------------------------------------------------- #


input_folder=$1
output_folder=$2


if [[ ! -d $output_folder ]]; 
then
    mkdir $output_folder
fi


for ifile in $input_folder/*vcf;
do 
    sampleID=`basename ${ifile%%"."*}`;

    # step1 - annotation
    $vep --dir ${vep_dir} --cache --offline --cache_version ${cache_version} \
    -assembly GRCh37 --format vcf --fa ${ref} --force_overwrite --vcf \
    --vcf_info_field ANN -i ${ifile} -o ${output_folder}/${sampleID}.anno.vcf \

    # step2 - remove variants with MAF >= 0.01
    python3 $rm_common_variant ${output_folder}/${sampleID}.anno.vcf > ${output_folder}/${sampleID}.somatic.vcf

done 









