#！/usr/bin/bash
# ------------------------------------------ #
#             Factera Pipeline
#
# run bwa-aln + factera for all samples
# under a provided folder 
#
# Usage:
# [admin@kai]$ ./run_factera [input_folder]
#
# ------------------------------------------ #
# set parameters
samtools="/data/ngs/softs/samtools/samtools"
factera="/data/ngs/softs/factera/factera.pl"
bwa="/data/ngs/softs/bwa-0.7.17/bwa"
hg19="/data/ngs/database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
hg19bit="/data/ngs/database/GATK_Resource_Bundle/hg19/hg19.2bit"
exons="/data/ngs/softs/factera/exons.bed"
bed="/data/ngs/database/bed/"
thread="20"
# ------------------------------------------- #
all_files=`find $1 -name *.fusion.json`
wkdir=`pwd`
for ifile in $all_files
do
	sample_path=`dirname $ifile`
	sampleID=`basename $sample_path`
	mkdir $wkdir/$sample_path/factera_fusion
	cd $wkdir/$sample_path/factera_fusion
	echo "factera starts at: `date`" >> $wkdir/${sample_path}/${sampleID}.log
	$bwa aln $hg19 -t ${thread} ../$sampleID.R1_001.trim.fq.gz > $sampleID.R1.sai
	$bwa aln $hg19 -t ${thread} ../$sampleID.R2_001.trim.fq.gz > $sampleID.R2.sai
	$bwa sampe $hg19 $sampleID.R1.sai $sampleID.R2.sai \
	../$sampleID.R1_001.trim.fq.gz ../$sampleID.R1_001.trim.fq.gz > $sampleID.aln.sam
	$samtools view -bS $sampleID.aln.sam > $sampleID.aln.bam
	$factera $sampleID.aln.bam ${exons} ${hg19bit} ${bed}/${sampleID: -2}.bed
	rm $sampleID.aln.sam $sampleID.R1.sai $sampleID.R2.sai
	echo "factera ends at: `date`" >> $wkdir/${sample_path}/${sampleID}.log
done
