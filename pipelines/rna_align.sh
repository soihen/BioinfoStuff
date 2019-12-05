#/usr/bin/bash

######################################
#            RNA alignment           #
#       fastp + star + samtools      #
######################################


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ];
then
    echo "* You MUST provide three arguments!"
    echo "Usage:"
    echo "bash rna-align.sh [input_folder] [reference_fasta] [GTF file]"
    echo " --------- "
    echo "The input_folder could contains many paired fastq reads"
    exit 1
fi

all_fqs=`find $1 -name *R1.fastq.gz`;
for read1 in $all_fqs;
do
    if [[ -f ${read1/R1/R2} ]]; then
        basedir=`dirname $read1`
        trimdir=$basedir/trimmed
        if [[ ! -d $trimdir ]]; then
            mkdir $trimdir
        fi
        qcdir=$basedir/qc
        if [[ ! -d $qcdir ]]; then
            mkdir $qcdir
        fi
        trimout1=$trimdir/`basename ${read1/R1/R1_trimmed}`
        trimout2=$trimdir/`basename ${read1/R1/R2_trimmed}`
        echo 
        echo ============================================
        echo "`date` --- start running fastp"
        fastp --in1 $read1 --in2 ${read1/R1/R2} \
        --out1 $trimout1 \
        --out2 $trimout2 \
        -c --length_required 3 --detect_adapter_for_pe -p \
        --thread 16 \
        --html $qcdir/`basename ${read1/_R1.fastq.gz/.html}` \
        --json $qcdir/`basename ${read1/_R1.fastq.gz/.json}`
        
        gunzip $trimout1 $trimout2

        aligndir=$basedir/aligned
        if [ ! -d $aligndir ]; then
            mkdir $aligndir
        fi
        echo
        echo ============================================
        echo "`date` --- start running STAR + samtools"
        if [[ ! -d `dirname $2`/star_index ]]; then
            echo 'build reference index for STAR'
            mkdir `dirname $2`/star_index
            STAR --runMode genomeGenerate \
            --genomeFastaFiles $2 \
            --sjdbGTFfile $3 \
            --genomeDir `dirname $2` \
            --runThreadN 16
        fi
        STAR --genomeDir `dirname $2`/star_index \
        --runThreadN 16 \
        --readFilesIn ${trimout1::-3} ${trimout2::-3} \
        --outFileNamePrefix $aligndir/`basename ${read1/R1.fastq.gz/}` \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard 

        samtools index $aligndir/`basename ${read1/R1.fastq.gz/Aligned.sortedByCoord.out.bam}`
        echo ============================================
    else
        echo "cannot find the paired fastq file for $read1 ！"
    fi
done
