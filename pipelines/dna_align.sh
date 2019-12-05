#/usr/bin/bash

######################################
#            DNA alignment           #
#     fastp + bwa-mem + samtools     #
######################################

if [ -z "$1" ] || [ -z "$2" ]
then
    echo "* You MUST provide two arguments!"
    echo "Usage:"
    echo "bash dna-align.sh [input_folder] [reference_fasta]"
    echo " --------- "
    echo "The input_folder could contains many paired fastq reads"
    exit 1
fi

all_fqs=`find $1 -name *_R1.fastq.gz`;
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
        echo 
        echo ============================================
        echo "`date` --- start running fastp"
        fastp --in1 $read1 --in2 ${read1/R1/R2} \
        --out1 $trimdir/`basename ${read1/R1/R1_trimmed}` \
        --out2 $trimdir/`basename ${read1/R1/R2_trimmed}` \
        -c --length_required 3 --detect_adapter_for_pe -p \
        --thread 16 \
        --html $qcdir/`basename ${read1/_R1.fastq.gz/.html}` \
        --json $qcdir/`basename ${read1/_R1.fastq.gz/.json}`

        aligndir=$basedir/aligned
        if [ ! -d $aligndir ]; then
            mkdir $aligndir
        fi
        echo
        echo ============================================
        echo "`date` --- start running bwa-mem + samtools"
        if [[ ! -f $2.bwt ]]; then
            echo 'build reference index for BWA'
            bwa index $2
        fi
        file_name=`basename $read1`
        id=${file_name%%"_R"*}
        bwa mem -t 16 -M \
        -R "@RG\\tID:$id\\tSM:${id}\\tPL:ILLUMINA\\tLB:library" \
        $2 \
        $trimdir/`basename ${read1/R1/R1_trimmed}` \
        $trimdir/`basename ${read1/R1/R2_trimmed}` \
        | samtools view -bS | samtools sort -m 2G -o $aligndir/`basename ${read1/_R1.fastq.gz/.sorted.bam}`
        samtools index $aligndir/`basename ${read1/_R1.fastq.gz/.sorted.bam}`
        samtools stats $aligndir/`basename ${read1/_R1.fastq.gz/.sorted.bam}` > $aligndir/`basename ${read1/_R1.fastq.gz/.stats}`
        echo ============================================
    else
        echo "cannot find the paired fastq file for $read1 ！"
    fi
done
