# Tertiary Analysis

As next-generation sequencing becoming more and more common in both clinics and medical research, interpretation of those sequencing data into biological and clinical useful information becomes the key to unraveling genome mystery. Generally, NGS data analysis usually consists of three major steps: 1) primary analysis, which converted raw sequencing signal into digital data (fastq file); 2) secondary analysis, where the most computationally intensive works occurs, perform alignment and variant calling, thus converted fastq file into varaint calling format file (VCF); and 3) tertiary analysis, interpretation of the VCF file, which involves biological classification of identified variants, determination of clinical relevance and provide actionable clinical advices such as drugs or therapies. 

This article will briefly address some key steps involved in tertiary analysis in cancer, and provide some plausible solution to some of issues.



## Step1: Normalisation of VCF

Variant calling software like Mutect2 generates VCF file that slightly differ from those used in majority of database. To be specific, in two ways:

1. variants info for a same genomic position are placed in one record, shown as in below figure

   ![image1](https://raw.githubusercontent.com/ZKai0801/BioinfoStuff/figures/figures/multiallelic-sites.png)

   In the above figure, allele TC is believed to mutate into both CT and CC allele. Sites like this, is so called multiallelic sites, which should be splitted into two VCF records. 

2. InDels are called based on position of alignment, which may not be unique. 

   (See indel normalisation section)

```bash
# split multiallelic sites
bcftools norm -m -both ${sampleID}.twicefiltered.vcf --threads ${threads} -o ${sampleID}.somatic.raw.split.vcf

# shift indel to the leftmost position
bcftools norm -f ${hg19fa} ${sampleID}.somatic.raw.split.vcf --threads ${threads} -o ${sampleID}.somatic.split.norm.vcf

```



## Step2: Filter-based annotation

Three major annotation softwares are: Annovar, ensembl-VEP and snpEFF

This blog here compare there three annotation tools in detail:

https://blog.goldenhelix.com/the-sate-of-variant-annotation-a-comparison-of-annovar-snpeff-and-vep/



Here, I used both Annovar and VEP to annotate my vcf files and merge information from them together

Annovar:

```bash
perl table_annovar.pl ${sampleID}.somatic.split.norm.vcf ${humandb} \
--outfile ${sampleID}.anno.somatic --buildver hg19  \
--thread ${threads} \
--protocol refGene,1000g2015aug_all,esp6500siv2_all,exac03nonpsych,cosmic70,clinvar_20191013,intervar_20180118,dbnsfp35a,avsnp150  \
--operation g,f,f,f,f,f,f,f,f  --vcfinput --remove --intronhgvs 50
```



VEP:

```
vep  --dir /data/ngs/ltj_data/vep_grch37_bak --cache --offline --cache_version 98  --refseq --assembly GRCh37  --fa  /data/ngs/database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta --force_overwrite --vcf --variant_class --gene_phenotype --vcf_info_field ANN --hgvs --hgvsg --transcript_version -i kaohe.somatic.split.norm.vcf -o kaohe.somatic.split.norm.vep_refseq_1.vcf
```





## Step3: Filter

Filter out benign variants 

So the key step here is to determine if the variant is benign. 





```bash
python3 snv_filter.py ${sampleID}.anno.somatic.hg19_multianno.vcf 
```



**

**How to get reference sequence based on genomic coordinate ?**

1. R package -- `GenomicRanges::GRanges`
2. Python package -- `pyfaidx` 



This is how samtools make index for Fasta file allowing fast sequence retrieving :

http://www.htslib.org/doc/faidx.html

> |  Colname  |                           Details                            |
> | :-------: | :----------------------------------------------------------: |
> |   NAME    |               Name of this reference sequence                |
> |  LENGTH   |      Total length of this reference sequence, in bases       |
> |  OFFSET   | Offset in the FASTA/FASTQ file of this sequence's first base |
> | LINEBASES |               The number of bases on each line               |
> | LINEWIDTH |   The number of bytes in each line, including the newline    |



Below is a script I wrote to retrieve reference sequence (before I found those two packages... )

```python
__doc__ = """
Use genomic coordinate to extract sequence from FASTA file
The fasta file need to be indexed in advanced using samtools
See "main" for example usage
"""
__author__ = "Zhentian Kai"
__date__ = "2019/11/04"

import pandas as pd
from collections import defaultdict
import os


def extract_seq(genomic_coordinate, reference_genome="hg19.fasta"):
    """
    :param: genomic_coordinate: a tuple, e.g. (chr1, 100, 1000)
    :param: reference_genome: path to the FASTA file
    ** fadix_dict = {'seq_name': (length, offset, linebase), ...}
    """
    fadix = reference_genome + ".fai"
    fadix_dict = read_fadix(fadix)

    if (genomic_coordinate[1] > fadix_dict[genomic_coordinate[0]][0] or 
        genomic_coordinate[2] > fadix_dict[genomic_coordinate[0]][0] or 
        genomic_coordinate[1] < 0 or 
        genomic_coordinate[2] < 0):
        raise ValueError("Genomic coordinate fall outside of FASTA boundaries")
    elif genomic_coordinate[1] > genomic_coordinate[2]:
        raise ValueError("End index MUST be bigger than or equal with start index")

    start_line_index = genomic_coordinate[1] // fadix_dict[genomic_coordinate[0]][2]
    start_index = genomic_coordinate[1] % fadix_dict[genomic_coordinate[0]][2]
    end_line_index = genomic_coordinate[2] // fadix_dict[genomic_coordinate[0]][2]
    end_index = genomic_coordinate[2] % fadix_dict[genomic_coordinate[0]][2]
    
    # ----------- get sequence ------------- #
    if start_line_index == end_line_index:
        with open(reference_genome, "r") as fh:
            # set offset 
            fh.seek(fadix_dict[genomic_coordinate[0]][1])
            for index, line in enumerate(fh):    
                    if index == start_line_index:
                        seq = line[start_index-1: end_index]
                        return seq

    elif end_line_index > start_line_index:
        seq = ''
        with open(reference_genome, "r") as fh:
            # set offset 
            fh.seek(fadix_dict[genomic_coordinate[0]][1])
            for index, line in enumerate(fh):
                    if index == start_line_index:
                        seq += line[start_index-1: fadix_dict[genomic_coordinate[0]][2]]
                    elif start_line_index < index < end_line_index:
                        seq += line[:fadix_dict[genomic_coordinate[0]][2]]
                    elif index == end_line_index:
                        seq += line[:end_index]
                        break
        return seq
    
    else:
        raise ValueError("End index MUST be bigger than or equal with start index")
    


def read_fadix(fadix):
    """
    fadix for FASTA file has five columns,
    see samtools for details:
    http://www.htslib.org/doc/faidx.html
    :param: fadix -- path to the fadix file
    :return: fadix_dict = {'seq_name':(length, offset, linebase), ...}
    """
    fadix_dict = {}
    df = pd.read_csv(fadix, sep="\t", header=None)
    for index, row in df.iterrows():
        fadix_dict[row[0]] = (row[1], row[2], row[3])
    return fadix_dict


if __name__ == "__main__":
    print(extract_seq(("chr11_gl000202_random", 1000, 1500)))
```





## Notes: Indel Normalisation

The position of InDel is not unique in the reference. For example:

```text
ref            : CGTATGATCTA [GCGCGC] TAGCTAGCTAGC
left-alignment : CGTATGATCTA [--GCGC] TAGCTAGCTAGC
right-alignment: CGTATGATCTA [GCGC--] TAGCTAGCTAGC
```



This issue does not only affect on annotation, but also affects varaint calling step. Here is a question posted on GATK forum this year:

> Hello,
> Could you please inform me whether it is necessary to `normalize + left align` INDELs before performing VQSR. Since all the training data sets are normalized and left aligned, how does it effect when the input VCF is not normalized and left aligned ?

Sadly, GATK team haven't solve this issue yet.



Here is a detailed slides addressing this issue:

https://genome.sph.umich.edu/w/images/b/b4/Variant_Calling_and_Filtering_for_INDELs.pdf



In the standards and guidelines for interpretation and reporting of sequence variant made by AMP-ACMG-CAP, they mentioned:

> Although the HGVS system recommends right-aligned (shifting the start
> position of the variant to the right until it is no longer possible to do so) representation of sequence variants, VCF specifications require left-aligned representation.



VCF format also specifiy denotion of indel:

>For InDels, the reference String must include the base before the event (which must be reflected in the POS field). (String, Required).



### Left-alignment

Most of database use left-alignment system, for example:

- avsnp138
- avsnp142
- clinvar_20150330
- 1000g2014oct
- exac03
- esp6500siv2

> The standard convention is to place an indel at the left-most position possible. —GATK

> Since left-normalization is gaining more and more popularity, my suggestion is to just use left-normalization, and that database curators as well as users both use this practice, so that we can compare apples to apples. —ANNOVAR



Therefore, in order to annotate correct information for a given variant, we need to left align all indels:

```bash
bcftools norm -f ${reference} [input] -o [output]
```



### Right-alignment

Right-alignment is much harder to achieve than left-alignment. Left-alignment is simly move indel to the 5' end of reference until it is not possible to do so, whereas in right-alignment, we need to consider which transcript the variant lies on. If the transcript lie on forward strand, then just move indels to 5'; else, move indels to the opposite direction. 

See HGVS annotation section for solution




## Step4: merge MNV 

MNVs (multi-nucleotide variants) defined as two or more nearby variants existing on the same haplotype in an individual, are a clinically and biologically important class of genetic variation. (Wang et al., 2019)

![image2](https://raw.githubusercontent.com/ZKai0801/BioinfoStuff/figures/figures/mnv.png)

 For instance, the two variants depicted in [Fig. 1b](https://www.biorxiv.org/content/10.1101/573378v2.full#F1) are each predicted individually to have missense consequences, but in combination result in a nonsense variant. (Wang et al., 2019)

 

***in cis***  -- variants occur on same haplotype

***in trans*** -- variants occur on different haplotype



**Phasing methods:**

1. read-based phasing

   assesses whether nearby variants co-segregate on the same reads in DNA sequencing data

2. family-based phasing

   assesses whether pairs of variants are co-inherited within families

3. population-based phasing

   leverages haplotype sharing between members of a large genotyped population to make a statistical inference of phase



See script [merge_mnv.py](https://github.com/ZKai0801/BioinfoStuff/blob/master/tertiary_analysis/merge_mnv.py) for solution





## Step5: HGVS annotation

https://mutalyzer.nl/





## Step6: Annotation of clinical significance



```bash
# 本脚本注释三种信息：
# LRG.record -- LRG号
# LRG.transcript -- LRG号所对应的转录本号
# LRG.strand -- LRG转录本的正负链信息
python3 annotate_lrg.py [input_vcf] [parsed_LRG_database]
```







