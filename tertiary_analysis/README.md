# Tertiary Analysis



## Normalisation

```bash
# 拆分 multiallelic sites
bcftools norm -m -both ${sampleID}.twicefiltered.vcf --threads ${threads} -o ${sampleID}.somatic.raw.split.vcf

# 左移
bcftools norm -f ${hg19fa} ${sampleID}.somatic.raw.split.vcf --threads ${threads} -o ${sampleID}.somatic.split.norm.vcf

```



## Filter-based annotation



```bash
perl table_annovar.pl ${sampleID}.somatic.split.norm.vcf ${humandb} \
--outfile ${sampleID}.anno.somatic --buildver hg19  \
--thread ${threads} \
--protocol refGene,1000g2015aug_all,esp6500siv2_all,exac03nonpsych,cosmic70,clinvar_20191013,intervar_20180118,dbnsfp35a,avsnp150  \
--operation g,f,f,f,f,f,f,f,f  --vcfinput --remove --intronhgvs 50
```





## Filter

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



Below is a script I wrote to retrieve reference sequence (before I knew those two packages.. )

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





## Indel Normalisation

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




## merge MNV 



MNVs (multi-nucleotide variants) defined as two or more nearby variants existing on the same haplotype in an individual, are a clinically and biologically important class of genetic variation. (Wang et al., 2019)



![image-20191211133824859](C:\Users\Topgen191017\AppData\Roaming\Typora\typora-user-images\image-20191211133824859.png)

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





## HGVS annotation

https://mutalyzer.nl/





## Annotation of clinical significance



```bash
# 本脚本注释三种信息：
# LRG.record -- LRG号
# LRG.transcript -- LRG号所对应的转录本号
# LRG.strand -- LRG转录本的正负链信息
python3 annotate_lrg.py [input_vcf] [parsed_LRG_database]
```







