# Tertiary Analysis (unfinished)

As next-generation sequencing becoming more and more common in both clinics and medical research, interpretation of those sequencing data into biological and clinical useful information becomes the key to unraveling genome mystery. Generally, NGS data analysis usually consists of three major steps: 1) primary analysis, which converted raw sequencing signal into digital data (fastq file); 2) secondary analysis, where the most computationally intensive works occurs, perform alignment and variant calling, thus converted fastq file into varaint calling format file (VCF); and 3) tertiary analysis, interpretation of the VCF file, which involves biological classification of identified variants, determination of clinical relevance and provide actionable clinical advices such as drugs or therapies. 

This note will briefly address some key steps involved in tertiary analysis in **cancer**, and provide some plausible solution to some of issues.

## Overview

Tertiary Analysis involves two key questions: **classification of variants** and **representation of variants**. Unlike in Mendelian diseases , where variants are classified based on pathogenicity ([Richards et al., 2015](https://www.nature.com/articles/gim201530)), interpretation of cancer somatic variation should be focused on their impact on clinical care ([Li et al., 2017](https://www.ncbi.nlm.nih.gov/pubmed/27993330/)). i.e. four-tier classification: Variants of strong clinical significance, variants of potential clinical significance, variants of unknown clinical significance and benign or likely benign variants. Shown in figure below: 

![image-20191219093850372](https://github.com/ZKai0801/BioinfoStuff/blob/figures/figures/cancer_variants_classification.png?raw=true)

As you may noticed, the tier IV seems a little bit discordant with other three tiers. The word "benign" is actually a description of variants' pathogenic feature. And that is exactly how the guideline told us. 

> Categorization and interpretation of tier IV variants may refer to recently published ACMG/AMP standards and guidelines for the interpretation of
> germline sequence variants.

Below is the criteria for classifying benign variants.

![image-20191219092415015](https://github.com/ZKai0801/BioinfoStuff/blob/figures/figures/benign_classification.png?raw=true)

So the overall logic procedures for tertiary analysis here involves three steps. Firstly, normalise VCF file and filter out any deemed "false positive" variants. (This step actually belongs to secondary analysis) Secondly, merge MNVs and annotate VCF file. And finally, remove all benign variants. 

## Step1: Normalisation of VCF

Variant calling software like Mutect2 generates VCF file that slightly differ from those used in majority of database. To be specific, in two ways:

1. variants info for a same genomic position are placed in one record, shown as in below figure

   ![image-20191213095536400](https://github.com/ZKai0801/BioinfoStuff/blob/figures/figures/multiallelic-sites.png?raw=true)

   In the above figure, allele TC is believed to mutate into both CT and CC allele. Sites like this, is so called multiallelic sites, which should be splitted into two VCF records. 

2. InDels are called based on position of alignment, which may not be unique. 

   (See indel normalisation section)

```bash
# split multiallelic sites
bcftools norm -m -both ${sampleID}.twicefiltered.vcf --threads ${threads} -o ${sampleID}.somatic.raw.split.vcf

# shift indel to the leftmost position
bcftools norm -f ${hg19fa} ${sampleID}.somatic.raw.split.vcf --threads ${threads} -o ${sampleID}.somatic.split.norm.vcf

```

## Step2: Filter false positive variants

To confidently call a variant, people usually filter variants based on three points:

- a minimal allelic depth
- a minimal allele frequency 
- sites listed in in-house built blacklist (variants that frequently seen in multiple sample, maybe caused by problem in panel)

## Step3: merge MNVs

MNV (multi-nucleotide variant) is defined as two or more nearby variants existing on the same haplotype in an individual, is a clinically and biologically important class of genetic variation. (Wang et al., 2019)

![image2](https://github.com/ZKai0801/BioinfoStuff/blob/figures/figures/mnv.png?raw=true)

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

Therefore, it is essential to correctly identify and annotate MNVs. GATK, contains phasing information on its output file, which provides a convenient way to identifying MNVs. For further details, please see script [merge_mnv.py](https://github.com/ZKai0801/BioinfoStuff/blob/master/tertiary_analysis/merge_mnv.py)  in this package.

## Step4: Annotation

Bear in mind that what we really want to do in annotation step: 1) understand clinical significance of each variant and 2) understand benignity of each variant. Of course, the more annotation softwares/databases we used, the more likely we can get our desired information. But on the hand, more tools/database we used, it more likely we get discrepancy between each tool/databases, and increasingly difficulty to utilise those information. In order to solve these discordance, we need to prioritise variants with clinical significance, then look at rest of variants and tell which one is VUS or benign. Information regarding each variant should also be categorize into different categories. For example, no matter how many functional studies can be found to support one variant is benign, this can only be treated as one score for BS4 (see last section for details).

Currently, there are three major annotation softwares: Annovar, ensembl-VEP and snpEFF

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
vep  --dir ${vep_dir} --cache --offline --cache_version 98  --refseq --assembly GRCh37 --fa ${fasta} --force_overwrite --vcf --variant_class --gene_phenotype --vcf_info_field ANN --hgvs --hgvsg --transcript_version -i ${input} -o ${output}
```

## Key problem1: HGVS annotation

HGVS nomenclature is a system that describes sequence variation. This website provide a way to correct HGVS annotation: https://mutalyzer.nl/. As HGVS system not only annotates variant based on genomic coordinate, but also in resulted transcripts changes and amino acid sequences changes, the correct annotation on which transcript become vital. After tried a number of annotation  tools, Ensembl-VEP (version 98.3 and later version) can correctly annotate HGVS name to variants. (Well, I can't say one hundred percent correct, but so far they all looks fine) Below section shows a key problem I encountered during HGVS annotation. 

## Key problem2: Indel Normalisation

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

```bash
vep --dir <dataset> --cache --offline --cache_version 98 --use_given_ref --refseq --assembly GRCh37  --fa  <reference_fasta> --force_overwrite --vcf --variant_class --gene_phenotype --vcf_info_field ANN --hgvs --hgvsg --transcript_version -i input.vcf -o output.vcf 
```

## Step5: Remove benign variants

Filter out benign variants 

![benign_classification2](https://github.com/ZKai0801/BioinfoStuff/blob/figures/figures/benign_classification2.png?raw=true)

But as we can see from figure 2, not many criteria can be used to determine the benign status of a somatic variants. 

1. BA -- population freq > 1% in ESP, 1000 genome project or ExAC (dbSNP use gnomAD freq)
2. BS1 -- population freq > freq in rare disease (orphatNet) by 1%
3. BS2 -- For early on-set cancer, if the variants found in 1000g project, then viewed them as benign
4. BS3 -- build an internal database that show no damaging effect made by the variant
5. BS4 --  for hereditary cancer, no segregation in affected members of a family
6. BP1 -- missense variant in a gene (tumor suppressor genes) for which primarily truncating variants are known to cause disease
7. BP2 -- dropped
8. BP3 -- In-frame deletions/insertions in a repetitive region without a known function
9. BP4 -- multiple computational prediction software shows no impact
10. BP5 -- dropped
11. BP6 -- reputable sources (database such as clinvar, HGMD) show them as benign
12. BP7 -- synonymous variant that are predicted have no impact on splicing


