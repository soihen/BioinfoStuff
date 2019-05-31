# 面试题

面试本身的意义无非是测量出面试者的水平高低，因此一味的追求难题、复杂的题并不符合我们的需求; 呈梯度装难度的题目才有必要

另外，如果我们将查阅资料的能力也记在面试范围内的话我觉得有必要在面试题开头写一句“本测试为开卷”

再者就是，鉴于我们现在脚本都没有用perl实现的，我在想要不要把perl这项去除掉，毕竟专精比广而浅更重要



## Coding (本部分为开卷)

### Linux Bash

1. 获取文件第10行
2. grep 出日志文件中包含 Error 的所有行
3. 利用awk打印文件BED文件中第四列（基因信息）中所有独特的基因
4. 利用awk从BED文件中计算panel的大小
5. 严博你之前说了一个什么来着，我当时觉得挺好，但是现在忘了。。。



### Python

1. 将一段DNA序列转换成其互补序列



### 综合题

MAF文件 是由 TCGA 所定义的一种格式：

> Mutation Annotation Format (MAF) is a tab-delimited text file with aggregated mutation information from [VCF Files](https://docs.gdc.cancer.gov/Data/File_Formats/VCF_Format/) and are generated on a project-level. 

已知：

- MAF文件中每一行为一个变异
- 来自所有样本的变异均以行的形式叠加在一起
- 第16列名为`Tumor_Sample_Barcode，`可用于区分变异是否来自于同一样本

假设我们有一个MAF文件，现在想将突变个数小于5的所有样本都过滤掉，请问可以如何实现？

（Python/R 任选其一）



## Bioinformatics

这方面的题目我觉得比较难出，因为如果前来面试的人之前并没有接触过肿瘤方面的话，如果考察对于知识的掌握就无法准确的反映其个人的能力。毕竟，能力才能决定他能走多远多快，能力强，即使不了解这个行业，也能在段时间内上手

这部分需要闭卷，不然就全google出来了。



1. 请简述二代测序的原理 （考察理论基础）
2. 请简述  fastq, sam, bed, vcf 文件所包含的信息 （考察有多少实际经验）
3. 从拿到原始的测序文件开始，到肿瘤分析，大致可以分为几步、几个方面 （考察对肿瘤的了解）
4. 大致阐述 MSI（microsatelite instability）的表现、产生原因以及检测手段 （考察对于细分领域了解的深浅）





