#!/usr/bin env python3
__doc__ = """
annotate CNV vcf files from CNVkit;
and export them into the official format 
"""
__author__ = "Kai"
__date__ = "17/07/2019"

from collections import defaultdict
import pandas as pd
import numpy as np
import argparse
import pysam
import os


def convert_vcf(vcf_file, ofname, refflat):
    """

    :param vcf_file: the vcf file that generated from CNVkit
    :param ofname: output filename 
    :param refflat: refGene.txt -- see https://genome.ucsc.edu/FAQ/FAQformat.html for details
    :return:
    """
    vcf_info, sampleID = parse_cnv_vcf(vcf_file)
    annotated_cnv = compare_refflat(refflat, vcf_info)
    header = ("sampleID", "gene", "transcript", "chromosome", "regions", "type", "copy_number")

    cnv_data = []
    for seg_key in annotated_cnv:
        for cnv_gene in annotated_cnv[seg_key]:
            transcript, chrom, regions, cn = annotated_cnv[seg_key][cnv_gene][0], "chr" + annotated_cnv[seg_key][cnv_gene][2], \
                                             annotated_cnv[seg_key][cnv_gene][3], annotated_cnv[seg_key][cnv_gene][4]
            if cn > 2:
                type = "amplification"
            elif cn == 2:
                type = "neutral"
            else:
                type = "deletion"
            cnv_data.append((sampleID, cnv_gene, transcript, chrom, regions, type, cn))

    df = pd.DataFrame(np.array(cnv_data), columns=header)
    df.to_excel(ofname, index=None)


def compare_refflat(refflat, vcf_info):
    """
    compare CNV segments with refFlat.txt, to find:
    1) harboured gene;
    2) corresponding transcript to that gene (longest)
    3) exons that was influenced by that CNV (possible breakpoint)
    ***
    Noted that an individual gene may harbour more than one CNVs
    i.e. breakpoint occurs inside the gene
    Therefore each CNV segment must be consider separately
    :param refflat:
    :param vcf_info: obtained from parse_cnv_vcf()
    :return:
    """
    # annotation information for each CNV segment are stored separately in annotated_cnv
    # annotated_cnv[seg_key] looks like this:
    #   -- {gene: (transcript, exon_count, chromosome, regions, cn),...}
    annotated_cnv = {}

    with open(refflat, "r") as fh:
        for line in fh:
            line = line.strip().split("\t")

            # skip genes on different chromosomes
            chrom = line[2].lower().replace("chr", "")
            if chrom not in vcf_info:
                continue

            # skip genes that completely not overlapped with CNV segment
            transcript_start, transcript_end = int(line[4]), int(line[5])
            for seg in vcf_info[chrom]:
                seg_key = "{}-{}".format(seg[0], seg[1])

                if seg[1] <= transcript_start or seg[0] >= transcript_end:
                    continue

                gene_name, transcript_name, strand = line[0], line[1], line[3]
                exon_count, exon_starts, exon_ends = int(line[-3]), line[-2], line[-1]

                # only use the longest transcript for annotation
                if seg_key in annotated_cnv:
                    if gene_name in annotated_cnv[seg_key]:
                        if annotated_cnv[seg_key][gene_name][1] >= exon_count:
                            continue

                exon_starts = list(map(int, exon_starts.split(",")[:-1]))
                exon_ends = list(map(int, exon_ends.split(",")[:-1]))
                boundaries = np.array(list(zip(exon_starts, exon_ends)))
                boundaries = list(boundaries.flatten())
                overlaps = find_overlap(seg, boundaries, strand, exon_count)
                if seg_key not in annotated_cnv:
                    annotated_cnv[seg_key] = {}
                annotated_cnv[seg_key][gene_name] = (transcript_name, exon_count, chrom, overlaps, seg[2])
    return annotated_cnv


def find_overlap(seg, boundaries, strand, exon_count):
    """
    find overlapped exons
    :param seg: each CNV segment -- (start, end, CN)
    :param boundaries: boundaries for each exon
    :param strand:
    :param exon_count:
    :return:
    """
    for index, pos in enumerate(boundaries):
        if index == 0:
            if pos >= seg[0]:
                start = 0

        if pos <= seg[0]:
            start = index

        if pos < seg[1]:
            end = index

    # exons are ranked from low to high in "+" strand,
    # but are ranked from high to low in "-" strand
    exon_start = int(start/2) + 1
    exon_end = int(end/2) + 1
    if strand == '-':
        exon_start = exon_count - exon_start + 1
        exon_end = exon_count - exon_end + 1
        exon_start, exon_end = exon_end, exon_start

    overlap = "Exon{}-{}".format(exon_start, exon_end)
    return overlap


def parse_cnv_vcf(vcf_file):
    """
    parse vcf file into dict
    :param vcf_file:
    :return: {'chrom': [(start, end, CN), ...], ...}
    """
    vcf_info = defaultdict(list)
    with pysam.VariantFile(vcf_file, 'r') as vcf:
        sampleID = vcf.header.samples[0]
        for record in vcf:
            if 'CN' in record.samples[sampleID]:
                # chromosome name sometimes startswith Chr, sometimes startswith chr
                vcf_info[record.chrom.lower().replace("chr", "")].append((record.start, record.stop, record.samples[sampleID]['CN']))
    return vcf_info, sampleID


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", help="the cnvkit vcf output")
    parser.add_argument("-r", "--refflat", help="UCSC refFLat file")
    parser.add_argument("-o", "--ofname", help="the output filename")
    args = parser.parse_args()

    if args.ofname:
        ofname = args.ofname
    else:
        ofname = os.path.join(os.path.dirname(args.input), os.path.basename(args.input).split(".")[0] + ".annotated.xlsx")
    convert_vcf(args.input, ofname, args.refflat)
