__doc__ = """
Get regions with off-target reads
N.B. `samtools depth` has to be run in advance to get read depth for all genomic coordinates 
You can use my script "annotate_bed.py" on the output file to see what are those off-target regions.
"""
__date__ = "2019/09/04"
__author__ = "Kai"

from collections import defaultdict
import pandas as pd
import argparse
import os


def output_off_target_bed(samtools_depth, bed_file, ofname, threshold):
    """

    :param samtools_depth:
    :param bed_file:
    :param ofname:
    :return:
    """
    off_target_interval = merge_off_target_region(samtools_depth, bed_file, threshold)
    with open(ofname, "w") as fh:
        for each_interval in off_target_interval:
            each_line = "{}\t{}\t{}\n".format(each_interval[0], each_interval[1], each_interval[2])
            fh.write(each_line)


def merge_off_target_region(samtools_depth, bed_file, threshold):
    """

    :param samtools_depth:
    :param bed_file:
    :return:
    """
    off_target_region = get_off_target_region(samtools_depth, bed_file, threshold)
    bed_interval = set()
    for chrom in off_target_region:
        bases = off_target_region[chrom]
        start, end = bases[0], bases[0]
        for index, pos in enumerate(bases):
            if index == 0:
                continue
            if pos == end + 1:
                end = pos
            elif pos > end + 1:
                bed_interval.add((chrom, start, end))
                start = pos
                end = pos
    return bed_interval


def get_off_target_region(samtools_depth, bed_file, threshold):
    """

    :param samtools_depth:
    :param threshold:
    :return:
    """
    bed_intervals = read_bed(bed_file)
    off_target_reads = defaultdict(list)
    df = pd.read_csv(samtools_depth, sep="\t", header=None)
    for index, row in df.iterrows():
        # index 0, 1, 2 -- chr, pos, depth
        # ignore bases located inside target regions
        if row[0] in bed_intervals:
            if row[1] in bed_intervals[row[0]]:
                continue
        if row[2] < threshold:
            continue
        off_target_reads[row[0]].append(row[1])
    return off_target_reads


def read_bed(bed_file):
    """

    :param bed_file:
    :return:
    """
    bed_interval = {}
    with open(bed_file, "r") as fh:
        for line in fh:
            line = line.strip().split()
            bed_interval[line[0]] = range(int(line[1]), int(line[2])+1)
    return bed_interval


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", help="output from samtools depth", required=True)
    parser.add_argument("-b", "--bed", help="target BED file", required=True)
    parser.add_argument("-o", "--ofname", help="output filename")
    parser.add_argument("-t", "--threshold", help="threshold for minimal depth coverage", default=100)
    args = parser.parse_args()

    input_fname = os.path.abspath(args.input)
    if not args.ofname:
        ofname = os.path.join(os.path.dirname(input_fname), os.path.basename(input_fname).split(".")[0] + ".off_target.bed")
    else:
        ofname = args.ofname
    output_off_target_bed(input_fname, args.bed, ofname, args.threshold)