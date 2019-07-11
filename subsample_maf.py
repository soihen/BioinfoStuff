__doc__ = """simulate panel-sized data by sub-sampling WGS/WES data

Required file:
    MAF file for WGS/WES data
    BED file for a targeted panel
"""
__author__ = "Kai"
__date__ = "11/07/2019"

from collections import defaultdict
import pandas as pd
import argparse
import os


def read_bed(bed):
    """
    read BED file into a list
    :param bed: path to a BED file
    :return: a dictionary --- {chr: [(start, end)...] ... ]
    """
    bedinfo = defaultdict(list)
    with open(bed, 'r') as f:
        for line in f:
            line = line.strip().split()
            chrom = line[0].lower().replace("chr", "")
            bedinfo[chrom].append((int(line[1]), int(line[2])))
    return bedinfo


def subsample_maf(maf, bedinfo, outfile):
    """
    subsample maf data
    :param maf: path to a maf file
    :param bedinfo: a list returned from the function read_bed()
    :param outfile: the output MAF filename
    :return:
    """
    drops = []
    df = pd.read_csv(maf, sep="\t", comment="#", low_memory=False)
    for index, row in df.iterrows():
        chrom = str(row['Chromosome']).lower().replace("chr", "")
        # catch mutants that locate inside bed intervals
        for start, end in bedinfo[chrom]:
            if row['Start_Position'] >= start and row['End_Position'] <= end:
                break
        else:
            drops.append(index)

    df = df.drop(drops)
    df.to_csv(outfile, sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("maf", help="path to the WGS/WES MAF file")
    parser.add_argument("bed", help="path to the simulated BED file")
    parser.add_argument("-o", "--output", help="output filename")
    args = parser.parse_args()

    maf = os.path.abspath(args.maf)
    if args.output:
        output = args.output
    else:
        output = os.path.join(os.path.dirname(maf), os.path.basename(maf).split(".")[0]+".subsampled.maf")

    bedinfo = read_bed(args.bed)
    subsample_maf(maf, bedinfo, output)

