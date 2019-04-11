__author__ = 'Kai'
__date__ = '11/04/2019'
__contact__ = 'zhentian.kai@outlook.com'

__doc__ = '''
    Calculate target BED file from bait BED file.

    target BED -- exonic regions of gene of interest
    bait BED -- the actual captured regions,
                usually including 50bp flanking either side of exons of interest
    
    This script will use exon information provided on UCSC refFlat.txt file.
    If the exon is completely encompassed by a baited interval,
    then viewed this exon as a targeted interval.
    If an additional gene_list is provided with -l/--list,
    then only this subset of gene will be included in the targeted BED.
    
    N.B. you may want to sort the BED file afterwards:
    [kai@admin]$ sort -V -k1,1 -k2,2 -k3,3 [xxx.bed]

    Usage:
    [kai@admin]$ python3 bait2target.py -i [bait_BED] -f [refFlat] -o [target_BED] -l [gene_list]

    Required:
        1) bait BED
        2) refFlat.txt
        
    Optional:
        1) gene_list
        2) target_BED
'''

import argparse
import os
import pandas as pd
from collections import defaultdict
import time


def read_bed(bait_bed, annotation = True):
    '''
        read bait BED file and yield tuples: (chr, start, end)...
        The bait BED could be either annotated or not annotated
        @param: annotation -- True or False
        As the information located on the fourth field may not be what we want,
        turing annotation to False provide a way to ignore this information.
    '''
    try:
        df = pd.read_csv(bait_bed, header = None, sep = '\t')
    except:
        raise IOError("The provided bait_bed is in wrong format!")
    if len(df.columns) == 4 and annotation == True:
        for index, row in df.iterrows():
            # keep in mind that row[3], the annotation, could be either empty or contain multiple genes
            # if it is empty, it will be annotated either as a dot, or a empty field
            yield row[0], row[1], row[2], row[3]
    elif len(df.columns) == 3 or annotation == False:
        for index, row in df.iterrows():
            # the reason why I yield 000 on the fourth column is that we need to distinguish
            # the concept of lack of annotation and no gene could be annotated on that interval
            yield row[0], row[1], row[2], '000'



def read_gene_list(gene_list):
    '''
        gene_list should be a plain text file, where each line is a single gene name
    '''
    genelist = set()
    try:
        with open(gene_list, 'r') as fh:
            for line in fh:
                genelist.add(line.strip())
        return genelist
    except:
        raise IOError("The provided gene_list file does not exist!")



def read_refFlat(refFlat, genelist):
    '''
        read gene and corresponding exon coordinates into a dictionary
        if genelist is provided, then only read a subset of refFlat into the dictionary
        @param: genelist -- returned by read_gene_list()
        @return: exons -- a dictionary where each key is a gene; and key is a list of tuples: (chr, start, end)
    '''
    try:
        df = pd.read_csv(refFlat, header = None, sep = '\t')
    except:
        raise IOError("The provided refFlat.txt is in wrong format!")

    dict_exons = defaultdict(set)
    for index, row in df.iterrows():
        # start and end position of exons are separated by comma and has an empty column on the end
        gene, chr, starts, ends = row[0], row[2], row[9].split(',')[:-1], row[10].split(',')[:-1]
        # genelist existed -- selectively add genes that included in the genelist
        if genelist:
            if gene in genelist:
                # making a list of tuples: [(chr, start, end)...]
                for exon in list(zip([chr]*len(starts), starts, ends)):
                    dict_exons[gene].add(exon)
        else:
            for exon in list(zip([chr]*len(starts), starts, ends)):
                dict_exons[gene].add(exon)
    return dict_exons



def find_targeted_interval(chr, start, end, genename, dict_exons, fw):
    '''
        identify the targeted interval from baited interval
        baited interval may be so large that span across multiple gene's exons,
        but gene's exons unlikely have identical coordinates.
        If overlap can be found between exons, we will write all of exons to the output.
        @param: chr, start, end, genename -- yielded from read_bed()
        @param: fw -- output file handler that need to be opened in advance
    '''
    # use the genename to locate exon coordinates if genename is provided
    if genename and genename != '.':
        # genename is directly derived from BED files
        # therefore may have multiple annotation existed
        if ',' in genename:
            for g in genename.split(','):
                # g -- a single gene name
                if g in dict_exons:
                    for exon in dict_exons[g]:
                        # each exon would be (chr, start, end)
                        # The bait region MUST compleletly covers the exonic regions of interests
                        if chr == exon[0] and start <= int(exon[1]) < int(exon[2]) <= end:
                            fw.write('{}\t{}\t{}\t{}\n'.format(exon[0], exon[1], exon[2], g))
        else:
            if genename in dict_exons:
                for exon in dict_exons[genename]:
                    if chr == exon[0] and start <= int(exon[1]) < int(exon[2]) <= end:
                        fw.write('{}\t{}\t{}\t{}\n'.format(exon[0], exon[1], exon[2], genename))

    # lack of annotation
    elif genename == '000':
        # loop over all the dict_exons
        for gene, exons in dict_exons.items():
            for exon in exons:
                if chr == exon[0] and start <= int(exon[1]) < int(exon[2]) <= end:
                    fw.write('{}\t{}\t{}\t{}\n'.format(exon[0], exon[1], exon[2], gene))

    # no exon is expected to be found within baited interval
    elif genename == '.' or genename == '' or genename == None or pd.isnull(genename):
        pass



def wrapper():
    # argparsing...
    parser = argparse.ArgumentParser(description=__doc__, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', help = 'Path to the bait BED file', required = True)
    parser.add_argument('-o', '--output', help = 'The output file name')
    parser.add_argument('-f', '--refFlat', help = 'Path to the refFlat.txt', required = True)
    parser.add_argument('-l', '--list', help = 'Path to a text file containing genes that wanted to be included in the target BED')
    parser.add_argument('-ignore_anno', action = 'store_true', help = 'ignore the annotation information provided on the bait BED if this flag is used')
    args = parser.parse_args()
    # preparing params...
    inputfile = args.input
    if args.output:
        outputfile = args.output
    else:
        outputfile = os.path.join(os.path.dirname(inputfile), os.path.basename(inputfile).split('.')[0]) + '.target.bed'
    if args.list:
        genelist = read_gene_list(args.list)
    else: genelist = []
    dict_exons = read_refFlat(args.refFlat, genelist)
    # wrapping...
    fw = open(outputfile, 'w')
    if args.ignore_anno:
        for chr, start, end, genename in read_bed(args.input, False):
            find_targeted_interval(chr, start, end, genename, dict_exons, fw)
    else:
        for chr, start, end, genename in read_bed(args.input):
            find_targeted_interval(chr, start, end, genename, dict_exons, fw)
    fw.close()


if __name__ == '__main__':
    wrapper()
