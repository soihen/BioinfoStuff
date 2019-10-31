__doc__ = """
parse LRG database in advance, to facilitate downstream analysis
output file is a tab-delimited file, arranging like this:
    gene_name   transcript_name
    ...

Usage:
[kai@admin]$ python3 parse_LRG_database.py [LRG_GRCh37.bed] [output_filename]

The LRG database was downloaded from https://www.lrg-sequence.org/data/
"""
__author__ = "Kai"
__date__ = "2019/10/31"


import sys
import os
import re


if __name__ == "__main__":
    LRG = os.path.abspath(sys.argv[1])
    ofname = os.path.abspath(sys.argv[2])

    # allocate database into two dictionaries
    LRG_genes = {}
    LRG_transcripts = {}
    with open(LRG, "r") as fh:
        for line in fh:
            if line.startswith('track'):
                continue

            match = re.search(r"(LRG_\d+)\((\S+)\)", line.split()[3])
            if match:
                LRG_genes[match.group(2)] = match.group(1)
                continue

            match = re.search(r'(LRG_\d+)t1\((N[MR]_\S+)\)', line.split()[3])
            if match:
                LRG_transcripts[match.group(1)] = match.group(2).split("|")[0]
            
           
    # print some logging info 
    print('** The number of genes identified: {}'.format(len(LRG_genes)))
    print('** The number of transcripts identified: {}'.format(len(LRG_transcripts)))
    
    # write to output file
    with open(ofname, "w") as fh:
        for gene in LRG_genes:
            fh.write("{}\t{}\n".format(gene, LRG_transcripts[LRG_genes[gene]]))
    print("** Output has been written to {}".format(ofname))
