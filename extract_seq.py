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
                        seq = line[start_index: end_index+1]
                        return seq

    elif end_line_index > start_line_index:
        seq = ''
        with open(reference_genome, "r") as fh:
            # set offset 
            fh.seek(fadix_dict[genomic_coordinate[0]][1])
            for index, line in enumerate(fh):
                    if index == start_line_index:
                        seq += line[start_index: fadix_dict[genomic_coordinate[0]][2]]
                    elif start_line_index < index < end_line_index:
                        seq += line[:fadix_dict[genomic_coordinate[0]][2]]
                    elif index == end_line_index:
                        seq += line[:end_index+1]
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
