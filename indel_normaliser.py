__doc__ = """
Normalise indel position

VCF specifications and HGVS system requires different representation of indel:
HGVS requires right-aligned indel, while VCF requires left-aligned.
therefore, sometimes we need to shift indels to the right or left until it is no
longer possible to do so.

**
So far only support vcf file with one sample 

Usage example:
[kai@admin]$ python3 indel_normaliser.py --left [vcf_file] [reference_fasta]
"""
__author__ = "Kai"
__date__ = "18/11/2019"


import argparse
import pandas as pd
import logging
import pysam
import os


logger = logging.getLogger()
logging.basicConfig(level=logging.INFO, format="%(asctime)s -- %(message)s")


def normaliser(vcf_file, reference, direction, ofname):
    """
    * vcf uses 0-based, inclusive system
    """
    logger.info("shifting indel...")
    with pysam.VariantFile(vcf_file, "r") as vcfin, pysam.VariantFile(ofname, "w", header=vcfin.header) as vcfout:        
        for record in vcfin:
            # filter out SNVs
            if (record.rlen == 1) and (len(record.alts[0]) == 1):
                vcfout.write(record)
                continue
            
            # deletion
            if record.rlen != 1:
                chrom = record.chrom
                start = record.start
                stop = record.stop
                while True:
                    if direction == 'left':
                        start -= record.rlen
                        stop -= record.rlen
                        genomic_coordinate = (chrom, start, stop)
                    elif direction == 'right':
                        start += record.rlen
                        stop += record.rlen
                        genomic_coordinate = (chrom, start, stop)

                    shifted_ref = extract_seq(genomic_coordinate, reference).upper()
                    if shifted_ref == record.ref:
                        record.start, record.stop = start, stop
                    else:
                        vcfout.write(record)
                        break
            
            # insertion
            elif len(record.alts[0]) != 1:
                chrom = record.chrom
                start = record.start
                stop = record.stop
                while True:
                    if direction == 'left':
                        start -= record.rlen
                        stop -= record.rlen
                        genomic_coordinate = (chrom, start, stop)
                    elif direction == 'right':
                        start += record.rlen
                        stop += record.rlen
                        genomic_coordinate = (chrom, start, stop)

                    shifted_ref = extract_seq(genomic_coordinate, reference).upper()
                    if shifted_ref == record.ref:
                        record.start, record.stop = start, stop
                    else:
                        vcfout.write(record)
                        break
    logger.info("The output file was written to %s", ofname)


def extract_seq(genomic_coordinate, reference_genome):
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
                        seq = line[start_index: end_index]
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
    if not os.path.exists(fadix):
        raise FileNotFoundError("You need to index your fasta file using samtools!")
    df = pd.read_csv(fadix, sep="\t", header=None)
    for index, row in df.iterrows():
        fadix_dict[row[0]] = (row[1], row[2], row[3])
    return fadix_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left", action='store_true', help="shift indel to the left")
    parser.add_argument("--right", action='store_true', help="shift indel to the right")
    parser.add_argument("vcf", help="the input vcf file")
    parser.add_argument("reference", help="the reference fasta file")
    args = parser.parse_args()

    vcfin = os.path.abspath(args.vcf)
    if args.left:
        ofname = os.path.join(os.path.dirname(vcfin), 
                              os.path.basename(vcfin).split(".")[0] + ".left_shifted.vcf")
        normaliser(vcfin, args.reference, 'left', ofname)
    elif args.right:
        ofname = os.path.join(os.path.dirname(vcfin), 
                              os.path.basename(vcfin).split(".")[0] + ".right_shifted.vcf")
        normaliser(vcfin, args.reference, 'right', ofname)
    