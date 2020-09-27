"""
Convert between VCF format to TXT format (e.g. Annovar output txt file);
VCF follows 'leading base rule', which means one leading base should be included in delins variants; whereas TXT does not
Both file format use 1-base coordinate system, TXT file format contains a 'Stop' column in particular
"""

import argparse
import pysam
import sys
from pyfaidx import Fasta


VCF_HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chrM,length=16571,assembly=hg19>
##contig=<ID=chr1,length=249250621,assembly=hg19>
##contig=<ID=chr2,length=243199373,assembly=hg19>
##contig=<ID=chr3,length=198022430,assembly=hg19>
##contig=<ID=chr4,length=191154276,assembly=hg19>
##contig=<ID=chr5,length=180915260,assembly=hg19>
##contig=<ID=chr6,length=171115067,assembly=hg19>
##contig=<ID=chr7,length=159138663,assembly=hg19>
##contig=<ID=chr8,length=146364022,assembly=hg19>
##contig=<ID=chr9,length=141213431,assembly=hg19>
##contig=<ID=chr10,length=135534747,assembly=hg19>
##contig=<ID=chr11,length=135006516,assembly=hg19>
##contig=<ID=chr12,length=133851895,assembly=hg19>
##contig=<ID=chr13,length=115169878,assembly=hg19>
##contig=<ID=chr14,length=107349540,assembly=hg19>
##contig=<ID=chr15,length=102531392,assembly=hg19>
##contig=<ID=chr16,length=90354753,assembly=hg19>
##contig=<ID=chr17,length=81195210,assembly=hg19>
##contig=<ID=chr18,length=78077248,assembly=hg19>
##contig=<ID=chr19,length=59128983,assembly=hg19>
##contig=<ID=chr20,length=63025520,assembly=hg19>
##contig=<ID=chr21,length=48129895,assembly=hg19>
##contig=<ID=chr22,length=51304566,assembly=hg19>
##contig=<ID=chrX,length=155270560,assembly=hg19>
##contig=<ID=chrY,length=59373566,assembly=hg19>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

TXT_HEADER = "Chrom\tStart\tEnd\tRef\tAlt\n"
TSV_HEADER = "Chrom\tPOS\tRef\tAlt\n"
CSV_HEADER = "Chrom,POS,Ref,Alt\n"


def rm_leading_base(chrom, pos, ref, alt):
    # does not contain a leading base
    # MUST be a substitution
    if ref[0] != alt[0]:
        return chrom, str(pos), str(pos+len(ref)-1), ref, alt
    
    # contain a leading base
    # deletion
    if len(ref) > len(alt):
        alt = '-' if len(alt) == 1 else alt[1:]
        return chrom, str(pos+1), str(pos+len(ref)-2), ref[1:], alt
    
    # insertion
    ref = '-' if len(ref) == 1 else ref[1:]
    return chrom, str(pos), str(pos), ref, alt[1:]


def add_leading_base(chrom, pos, ref, alt, fasta):
    # substitution
    if len(ref) == len(alt) and (ref != '-' and alt != '-'):
        return chrom, str(pos), ref, alt
    
    # delins or MNV
    # deletion
    if alt == '-':
        leading_base = fasta[chrom][pos-2: pos-1].seq.upper()
        return chrom, str(pos-1), leading_base+ref, leading_base
    
    # MNV deletion
    if len(ref) > len(alt):
        leading_base = fasta[chrom][pos-2: pos-1].seq.upper()
        return chrom, str(pos-1), leading_base+ref, leading_base+alt
    
    # insertion
    if ref == '-':
        leading_base = fasta[chrom][pos-1: pos].seq.upper()
        return chrom, str(pos), leading_base, leading_base+alt
    
    # MNV insertion
    if len(ref) < len(alt):
        leading_base = fasta[chrom][pos-1: pos].seq.upper()
        return chrom, str(pos), leading_base+ref, leading_base+alt


def vcf2txt(fname, ofname, fmat, index_chr, index_pos, index_ref, index_alt):
    """ rm one leading base from VCF or VCF-like format file """
    contents = []

    if fmat == 'vcf':
        with pysam.VariantFile(fname, "r") as vcfin:
            for rec in vcfin:
                contents.append(rm_leading_base(rec.chrom, rec.pos, rec.ref, rec.alts[0]))
    else:
        if fmat == 'csv':
            sep = ","
        elif fmat == 'tsv':
            sep = "\t"
        with open(fname, "r") as fh:
            for index, line in enumerate(fh):
                if index == 0:
                    continue
                line = line.strip().split(sep)
                # skip empty line (sometimes file contain an empty line in the end  )
                if not line[0]: 
                    continue
                chrom = line[index_chr]
                pos = int(line[index_pos])
                ref = line[index_ref]
                alt = line[index_alt]
                contents.append(rm_leading_base(chrom, pos, ref, alt))
    
    with open(ofname, "w") as fw:
        fw.write(TXT_HEADER)
        for row in contents:
            fw.write("\t".join(row)+"\n")


def txt2vcf(fname, ofname, fmat, index_chr, index_start, index_ref, index_alt, reference):
    """ add one leading base to txt file """
    contents = []
    # beware Fasta object is a zero-base and half-closed system
    fasta = Fasta(reference)

    with open(fname, "r") as fh:
        for index, line in enumerate(fh):
            if index == 0:
                continue
            line = line.strip().split("\t")
            chrom = line[index_chr]
            pos = int(line[index_start])
            ref = line[index_ref]
            alt = line[index_alt]
            contents.append(add_leading_base(chrom, pos, ref, alt, fasta))
    
    # write output
    if fmat == 'vcf':
        with open(ofname, "w") as fw:
            fw.write(VCF_HEADER)
            for row in contents:
                outrow = [row[0], row[1], ".", row[2], row[3], ".", "PASS", "."]
                fw.write("\t".join(outrow)+"\n")
    elif fmat == 'tsv':
        with open(ofname, "w") as fw:
            fw.write(TSV_HEADER)
            for row in contents:
                fw.write("\t".join(row)+"\n")
    elif fmat == 'csv':
        with open(ofname, "w") as fw:
            fw.write(CSV_HEADER)
            for row in contents:
                fw.write(",".join(row)+"\n")


def main(args):
    if args.cmd == 'vcf2txt':
        vcf2txt(args.input, args.output, args.format, 
                args.chr, args.pos, args.ref, args.alt)
    elif args.cmd == 'txt2vcf':
        txt2vcf(args.input, args.output, args.format,
                args.chr, args.start, args.ref, args.alt, 
                args.reference)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)

    subparsers = parser.add_subparsers(dest='cmd')
    vcf2txt_parser = subparsers.add_parser('vcf2txt', help="convert from vcf (1-based coordinate system) to txt (0-based coordinate system)")
    vcf2txt_parser.add_argument("-i", "--input", required=True, help="Input filename")
    vcf2txt_parser.add_argument("-o", "--output", required=True, help="Output filename")
    vcf2txt_parser.add_argument("-f", "--format", default="vcf", choices=['vcf', 'csv', 'tsv'], help="input file format")
    vcf2txt_parser.add_argument("--chr", type=int, help="index for chromosome column")
    vcf2txt_parser.add_argument("--pos", type=int, help="index for position column")
    vcf2txt_parser.add_argument("--ref", type=int, help="index for reference seq column")
    vcf2txt_parser.add_argument("--alt", type=int, help="index for variant seq column")

    txt2vcf_parser = subparsers.add_parser('txt2vcf', help="convert from txt (0-based coordinate system) to vcf (1-based coordinate system)")
    txt2vcf_parser.add_argument("-i", "--input", required=True, help="Input filename")
    txt2vcf_parser.add_argument("-o", "--output", required=True, help="Output filename")
    txt2vcf_parser.add_argument("-f", "--format", default="vcf", choices=['vcf', 'tsv', 'csv'], help="output file format")
    txt2vcf_parser.add_argument("--chr", type=int, default=0, help="index for chromosome column")
    txt2vcf_parser.add_argument("--start", type=int, default=1, help="index for start position column")
    txt2vcf_parser.add_argument("--ref", type=int, default=3, help="index for reference seq column")
    txt2vcf_parser.add_argument("--alt", type=int, default=4, help="index for variant seq column")
    txt2vcf_parser.add_argument("-r", "--reference", required=True, help="path to the reference Fasta file")

    args = parser.parse_args()

    main(args)
