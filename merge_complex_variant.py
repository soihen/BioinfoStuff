__doc__ = """
merge variants belonging to a same phasing group

GATK can give a PID (physical phasing ID) to each unique variant.
Variants with the same PID can be merged into one complex variant

Usage:
[kai@admin]$ python3 merge_complex_variant.py [input_vcf] [reference]

***
Currently only support for vcf file with one sample only
"""
__author__ = "Kai"
__version__ = "v0.0"
__date__ = "21/11/2019"


import pysam
from pyfaidx import Fasta
import os
from collections import defaultdict
from collections import Counter
import argparse


def main(vcf_file, ofname, reference):
    pids = identify_phasing_groups(vcf_file)
    reference = Fasta(reference)
    counter = Counter()
    with pysam.VariantFile(vcf_file, "r") as vcfin, pysam.VariantFile(ofname, "w", header=vcfin.header) as vcfout:
        sampleID = vcfin.header.samples[0]
        for record in vcfin:
            if 'PID' not in record.samples[sampleID]:
                vcfout.write(record)
                continue
            # make sure only one record will be written to outfile
            counter[record.samples[sampleID]['PID']] += 1
            if counter[record.samples[sampleID]['PID']] == len(pids[record.samples[sampleID]['PID']]):
                merged_record = merge_variants(pids[record.samples[sampleID]['PID']], reference)
                vcfout.write(merged_record)
            

def merge_variants(variants_list, reference):
    """
    merge all variants in the list into a complex variant
    """
    merged_record = variants_list[0]
    for variant_record in variants_list[1:]:
        ref_seq = reference[variant_record.chrom][merged_record.start: variant_record.stop].seq
        gap_ref = reference[variant_record.chrom][merged_record.stop: variant_record.start].seq
        alt_seq = merged_record.alts[0] + gap_ref + variant_record.alts[0]
        merged_record.stop = variant_record.stop
        merged_record.ref = ref_seq
        merged_record.alts = tuple([alt_seq])
    return merged_record


def identify_phasing_groups(vcf_file):
    """
    identify variants on the same phasing groups
    """
    pids = defaultdict(list)
    with pysam.VariantFile(vcf_file, "r") as vcf:
        sampleID = vcf.header.samples[0]
        for record in vcf:
            if 'PID' in record.samples[sampleID]:
                pids[record.samples[sampleID]['PID']].append(record)
    return pids



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", help="the input vcf file")
    parser.add_argument("reference", help="the reference fasta file")
    parser.add_argument("-o", "--output", help="the output filename")
    args = parser.parse_args()

    if args.output:
        main(args.vcf, args.output, args.reference)
    else:
        infame = os.path.abspath(args.vcf)
        ofname = os.path.join(os.path.dirname(infame), os.path.basename(infame).split(".")[0]+".phasing_merged.vcf")
        main(infame, ofname, args.reference)