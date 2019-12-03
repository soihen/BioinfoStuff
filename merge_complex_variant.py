__doc__ = """
merge variants belonging to a same phasing group

GATK can give a PID (physical phasing ID) to each unique variant.
Variants with the same PID can be merged into one complex variant

***
multiallelic sites in the input VCF file should be splitted in advance:
`bcftools norm -m both [input] -o [output]`
The VCF file should also be sorted in advance

Usage:
[kai@admin]$ python3 merge_complex_variant.py [input_vcf] [reference]

***
Currently only support for vcf file with one sample only

"""
__author__ = "Kai"
__version__ = "v0.2"
__date__ = "22/11/2019"


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
                if merged_record:
                    print("* {}:{}-{} {}->{}".format(merged_record.chrom, merged_record.start, merged_record.stop,
                                                     merged_record.ref, merged_record.alts[0]))
                    vcfout.write(merged_record)
    print("--------------")
    print("* The output file was written to {}".format(ofname))


def merge_variants(variants_list, reference):
    """
    merge all variants in the list into a complex variant    
    """
    merged_record = variants_list[0]
    for variant_record in variants_list[1:]:
        # two variants lists in order, and not overlap with each other
        if merged_record.stop < variant_record.start:
            ref_seq = reference[variant_record.chrom][merged_record.start: variant_record.stop].seq.upper()
            gap_ref = reference[variant_record.chrom][merged_record.stop: variant_record.start].seq.upper()
            alt_seq = merged_record.alts[0] + gap_ref + variant_record.alts[0]
            merged_record.stop = variant_record.stop
            merged_record.ref = ref_seq
            merged_record.alts = tuple([alt_seq])
            
        # two variants next to each other
        elif merged_record.stop == variant_record.start:
            ref_seq = reference[variant_record.chrom][merged_record.start: variant_record.stop].seq.upper()
            alt_seq = merged_record.alts[0] + variant_record.alts[0]
            merged_record.stop = variant_record.stop 
            merged_record.ref = ref_seq
            merged_record.alts = tuple([alt_seq])

        # two variants overlapped
        elif merged_record.start < variant_record.start < merged_record.stop < variant_record.stop:
            # if two alt seqs are same, then merge, otherwise, discard them both
            var1_alt = merged_record.alts[0][variant_record.start-merged_record.start: variant_record.stop-merged_record.start]
            var2_alt = variant_record.alts[0][:merged_record.stop-variant_record.start]
            if var1_alt != var2_alt:
                print("** Warning: those two variants mutated to different alleles, abondon them both")
                print("* {}:{}-{} {}->{}".format(merged_record.chrom, merged_record.start, merged_record.stop,
                                                 merged_record.ref, merged_record.alts[0]))
                print("* {}:{}-{} {}->{}".format(variant_record.chrom, variant_record.start, variant_record.stop,
                                                 variant_record.ref, variant_record.alts[0]))
                return None
            merged_record.stop = variant_record.stop
            ref_seq = reference[variant_record.chrom][merged_record.start: variant_record.stop].seq.upper()
            merged_record.ref = ref_seq
            alt_seq = merged_record.alts[0] + variant_record.alts[0][merged_record.stop-variant_record.start:]

        # the first variant encomprise the second variant
        elif merged_record.start < variant_record.start < variant_record.stop <= merged_record.stop:
            ref_seq = reference[variant_record.chrom][merged_record.start: variant_record.stop].seq.upper()
            # if two alt seqs are same, then merge, otherwise, discard them both
            var1_alt = merged_record.alts[0][variant_record.start-merged_record.start: variant_record.stop-merged_record.start]
            var2_alt = variant_record.alts[0]
            if var1_alt != var2_alt:
                print("** Warning: those two variants mutated to different alleles, abondon them both")
                print("*** {}:{}-{} {}->{}".format(merged_record.chrom, merged_record.start, merged_record.stop,
                                                   merged_record.ref, merged_record.alts[0]))
                print("*** {}:{}-{} {}->{}".format(variant_record.chrom, variant_record.start, variant_record.stop,
                                                   variant_record.ref, variant_record.alts[0]))
                return None
        
        else:
            raise Exception("Unexpected variant position, sort your vcf before use this scirpt")
   
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
    print("* {} phasing group(s) have been identified".format(len(pids)))  
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
