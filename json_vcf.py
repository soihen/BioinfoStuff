""" convert VCF to JSON or JSON to VCF """

import pysam
import json
import copy
import argparse



VARIANT_DICT = {"chrom": None, "start": None, "end": None, "ref": None, "alt": None,
                "ID": None, "FILTER": None, "QUAL": None, "VCF_FORMAT": {}, "VCF_INFO": {}}


def vcf2json(fname, ofname):
    """
    VCF -> JSON
    Noted that the generate JSON file contains VCF header as well,
    which is required to convert back to VCF.
    """
    with pysam.VariantFile(fname, "r") as vcf:
        header = vcf.header
        sampleID = header.samples[0]
        vcf_dict = {"sampleID": sampleID, "vcf_header": str(header), "variants": []}
        for rec in vcf:
            variant_dict = copy.deepcopy(VARIANT_DICT)
            variant_dict['chrom'] = rec.chrom
            variant_dict['start'] = rec.start
            variant_dict['end'] = rec.stop
            variant_dict['ref'] = rec.ref
            variant_dict['alt'] = rec.alts[0]
            variant_dict['ID'] = rec.id if rec.id else "."
            variant_dict['FILTER'] = rec.filter.keys()
            variant_dict['QUAL'] = rec.qual if rec.qual else "."

            # Below two lines are used to round all digits 
            # that the dicts contained to three decimal places
            format_dict = {k: list(map(lambda x: round(x, 3) \
                if isinstance(x, float) else x, v)) \
                for k, v in dict(rec.samples[sampleID]).items()}

            info_dict = {k: list(map(lambda x: round(x, 3) \
                if isinstance(x, float) else x, v)) \
                if isinstance(v ,tuple) else v \
                for k, v in dict(rec.info).items()}
            
            variant_dict['VCF_FORMAT'] = format_dict
            variant_dict['VCF_INFO'] = info_dict
            vcf_dict['variants'].append(variant_dict)
    
    with open(ofname, "w") as fw:
        json.dump(vcf_dict, fw)
 

def json2vcf(fname, ofname):
    with open(fname, "r") as fh:
        vcf_dict = json.load(fh)
    
    header = vcf_dict['vcf_header']
    sampleID = vcf_dict['sampleID']
    variants = vcf_dict['variants']

    with open(ofname, "w") as fw:
        fw.write(header)

        for var in variants:
            # be aware of the type and the number of INFO/FORMAT
            info_list = []
            for key, value in var['VCF_INFO'].items():
                if not value:
                    info_list.append(key)
                elif isinstance(value, list) and len(value) == 1:
                    info_list.append(key+"="+str(value[0]))
                elif isinstance(value, list) and len(value) > 1:
                    info_list.append(key + "=" + ",".join(str(value)))
                else:
                    info_list.append(key + "=" + str(value))
            info_str = ";".join(info_list)

            format_list = []
            for key, value in var['VCF_FORMAT'].items():
                if key == 'GT':
                    format_list.append("/".join(map(str, value)))
                elif isinstance(value, list):
                    format_list.append(",".join(map(str, value)))
            format_str = ":".join(format_list)

            rec = [var['chrom'], str(var['start']), var['ID'], 
                    var['ref'], var['alt'], str(var['QUAL']), 
                    ",".join(var['FILTER']), info_str,
                    ":".join(var['VCF_FORMAT'].keys()), 
                    format_str]
            
            fw.write("\t".join(rec)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input filename, could be either .json or .vcf file")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    if args.input.endswith("vcf"):
        vcf2json(args.input, args.output)
    elif args.input.endswith("json"):
        json2vcf(args.input, args.output)
