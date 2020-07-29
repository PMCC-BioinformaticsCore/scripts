#!/usr/bin/python3

"""
Title: calculate_strelka_somatic_ad_af
Author: Jiaan Yu
Date: 07-07-2020

This script
 - Extract and calculate AD and AF value for each variant (both SNVs and INDELs)
Based on https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#somatic


"""
import re
import argparse
import sys
from variant import Variant

AD_HEADER_LINE='##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Extracted allelic \
depths for the ref and alt alleles based on strelka recommendation \
- Calculated By Bioinformatics Dept">\n'
AF_HEADER_LINE='##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele \
Frequency, for each ALT allele, in the same order as listed - Calculated By \
Bioinformatics Dept">\n'

parser = argparse.ArgumentParser(description="Calculate AF for strelka somatic vcf")
parser.add_argument("-i", dest="input", help="Input strelka somatic vcf", required=True)
parser.add_argument("-o", dest="output", help="Output vcf", required=True)
args = parser.parse_args()

with open(args.input, "r") as f_vcf, open(args.output, "w") as f_vcf_out:
    for line in f_vcf:
        # Write the header lines
        if line.startswith("##"):
            f_vcf_out.write(line)
        elif line.startswith("#"):
            f_vcf_out.write(AD_HEADER_LINE)
            f_vcf_out.write(AF_HEADER_LINE)
            f_vcf_out.write(line)
            normal_index, tumor_index = (line.split().index("NORMAL"),
                                         line.split().index("TUMOR"))

        else:
            variant = Variant().process_somatic_variant(line, "strelka",
                                                        normal_index, tumor_index)
            # Raise expection for multi-allelic variants
            if "," in variant.ref or "," in variant.alt:
                print("The following varinat is a multi-allelic variant:\n" +
                      variant.write())
                sys.exit(0)

            else:
                for vtype in ['normal', 'tumor']:
                    if len(variant.ref) == len(variant.alt) == 1:
                        refCounts = int(variant.format[vtype][variant.ref+"U"].split(",")[0])
                        altCounts = int(variant.format[vtype][variant.alt+"U"].split(",")[0])
                    else:
                        refCounts = int(variant.format[vtype]["TAR"].split(",")[0])
                        altCounts = int(variant.format[vtype]["TIR"].split(",")[0])
                    variant.format[vtype]["AD"] = str(refCounts) + "," + str(altCounts)
                    if refCounts + altCounts != 0:
                        variant.format[vtype]["AF"] = str(round(altCounts / float(refCounts + altCounts), 2))
                    else:
                        variant.format[vtype]["AF"] = "0.00"
                f_vcf_out.write(variant.write(somatic=True))
