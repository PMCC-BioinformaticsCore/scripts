#!/usr/bin/python3

"""
Title: calculate_strelka_germline_af
Author: Jiaan Yu
Date: 03-07-2020

This script
 - Calculate AF value for each variant (both SNVs and INDELs)
 - Make DP for INDELs variant by coping over the value of DPI

"""
import re
import argparse
import sys
from variant import Variant

AF_HEADER_LINE='####INFO=<ID=AF,Number=A,Type=Float,Description="Allele \
Frequency, for each ALT allele, in the same order as listed - Calculated By \
Bioinformatics Dept">>\n'

parser = argparse.ArgumentParser(description="Calculate AF for strelka germline vcf")
parser.add_argument("-i", dest="input", help="Input strelka germline vcf", required=True)
parser.add_argument("-o", dest="output", help="Output vcf", required=True)
args = parser.parse_args()

with open(args.input, "r") as f_vcf, open(args.output, "w") as f_vcf_out:
    for line in f_vcf:
        # Write the header lines
        if line.startswith("##"):
            f_vcf_out.write(line)
        elif line.startswith("#"):
            f_vcf_out.write(AF_HEADER_LINE)
            f_vcf_out.write(line)

        # Calcualte AF for strelka variants
        else:
            variant = Variant().read_variant(line)
            if "DP" in variant.format.keys():
                if variant.format["DP"] == "0":
                    variant.format["AF"] = "."
                else:
                    variant.format["AF"] = ",".join([
                        str(round(float(i) / float(variant.format["DP"]), 2))
                        for i in variant.format["AD"].split(",")[1:]])
            else:
                variant.format["DP"] = variant.format["DPI"]
                if variant.format["DPI"] == "0":
                    variant.format["AF"] = "."
                else:
                    variant.format["AF"] = ",".join([
                        str(round(float(i) / float(variant.format["DPI"]), 2))
                        for i in variant.format["AD"].split(",")[1:]])
            f_vcf_out.write(variant.write())
