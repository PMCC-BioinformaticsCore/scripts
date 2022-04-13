#!/usr/bin/python3

"""
Title: replace_N_fusion_vcf.py
Author: Jiaan Yu

For fusion vcf, replace the "N" base in ALT from REF

Example:
'replace_N_fusion_vcf.py'
  --input [File, vcf file from bcftools +fill-from-fasta]
  > [String, output vcf]

Example to run bcftools
bcftools +fill-from-fasta
  [File, vcf file from MegaFusion.py]
  -- 
  --column 'REF'
  --fasta [File, reference Fasta file]
  > 'output.fill.vcf'

"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, type=str, help="path to input vcf")
args = parser.parse_args()

with open(args.input, 'r') as f:
    for line in f:
        if line.startswith("#"):
            print(line.strip('\n'))
        else:
            line=line.strip('\n').split('\t')
            ref=line[3]
            line[4]=line[4].replace("N", ref)
            print("\t".join(line))

