#!/usr/bin/python3

# For fusion vcf, replace the "N" base in ALT from REF

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

