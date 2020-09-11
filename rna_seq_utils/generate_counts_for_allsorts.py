#!/usr/bin/python3

import sys
import argparse

"""
Title: generate_counts_for_allsorts.py
Creator: Jiaan Yu
Date: 11-09-2020

This script generate csv files from featureCounts or htseq-count for ALLSorts input
"""

recognised_modes = ["featureCounts", "htseq-count"]

parser = argparse.ArgumentParser(description="Generate csv files from featureCounts \
    or htseq-count for ALLSorts input")
optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
required.add_argument("-i", action='append', help="one of more input files", required=True)
required.add_argument("-o", help="output csv", required=True)
required.add_argument("--type", help="must be either featureCounts or htseq-count",
                      required=True, choices=recognised_modes)
required.add_argument("--samples", action='append', help="one of more labels (samples) for input file(s) in the same order", required=True)
parser._action_groups.append(optional)
required.add_argument("--gene_column", help="Specify which column to use to obtain gene_id, htseq-count only", required=False)
required.add_argument("--count_column", help="Specify which column to use to obtain count for genes, htseq-count only", required=False)
args = parser.parse_args()


combined_count_list=[]

if args.type == "featureCounts":
    for i, args.i in enumerate(args.i):
        gene_list=[""]
        count_list=[args.samples[i]]
        file=args.i[i]

        with open(file, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    pass
                elif line.startswith("Geneid"):
                    pass
                else:
                    line=line.strip("\n").split("\t")
                    gene=line[0]
                    count=line[6]
                    gene_list.append(gene)
                    count_list.append(count)
        combined_count_list=combined_count_list.append(count_list)


elif args.type == "htseq-count":
    gene_column=int(args.gene_column)-1
    count_column=int(args.count_column)-1
    for i, file in enumerate(args.i):
        gene_list=[""]
        count_list=[args.samples[i]]

        with open(file, 'r') as f:
            for line in f:
                if line.startswith("__"):
                    pass
                else:
                    line=line.strip("\n").split("\t")
                    gene=line[gene_column]
                    count=line[count_column]
                    gene_list.append(gene)
                    count_list.append(count)
    combined_count_list=combined_count_list.append(count_list)


else:
    sys.exit("Unrecognise type entered. Please enter featureCounts or htseq-count")


with open(args.o, 'w') as out_f:
    out_f.write(",".join(gene_list)+"\n")
    for i in combined_count_list:
        out_f.write(",".join(i)+"\n")
