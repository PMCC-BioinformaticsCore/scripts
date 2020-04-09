#!/usr/bin/python3

"""
Title: addSymToDepthOfCoverage.py
Creator: Jason Li
Updated: Jiaan Yu
Date: 09-04-2020

Annotate gene to Gatk3 DepthOfCoverage interval_summary output
USAGE: 
    addSymToDepthOfCoverage.py -i <DepthOfCoverageOut> -b <annotatedBed> -o <annotatedDepthOfCoverageResult>

"""

import argparse

parser = argparse.ArgumentParser(description="Performance summary of bam")
parser.add_argument("-i", help="Gatk3 DepthOfCoverage interval_summary output", metavar = "INPUT", required=True)
parser.add_argument("-o", help="Output file name", metavar = "OUTPUT", required=True)
parser.add_argument("-bed", help="Annotated bed file", metavar = "BED", required=True)
args = parser.parse_args()


geneDict = {}
with open(args.bed, 'r') as bed_fh:
    for line in bed_fh:
        line_elems = line.rstrip("\n\r").split("\t")
        var_id = "%s:%d-%s" % (line_elems[0], int(line_elems[1])+1, line_elems[2])
        geneDict[var_id] = line_elems[3]

fh = open(args.i)
foh = open(args.o, 'w')
foh.write("GeneSym\t" + fh.readline())
for line in fh:
    line_elems = line.rstrip("\n\r").split("\t")
    var_id = line_elems[0]
    gene = geneDict[var_id]
    lineout = gene+"\t"+line
    foh.write(lineout)
fh.close()
foh.close()
