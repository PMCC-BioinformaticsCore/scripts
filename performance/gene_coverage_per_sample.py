#!/usr/bin/python3

import sys
import threading
import os
import re
from collections import defaultdict as dd
import csv
import argparse


"""
Title: gene_coverage_per_sample
Creator: Jason Li
Date: 27-02-2020
Modified: 23-06-2020 by Jiaan Yu

This script groups coverage by gene per sample,
it uses coverage stats (bedtools) from 
each region specified in the bed file, and 
averages(weighted) the coverage by gene 

Pre-requisite: run bedtools coverageBed -abam $bam -b $bed -hist

Example Bed File:
16	23614779	23614990	PALB2
16	23619184	23619333	PALB2
16	23625324	23625412	PALB2;TP53

Example bedtools output:

19	852271	852451	ELANE	9	6	180	0.0333333
19	852271	852451	ELANE	10	32	180	0.1777778
19	852271	852451	ELANE	11	26	180	0.1444445
19	852271	852451	ELANE	12	29	180	0.1611111
19	852271	852451	ELANE	13	13	180	0.0722222
19	852271	852451	ELANE	14	21	180	0.1166667


"""
lock = threading.Lock()


class Thread(threading.Thread):
    def __init__(self, t, *args):
        threading.Thread.__init__(self, target=t, args=args)
        self.start()


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def add_empty_row(matrix, numCols):
    matrix.append([0] * numCols)


def add_empty_row_of_lists(matrix, numCols, sizeList):
    matrix.append([[0.0] * sizeList] * numCols)


def proc(procId, sample_list, geneObject, regionsObject, folds):
    for sample in sample_list:

        # Initialise values
        region_depth, region_vals = dd(float), dd(list)
        region_len, region_folds, region_list = {}, {}, []

        gene_depth, gene_vals = dd(float), dd(list)
        gene_len, gene_folds, gene_list = dd(float), {}, []

        # folds = ["1", "10", "20", "100", "150", "200", "500", "1000"]
        folds = [float(i) for i in folds]

        print("Thread:", procId)
        print("Sample: ", sample[0])
        # print("Bedtools...")
        # print("Generating stats")
        dirt = os.getcwd()

        with open(sample[1], "r") as bed:
            for line in bed:
                if not line.startswith("all"):
                    line = line.strip().split("\t")
                    genes = line[3].split(";")
                    for gene_coln in genes:
                        gene = re.sub("_ex.*", "", gene_coln)
                        reg = line[3] + ":" + line[0] + ":" + line[1] + "-" + line[2]
                        depth, count, length, percent = (
                            float(line[4]),
                            int(line[5]),
                            float(line[6]),
                            float(line[7]),
                        )

                        if reg not in region_list:
                            region_len[reg] = length
                            region_folds[reg] = [
                                percent if depth >= f else 0 for f in folds
                            ]
                            region_list.append(reg)
                            gene_len[gene] += length
                        else:
                            region_folds_list = [
                                percent if depth >= f else 0 for f in folds
                            ]
                            region_folds[reg] = [
                                i + j for i, j in zip(region_folds_list, region_folds[reg])
                            ]
                        region_depth[reg] += depth * count
                        region_vals[reg].append(depth)

                        # print(region_depth, region_vals)
                        # print(region_len, region_folds, region_list)

                        if gene not in gene_list:
                            gene_folds[gene] = [
                                percent * length if depth >= f else 0 for f in folds
                            ]
                            gene_list.append(gene)
                        else:
                            gene_folds_list = [
                                percent * length if depth >= f else 0 for f in folds
                            ]
                            gene_folds[gene] = [
                                i + j for i, j in zip(gene_folds_list, gene_folds[gene])
                            ]
                        gene_depth[gene] += depth * count
                        gene_vals[gene].append(depth)

        # print(gene_depth, gene_vals)
        # print(gene_len, [f / gene_len[gene] for f in gene_folds[gene]], gene_list)

        # write output
        for gene in gene_list:
            row = [
                sample[0],
                gene,
                "{:0.2f}".format(gene_depth[gene] / gene_len[gene]),
                str(min(gene_vals[gene])),
                str(max(gene_vals[gene])),
            ] + ["{:0.2f}".format(f / gene_len[gene] * 100) for f in gene_folds[gene]]
            # print(",".join(row))
            with lock:
                geneObject.writerow(row)

        for reg in region_list:
            row = [
                sample[0],
                reg,
                "{:0.2f}".format(region_depth[reg] / region_len[reg]),
                str(min(region_vals[reg])),
                str(max(region_vals[reg])),
            ] + ["{:0.2f}".format(f * 100) for f in region_folds[reg]]
            # print(",".join(row))
            with lock:
                regionsObject.writerow(row)


def chunker(seq, size):
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


def main():
    procId = os.getpid()

# Rewrite the options using argparse
# Two options are deprecated: -b and -d
    parser = argparse.ArgumentParser(description="Gene or region coverage of bam")
    parser.add_argument(
        "-l", "--list", dest="list",
        help="List file: A tsv file contains SampleName\tPathToBedtoolsOutput on each line",
        required=False)
    parser.add_argument(
        "-n", "--name", dest="name", help="Sample name if list not used", required=False)
    parser.add_argument(
        "-p", "--path", dest="path", help="Path to bedtools output if list not used",
        required=False)
    parser.add_argument(
        "-b", "--bed", dest="bed", help="(Deprecated option) Bed file", required=False)
    parser.add_argument(
        "-g", "--gene", dest="gene", help="Output gene file", required=False)
    parser.add_argument(
        "-r", "--region", dest="region", help="Output region file", required=False)
    parser.add_argument(
        "-f", "--folds", dest="folds",
        help="Folds, quoted and commna sepparated, default 1,10,20,100",
        required=False)
    parser.add_argument(
        "-d", "--remove_duplicates", action="store_false", default=True, dest="dups",
        help="(Deprecated option) Remove marked duplicates in analysis, default:false",
        required=False)
    parser.add_argument(
        "-t", "--threads", dest="threads", default=32,
        type=int, help="number of threads, default:32", required=False)

    args = parser.parse_args()
    if args.gene is None or args.region is None:
        parser.print_help()
        sys.exit("ERROR: Missing argument")
    if args.folds is None:
        folds = ["1", "10", "20", "100", "150", "200", "500", "1000"]
    else:
        folds = args.folds.split(",")

    if args.list is None and (args.name is None and args.path is None):
        parser.print_help()
        sys.exit("ERROR: Missing argument, plese specify list file or sample name and path")

    if args.list is not None:
        num_samples = file_len(args.list)
        indexFile = open(args.list, "r")
        indexObject = csv.reader(indexFile, delimiter="\t")

    else:
        indexObject = [[args.name, args.path]]
        num_samples = 1
    geneFile = open(args.gene, "w+")
    geneObject = csv.writer(geneFile, delimiter="\t", lineterminator="\n")

    regionsFile = open(args.region, "w+")
    regionsObject = csv.writer(regionsFile, delimiter="\t", lineterminator="\n")

    # Write header
    geneObject.writerow(
        ["Sample", "Gene", "Mean coverage", "Min Coverage", "Max Coverage"]
        + [x + "x" for x in folds]
    )
    regionsObject.writerow(
        ["Sample", "Region", "Mean coverage", "Min Coverage", "Max Coverage"]
        + [x + "x" for x in folds]
    )
    # parse bed file for each sample
    count = 0
    div = num_samples / args.threads
    mod = num_samples % args.threads
    if num_samples <= args.threads:
        for sample in indexObject:
            Thread(
                proc,
                str(count) + "_" + str(procId),
                [sample],
                geneObject,
                regionsObject,
                folds,
            )
            count += 1
    else:
        size = div if mod is 0 else div + 1
        seq = []
        for sample in indexObject:
            seq.append(sample)
        for sample in chunker(seq, size):
            Thread(
                proc,
                str(count) + "_" + str(procId),
                sample,
                geneObject,
                regionsObject,
                folds,
            )
            count += 1


if __name__ == "__main__":
    main()
