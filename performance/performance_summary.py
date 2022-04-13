#!/usr/bin/python3

"""
Title: performance_summary
Author: Jiaan Yu
Date: 18-02-2020

Summary of bam

Examples:
For whole genome data
performance_summary.py
  --flagstat [samtools flagstat output on bam]
  --rmdup_flagstat [samtools flagstat on rmdup bam]
  --collect_insert_metrics [GATK CollectInsertSizeMetrics output text file]
  --coverage [bedtools genomeCoverageBed on bam]
  -o [output filename]

For targeted data
performance_summary.py
  --flagstat [samtools flagstat output on bam]
  --rmdup_flagstat [samtools flagstat on rmdup bam]
  --target_flagstat [samtools flagstat on rmdup + targeted bed bam]
  --collect_insert_metrics [GATK CollectInsertSizeMetrics output text file]
  --coverage [bedtools genomeCoverageBed on bam]
  -o [output filename]

"""

import os
import sys
import argparse
import errno

parser = argparse.ArgumentParser(description="Performance summary of bam")
optional = parser._action_groups.pop()
# required arguments
required = parser.add_argument_group("required arguments")
required.add_argument("--flagstat",  help="output of samtools flagstat on bam", required=True)
required.add_argument("--collect_insert_metrics",  help="output of CollectInsertMetrics (GATK or Picard) on bam", required=True)
required.add_argument("--coverage",  help="output of bedtools coverageBed for targeted bam; bedtools genomeCoverageBed for whole genome bam", required=True)
required.add_argument("-o", help="output summary csv name", required=True)
# optional arguments
parser._action_groups.append(optional)
parser.add_argument("--target_flagstat",  help="output of samtools flagstat of bam target on target bed. Only specified for targeted bam")
parser.add_argument("--rmdup_flagstat",  help="output of samtools flagstat of removed duplicates bam. File to be used to extract mapping infomation if specified, instead of the --flagstat file.")
parser.add_argument('--genome', action='store_true', help="calculate statistics for whole genome data. --target_flagstat must not be speicified")
args = parser.parse_args()

# Check if output directory exist
# Taken from https://stackoverflow.com/a/600612/119527
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def safe_open_w(path):
    ''' Open "path" for writing, creating any parent directories as needed.
    '''
    mkdir_p(os.path.dirname(path))
    return open(path, 'w')


# check input
if args.genome and args.target_flagstat:
    sys.exit("Both --genome and --target_flagstat are specified. Please check the input")
elif (not args.genome) and (not args.target_flagstat):
    sys.exit("One of --genome and --target_flagstat must be specified")
for f in [args.flagstat, args.collect_insert_metrics, args.coverage, args.target_flagstat, args.rmdup_flagstat]:
    if f and (not os.path.exists(f)):
        sys.exit("File %s could not be found on disk" % (f))


# Store the output into dictionary 
header=['Total reads','Mapped reads','% Reads mapped','% Mapped reads duplicates','Total reads minus duplicates','Mapped reads minus duplicates','Paired in sequencing','Read1','Read2','Properly paired', 'With itself and mate mapped','Singletons','With mate mapped to a different chr','With mate mapped to a different chr (mapQ>=5)','% Reads OnTarget','% OnTarget duplicates','OnTarget paired in sequencing','OnTarget read1','OnTarget read2','OnTarget properly paired', 'OnTarget With itself and mate mapped','OnTarget singletons','OnTarget with mate mapped to a different chr','OnTarget with mate mapped to a different chr (mapQ>=5)','% Target bases >=1-fold Coverage','% Target bases >=10-fold Coverage','% Target bases >=20-fold Coverage','% Target bases >=100-fold Coverage','Mean coverage for target bases','% Target bases with one fifth of mean coverage', 'Median Fragment Length']
output_dict = {}


# Get median insert size from GATK CollectInsertSizeMetrics
with open(args.collect_insert_metrics, "r") as f:
    for (i, line) in enumerate(f):
        if i == 7:
            output_dict['Median Fragment Length'] = line.split("\t")[0]
            #print(median_insert_size)
            break


# flagstat of bam
flagstat_total_names = ['Total reads', 'Secondary reads', 'Supplementary reads', 'Duplicate reads', 'Mapped reads']
flagstat_total_names2 = ['Paired in sequencing', 'Read1','Read2', 'Properly paired', 'With itself and mate mapped', 'Singletons', 'With mate mapped to a different chr', 'With mate mapped to a different chr (mapQ>=5)']

# Total / Duplicates / Mapped reads get from flagstat
flagstat_lines = open(args.flagstat, 'r').readlines()
for line, name in zip(flagstat_lines, flagstat_total_names):
    output_dict[name] = line.split(" ")[0]

if not args.rmdup_flagstat: # Use the flagstat with duplicates
    #print(flagstat_lines)
    for i, line in enumerate(flagstat_lines[5:]):
        output_dict[flagstat_total_names2[i]] = line.split(" ")[0]
    output_dict['Mapped reads minus duplicates'] = "NA"

else:
    rmdup_flagstat_lines = open(args.rmdup_flagstat, 'r').readlines()
    for i, line in enumerate(rmdup_flagstat_lines[5:]):
        output_dict[flagstat_total_names2[i]] = line.split(" ")[0]
    output_dict['Mapped reads minus duplicates'] = str(int(output_dict['Mapped reads'])\
        - int(output_dict['Duplicate reads']))

output_dict['Total reads minus duplicates'] = str(
    int(output_dict['Total reads']) - int(output_dict['Duplicate reads']))

if output_dict['Total reads'] != '0':
    output_dict['% Reads mapped'] = "{:0.2f}".format(
        float(output_dict['Mapped reads']) / float(output_dict['Total reads']) * 100)
else:
    output_dict['% Reads mapped'] = "NA"

if output_dict['Mapped reads'] != '0':
    output_dict['% Mapped reads duplicates'] = "{:0.2f}".format(
        float(output_dict['Duplicate reads']) / float(output_dict['Mapped reads']) * 100)
else:
    output_dict['% Mapped reads duplicates'] = "NA"

# when flagstat of target bam provided
flagstat_target_reads_names = ['OnTarget dupliate reads', 'OnTarget mapped reads', 'OnTarget paired in sequencing','OnTarget read1','OnTarget read2','OnTarget properly paired', 'OnTarget With itself and mate mapped','OnTarget singletons','OnTarget with mate mapped to a different chr','OnTarget with mate mapped to a different chr (mapQ>=5)']

if args.target_flagstat:
    target_flagstat_lines = open(args.target_flagstat, 'r').readlines()
    for i, line in enumerate(target_flagstat_lines[3:]):
        output_dict[flagstat_target_reads_names[i]] = line.split(" ")[0]

else:
    for i in flagstat_target_reads_names:
        output_dict[i] = "NA"
    output_dict['% Reads OnTarget'] = "NA"

# if --rmdup_flagstat provided
if output_dict['Mapped reads'] != '0' and output_dict['OnTarget mapped reads'] != 'NA':
    if args.rmdup_flagstat:
        output_dict['% Reads OnTarget'] = "{:0.2f}".format(
            float(output_dict['OnTarget mapped reads']) /
            (float(output_dict['Mapped reads']) -\
            float(output_dict['Duplicate reads'])) * 100)
    else:
        output_dict['% Reads OnTarget'] = "{:0.2f}".format(
            float(output_dict['OnTarget mapped reads']) /
            float(output_dict['Mapped reads']) * 100)
else:
    output_dict['% Reads OnTarget'] = "NA"
#print(output_dict['% Reads OnTarget'])


if output_dict['OnTarget mapped reads'] != '0' and output_dict['OnTarget mapped reads'] != 'NA' and output_dict['OnTarget dupliate reads'] != 'NA':
    output_dict['% OnTarget duplicates'] = "{:0.2f}".format(
        float(output_dict['OnTarget dupliate reads']) / 
        float(output_dict['OnTarget mapped reads']) * 100)
else:
    output_dict['% OnTarget duplicates'] = "NA"


# coverage
coverage_dict = {}
depth_list, numbases_list = [], []
total_coverage, depthNonzero, depth10, depth20, depth100 = 0, 0, 0, 0, 0
numbaseswithonefifthmean = 0

with open(args.coverage, 'r') as f:
    for line in f:
        # This part located at the end of the bedtools coverage input
        if line.startswith('all') or line.startswith('genome'):
            line = line.split('\t')
            depth = int(line[1])
            numbases = int(line[2])
            numtargetbases = float(line[3])
            coverage_dict[depth] = numbases
            total_coverage += depth * numbases
            depth_list.append(depth)
            numbases_list.append(numbases)
    # Safe guard in case the input is truncated at the end
    if not ( (line.startswith('all')) or (line.startswith('genome')) ):
        sys.exit("Wrong bedtools coverage input, or input is truncated.")

# Calculate mean coverage, 1/5 of coverage
mean_coverage = total_coverage / numtargetbases
output_dict['Mean coverage for target bases'] = "{:0.2f}".format(mean_coverage)
one_fifth_mean = mean_coverage * 0.2

for depth in sorted([i for i in coverage_dict.keys()]):
    if depth >= 1:
        depthNonzero += coverage_dict[depth]
        # coverage at 10-fold, 20-fold, 100-fold or more
        if depth >= 10:
            depth10 += coverage_dict[depth]
            if depth >= 20:
                depth20 += coverage_dict[depth]
                if depth >= 100:
                    depth100 += coverage_dict[depth]
        # 1/5 mean coverage
        if depth >= one_fifth_mean:
            numbaseswithonefifthmean += coverage_dict[depth]

# calculate
output_dict['% Target bases >=1-fold Coverage'] = "{:0.2f}".format(
    depthNonzero / numtargetbases * 100)
output_dict['% Target bases >=10-fold Coverage'] = "{:0.2f}".format(
    depth10 / numtargetbases * 100)
output_dict['% Target bases >=20-fold Coverage'] = "{:0.2f}".format(
    depth20 / numtargetbases * 100)
output_dict['% Target bases >=100-fold Coverage'] = "{:0.2f}".format(
    depth100 / numtargetbases * 100)
output_dict['% Target bases with one fifth of mean coverage'] = "{:0.2f}".format(
    numbaseswithonefifthmean / numtargetbases * 100)

#print(output_dict)
# for k,v in output_dict.items():
#     print(k, v)

if not "/" in args.o:
    f = open(args.o + ".csv", 'w')
else:
    f = safe_open_w(args.o + ".csv")
f.write(','.join(header) + "\n")
f.write(','.join([output_dict[i] for i in header]) + "\n")
f.close()
#print([output_dict[i] for i in header])
