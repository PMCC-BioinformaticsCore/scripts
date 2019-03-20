#!/usr/bin/python3

"""
Title: add_bam_stats.py
Author: Jiaan Yu
Date: 28-02-2019

Add variant statistics from bam file to vcf.
Works on both germline and somatic variants.
Users are advised to run pileup tools to generate pileup file. 
An example to run: 
samtools mpileup -A -B -Q 0 -d 10000 -l regions.tsv -f reference.fa bam_file > pileup_file
"""

###############################################################################

import re
import argparse
import sys
from variant import Variant, BAM_STATS_LINES


###############################################################################

def create_pileup_dict(pileup):
    """Use the pileup file to create a dictionary
       {'chr\tpos':pileup_line}
    """
    try:
        f = open(pileup, 'r')
    except:
        sys.exit('Failed to open file {}'.format(pileup))
    pileup_dict = {}
    with open(pileup, 'r') as f:
        for line in f:
            line = line.split()
            # clean up the start or end of read from pileup line
            line[4] = re.sub(r'(\^.|\$)', '', line[4])
            pileup_dict.update({'\t'.join(line[:2]): line[3:5]})
    return pileup_dict

###############################################################################

# Building API

# required and optional arguments
parser = argparse.ArgumentParser(description='Get stats from bam file and \
                                 write to vcf')
optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
required.add_argument("-i",  help="input vcf", required=True)
required.add_argument("-o", help="output vcf", required=True)
required.add_argument("--type", help="must be either germline or somatic",
                      required=True)
parser._action_groups.append(optional)
required.add_argument("--pileup", help="Pileup file extracted from bam \
                      file", required=False)
required.add_argument("--normal_pileup", help="Pileup file extracted from the normal sample bam, \
                      required if input is somatic vcf", required=False)
required.add_argument("--tumor_pileup", help="Pileup file extracted from the tumor sample, \
                      required if input is somatic vcf", required=False)
required.add_argument("--normal_id", help="Normal sample id, \
                      required if input is somatic vcf", required=False)
required.add_argument("--tumor_id", help="Tumor sample id, \
                      required if input is somatic vcf", required=False)
args = parser.parse_args()

# Input sanity check
if not args.i.endswith('.vcf'):
        sys.exit('{} is not a vcf'.format(args.i))
else:
    try:
        f = open(args.i, 'r')
    except:
        sys.exit('Failed to open file {}'.format(args.i))
if args.type != 'germline' and args.type != 'somatic':
    sys.exit('vcf type not recognised')
elif args.type == 'germline':
    if not args.pileup:
        sys.exit('pileup file not specified for germline vcf')
    elif not args.pileup.endswith('.pileup'):
        sys.exit('{} is not in pileup format'.format(args.pileup))
    else:
        try:
            f = open(args.pileup, 'r')
        except:
            sys.exit('Failed to open file {}'.format(args.pileup))
elif args.type == 'somatic':
    if not args.normal_pileup:
        sys.exit('normal pileup file is required for somatic vcf')
    elif not args.normal_pileup.endswith('.pileup'):
        sys.exit('{} is not in pileup format'.format(args.normal_pileup))
    elif not args.tumor_pileup:
        sys.exit('tumor pileup file is required for somatic vcf')
    elif not args.tumor_pileup.endswith('.pileup'):
        sys.exit('{} is not in pileup format'.format(args.tumor_pileup))
    elif not args.normal_id:
        sys.exit('normal sample id is required for somatic vcf')
    elif not args.tumor_id:
        sys.exist('tumor sample id is required for somatic vcf')


###############################################################################

vcf_in = args.i
vcf_out = args.o
pileup_in = args.pileup
vcf_type = args.type
normal_pileup = args.normal_pileup
tumor_pileup = args.tumor_pileup
normal_id = args.normal_id
tumor_id = args.tumor_id

###############################################################################

if vcf_type == 'germline':

    # Read the pileup file, and store into dictionary
    # Due to multi-allelic variants, need to read the whole pileup file in first
    pileup_dict = create_pileup_dict(pileup_in)

    with open(vcf_in, 'r') as f_vcf, open(vcf_out, 'w') as f_vcf_out:
        for line in f_vcf:

            # Write the header lines
            if line.startswith('##'):
                f_vcf_out.write(line)
            elif line.startswith('#'):
                for new_line in BAM_STATS_LINES:
                    f_vcf_out.write(new_line)
                f_vcf_out.write(line)

            # Calcualte bam stats for variants
            else:
                variant = Variant().read_variant(line)
                pileup = pileup_dict.get('\t'.join([variant.chr, variant.pos]),
                                         '')
                bam_stats = variant.cal_bam_stats(pileup)
                variant = variant.add_bam_stats(bam_stats)
                f_vcf_out.write(variant.write())

else:

    normal_pileup_dict = create_pileup_dict(normal_pileup)
    tumor_pileup_dict = create_pileup_dict(tumor_pileup)

    with open(vcf_in, 'r') as vcf_i, open(vcf_out, 'w') as vcf_o:
        for line in vcf_i:
            if line.startswith('#'):
                if line.startswith('##'):
                    vcf_o.write(line)
                else:
                    for new_line in BAM_STATS_LINES:
                        vcf_o.write(new_line)
                    line = line.split()
                    try:
                        normal_index = line.index(normal_id)
                    except:
                        sys.exit('Failed to match normal sample id')
                    try:
                        tumor_index = line.index(tumor_id)
                    except:
                        sys.exit('Failed to match tumor sample id')
                    line = '\t'.join(line[:-2] + [normal_id, tumor_id+'\n'])
                    vcf_o.write(line)
            else:
                variant = Variant().read_variant(line, somatic=True,
                                                 normal=normal_index,
                                                 tumor=tumor_index)
                normal_pileup = normal_pileup_dict.get('\t'.join([variant.chr,
                                                       variant.pos]), '')
                tumor_pileup = tumor_pileup_dict.get('\t'.join([variant.chr,
                                                     variant.pos]), '')
                normal_bam_stats = variant.cal_bam_stats(normal_pileup)
                tumor_bam_stats = variant.cal_bam_stats(tumor_pileup)
                variant = variant.add_bam_stats([normal_bam_stats,
                                                tumor_bam_stats],
                                                somatic=True)
                vcf_o.write(variant.write(somatic=True))
