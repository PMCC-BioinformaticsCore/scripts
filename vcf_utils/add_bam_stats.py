#!/usr/bin/python3

"""
Title: add_bam_stats.py
Author: Jiaan Yu
Date: 28-02-2019

Add variant statistics from bam file to vcf.
Works on both germline and somatic variants.
Users are advised to run mpileup tools to generate mpileup file. 
Example samtools command:
samtools mpileup -A -B -Q 0 -d 10000 -l regions.tsv -f reference.fa bam_file > mpileup_file

Examples
For Germline (or Tumour only) data
combine_vcf.py
add_bam_stats.py
  -i [File, vcf]
  --mpileup [File, mpileup output of bam]
  --type 'germline'
  -o [String, output filename]

For Somatic paired data
add_bam_stats.py
  -i [File, vcf]
  --normal_mpileup [File, mpileup output of normal bam]
  --tumor_mpileup [File, mpileup output of tumor bam]
  --tumor_id [String, tumor sample id]
  --type 'somatic'
  -o [String, output filename]

"""

###############################################################################

import re
import argparse
import sys
from variant import Variant, BAM_STATS_LINES


###############################################################################

def create_mpileup_dict(mpileup):
    """Use the mpileup file to create a dictionary
       {'chr\tpos':mpileup_line}
    """
    try:
        f = open(mpileup, 'r')
    except:
        sys.exit('Failed to open file {}'.format(mpileup))
    mpileup_dict = {}
    with open(mpileup, 'r') as f:
        for line in f:
            line = line.split()
            # clean up the start or end of read from mpileup line
            line[4] = re.sub(r'(\^.|\$)', '', line[4])
            mpileup_dict.update({'\t'.join(line[:2]): line[3:5]})
    return mpileup_dict

###############################################################################

# Building API
recognised_modes = ["germline", "somatic"]

# required and optional arguments
parser = argparse.ArgumentParser(description='Get stats from bam file and \
                                 write to vcf')
optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
required.add_argument("-i",  help="input vcf", required=True)
required.add_argument("-o", help="output vcf", required=True)
required.add_argument("--type", help="must be either germline or somatic",
                      required=True, choices=recognised_modes)
parser._action_groups.append(optional)
required.add_argument("--mpileup", help="mpileup file extracted from bam \
                      file", required=False)
required.add_argument("--normal_mpileup", help="mpileup file extracted from the normal sample bam, \
                      required if input is somatic vcf", required=False)
required.add_argument("--tumor_mpileup", help="mpileup file extracted from the tumor sample, \
                      required if input is somatic vcf", required=False)
required.add_argument("--normal_id", help="Normal sample id, \
                      required if input is somatic vcf", required=False)
required.add_argument("--tumor_id", help="Tumor sample id, \
                      required if input is somatic vcf", required=False)
args = parser.parse_args()

# Input sanity check
try:
    f = open(args.i, 'r')
except:
    sys.exit('Failed to open file {}'.format(args.i))
if args.type == 'germline':
    try:
        f = open(args.mpileup, 'r')
    except:
        sys.exit('Failed to open file {}'.format(args.mpileup))
elif args.type == 'somatic':
    if not args.normal_mpileup:
        sys.exit('normal mpileup file is required for somatic vcf')
    elif not args.tumor_mpileup:
        sys.exit('tumor mpileup file is required for somatic vcf')
    elif not args.normal_id:
        sys.exit('normal sample id is required for somatic vcf')
    elif not args.tumor_id:
        sys.exit('tumor sample id is required for somatic vcf')
    else:
        try:
            f = open(args.normal_mpileup, 'r')
        except:
            sys.exit('Failed to open file {}'.format(args.normal_mpileup))
        try:
            f = open(args.tumor_mpileup, 'r')
        except:
            sys.exit('Failed to open file {}'.format(args.tumor_mpileup))


###############################################################################

vcf_in = args.i
vcf_out = args.o
mpileup_in = args.mpileup
vcf_type = args.type
normal_mpileup = args.normal_mpileup
tumor_mpileup = args.tumor_mpileup
normal_id = args.normal_id
tumor_id = args.tumor_id

###############################################################################

if vcf_type == 'germline':

    # Read the mpileup file, and store into dictionary
    # Due to multi-allelic variants, need to read the whole mpileup file in first
    mpileup_dict = create_mpileup_dict(mpileup_in)

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
                mpileup = mpileup_dict.get('\t'.join([variant.chr, variant.pos]),
                                         '')
                bam_stats = variant.cal_bam_stats(mpileup)
                variant = variant.add_bam_stats(bam_stats)
                f_vcf_out.write(variant.write())

else:

    normal_mpileup_dict = create_mpileup_dict(normal_mpileup)
    tumor_mpileup_dict = create_mpileup_dict(tumor_mpileup)

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
                normal_mpileup = normal_mpileup_dict.get('\t'.join([variant.chr,
                                                       variant.pos]), '')
                tumor_mpileup = tumor_mpileup_dict.get('\t'.join([variant.chr,
                                                     variant.pos]), '')
                normal_bam_stats = variant.cal_bam_stats(normal_mpileup)
                tumor_bam_stats = variant.cal_bam_stats(tumor_mpileup)
                variant = variant.add_bam_stats([normal_bam_stats,
                                                tumor_bam_stats],
                                                somatic=True)
                vcf_o.write(variant.write(somatic=True))


