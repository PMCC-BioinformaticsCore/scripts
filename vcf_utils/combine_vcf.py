#!/usr/bin/python3

"""
Title: combine_vcf
Author: Jiaan Yu
Date: 28-02-2019
Modified: 29-05-2019 (Michael Franklin)

Extract (and normalise) information from input vcfs coming from different
callers, combine to one vcf, and calcualte mean and standard deviation of AD and DP if specified.
Works on germline and somatic vcfs.
Warning: the output vcf is not sorted by chromosome and position, users are advised to use other tools to sort this vcf. (Example: bcftools sort vcf -o sorted_vcf)
         Tested on HaplotypeCaller and Mutect2 (gatk 4.0.10.0), strelka (2.9.2) and vardict (1.5.1)

Notes on how the columns are being parsed:
-- AD/DP: INFO value overwritten by FORMAT value (if FORMAT value exists), in the process_(somatic_)variant
-- GT: The only enfored column in FORMAT
-- Other columns: extract_cols moves FORMAT column names to INFO; select_info finds column values from variants
"""

###############################################################################

import os
import sys
<<<<<<< HEAD
import argparse
from collections import OrderedDict, defaultdict as dd
=======
from collections import OrderedDict, defaultdict as dd, Counter
from itertools import combinations
>>>>>>> c1889550c3d9c941432671dd16940883b91fd285
from normalisedvcf import NormalisedVcf, sort_vcf
from variant import Variant
from vcfheader import STATS_HEADER, SOMATIC_STATS_HEADER, HEADER

###############################################################################

# Building API
recognised_modes = ["germline", "somatic"]


# required and optional arguments
parser = argparse.ArgumentParser(description="Extracts and combines the \
                                 information from germline / somatic vcfs \
                                 into one")
optional = parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
required.add_argument("-i",  help="input vcfs, the priority of the vcfs will be \
                      based on the order of the input. This parameter can be \
                      specified more than once", action="append", required=True)
required.add_argument("--columns", help="Columns to keep. This parameter can \
                      be specified more than once", action="append", required=True)
required.add_argument("-o", help="output vcf (unsorted)", required=True)
required.add_argument("--type", help="must be either germline or somatic", required=True,
                      choices=recognised_modes)
parser._action_groups.append(optional)
required.add_argument("--regions", help="Region file containing all the \
                      variants, used as samtools mpileup", required=False)
required.add_argument("--normal", help="Sample id of germline vcf, or normal \
                      sample id of somatic vcf", required=False)
required.add_argument("--tumor", help="tumor sample ID, required if inputs are \
                      somatic vcfs", required=False)
required.add_argument("--priority", help="The priority of the callers, must match \
                      with the callers in the source header", nargs="+", required=False)
args = parser.parse_args()

# Sanity check number of inputs
if len(args.i) < 1:
    sys.exit("There must be at least 1 VCF file")


# Check to make sure all of the files exist
input_files = args.i

invalid_files = [f for f in input_files if not os.path.exists(f)]
has_invalid_files = len(invalid_files) > 0
if has_invalid_files:
    ninvalid = len(invalid_files)
    nall = len(input_files)

    sys.exit("There were files that could not be found on disk (%s/%s): %" % (ninvalid, nall, ", ".join(invalid_files)))

has_duplicates = len(input_files) != len(set(input_files))

if has_duplicates:
    sys.exit("There are duplicates in the input vcfs, please specify again")
if args.priority and len(input_files) != len(args.priority):
    sys.exit("The number of vcfs (%s) does not match with the number of callers in priority (%s)" % (len(input_files), len(args.priority)))

if args.type not in recognised_modes:
    sys.exit("The mode '{s}' was not recognised, must be one of: %s" % (args.type, ", ".join(recognised_modes)))

if args.type == "somatic" and not (args.normal and args.tumor):
    sys.exit("normal and tumor ids are required for somatic vcfs")

###############################################################################

vcf_in = args.i
columns_to_keep = args.columns
vcf_out = args.o
regions = args.regions
vcf_type = args.type
normal_id = args.normal
tumor_id = args.tumor

###############################################################################

if vcf_type == "germline":

    # A list of cleaned vcf with extracted columns
    vcf_list = [NormalisedVcf(vcf).process_vcf(columns_to_keep) for vcf
                in vcf_in]

    # Combine the variants into a list
    combined_variants, variant_to_vcf_dict = [], dd(list)
    callers = [vcf.caller for vcf in vcf_list]

    # Sort the vcf
    vcf_list, callers = sort_vcf(args.priority, callers, vcf_list)

    for i, vcf in enumerate(vcf_list):
        combined_variants += list(vcf.variants.keys())
        # Dictionary that let varaint refer back to vcf
        for var in vcf.variants.keys():
            variant_to_vcf_dict[var].append(i)

    # Take the unquie variants, don"t sort them
    combined_variants = list(set(combined_variants))

    # Take the unquie header lines with preserved order
    # To do: Do we want to reorder the header lines?
    #        Protentially contain redundant lines
    combined_header = []
    for vcf in vcf_list:
        combined_header += vcf.meta_info

    # Write the combined vcf
    with open(vcf_out, "w") as combined_f:

        # Write the meta info and header lines
        for line in list(OrderedDict.fromkeys(combined_header)):
            combined_f.write(line)
        for line in STATS_HEADER:
            combined_f.write(line)
        # write the chr\tpos\t... line
        combined_f.write(vcf_list[0].header)

        # Write the variants
        for v_key in combined_variants:
            callers_indexes = variant_to_vcf_dict[v_key]
            callers_names = [callers[i] for i in callers_indexes]
            info_dict = OrderedDict()
            for i in callers_indexes:
                # Combine the selected information in the dictionary
                info_dict.update(vcf_list[i].variants[v_key].info)
            i = callers_indexes[0]
            combined_variant = Variant.combine_info(
                               vcf_list[i].variants[v_key], columns_to_keep,
                               callers_names, info_dict)
            combined_f.write(Variant.write(combined_variant))

else:

    # Process each vcf and extract the information from the selected columns
    vcf_list = [NormalisedVcf(vcf).process_somatic_vcf(columns_to_keep,
                normal_id, tumor_id) for vcf in vcf_in]

    # Combine the variants into a list
    combined_variants, variant_to_vcf_dict = [], dd(list)
    callers = [vcf.caller for vcf in vcf_list]

    # Sort the vcf
    vcf_list, callers = sort_vcf(args.priority, callers, vcf_list)

    for i, vcf in enumerate(vcf_list):
        combined_variants += list(vcf.variants.keys())
        # Dictionary that let varaint refer back to vcf
        for var in vcf.variants.keys():
            variant_to_vcf_dict[var].append(i)

    combined_variants = list(set(combined_variants))

    combined_header = []
    for vcf in vcf_list:
        combined_header += vcf.meta_info

    # Write the combined vcf
    with open(vcf_out, "w") as combined_f:

        # Write the meta info and header lines
        for line in list(OrderedDict.fromkeys(combined_header)):
            combined_f.write(line)
        for line in SOMATIC_STATS_HEADER:
            combined_f.write(line)
        # write the chr\tpos\t... line
        combined_f.write("\t".join([HEADER, normal_id, tumor_id+"\n"]))

        # Write the variants
        for v_key in combined_variants:
            callers_indexes = variant_to_vcf_dict[v_key]
            callers_names = [callers[i] for i in callers_indexes]
            info_dict = OrderedDict()
            for i in callers_indexes:
                info_dict.update(vcf_list[i].variants[v_key].info)
            i = callers_indexes[0]
            combined_variant = Variant.combine_info(
                               vcf_list[i].variants[v_key], columns_to_keep,
                               callers_names, info_dict, somatic=True)
            combined_f.write(Variant.write(combined_variant, somatic=True))

# Output combine varaints summary
# Count the number of variants 
summary = Counter([tuple(v) for v in variant_to_vcf_dict.values()])
summary_count = []
# Output the combinations of the callers
for i, caller in enumerate(callers):
    if i == 0:
        summary_header = [caller for caller in callers]
    else:
        for j in list(combinations(callers, i+1)):
            summary_header.append('-'.join(j))
for i in range(len(callers)+1):
    for j in list(combinations(range(len(callers)), i+1)):
        summary_count.append(str(summary[j]))
with open("Combine_variants_summary.tsv", "w") as f:
    f.write('\t'.join(summary_header) + '\n')
    f.write('\t'.join(summary_count) + '\n')

# Write variant location file for samtools pileup
if regions:
    with open(regions, "w") as loc_f:
        for v_key in combined_variants:
            loc_f.write("\t".join(v_key.split()[:2]) + "\n")
