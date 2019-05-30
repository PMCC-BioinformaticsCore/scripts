#!/usr/bin/python3

"""
Title: combine_vcf
Author: Jiaan Yu / Michael Franklin
Date: 29-05-2019
Modified: 30-05-2019 (Michael Franklin)

Substitute 'W|K|Y|R|S|M' with 'N'
"""

###############################################################################

import sys
import re
import os
import json

# We'll filter out the characters listed in this set:
# https://en.wikipedia.org/wiki/Nucleic_acid_notation
character_set = ["W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "Z"]

USAGE="USAGE: %s <inVcf> <outVcf>\n\nTake out ambiguity (IUPAC) codes in REF and ALT columns by converting them to Ns.\n\n" % sys.argv[0]

if len(sys.argv) <= 2:
    print(USAGE + "\033[91mError: An invalid number of parameters was provided (%s / 2).\033[0m\n" % (len(sys.argv) -1))
    sys.exit(1)

inpath = sys.argv[1]
outpath = sys.argv[2]

if not os.path.exists(inpath):
    print("\033[91mThe input file '%s' could not be found, or is not accessible\033[0m" % inpath)
    sys.exit(1)

if os.path.exists(outpath):
    print("\033[91mA file already exists at the output path '%s'\033[0m" % outpath)
    sys.exit(1)

replacements = 0
replacement_dict = {
    # $line: { indRef: prev, indAlt: prev }
}

line_number = 0

with open(inpath) as inputfp, open(outpath,'w') as outputfp:
    extrahead_elems = []

    indRef = None   # header.index("REF")
    indAlt = None   # header.index("ALT")

    p = re.compile("|".join(character_set))
    for line in inputfp:
        line_number += 1    # was started at 0, so lines will start at 1

        if line.startswith("##"):
            # just stream the header to new file
            outputfp.write(line)
            continue

        # probably write extrahead params here

        processed = line.rstrip("\n\r").split("\t")

        if indRef is None or indAlt is None:
            indRef = processed.index("REF")
            indAlt = processed.index("ALT")
        else:
            lIndRef, lIndAlt = processed[indRef], processed[indAlt]
            has_indRefMatch = p.match(lIndRef)
            has_indAltMatch = p.match(lIndAlt)
            if has_indRefMatch or has_indAltMatch:
                d = {}
                if has_indRefMatch:
                    d["indRef"] = lIndRef
                    processed[indRef] = p.sub("N", str(lIndRef))

                if has_indAltMatch:
                    d["indAlt"] = lIndAlt                
                    processed[indAlt] = p.sub("N", str(lIndAlt))

                replacement_dict[line_number] = d
                replacements += len(d)

        outputfp.write("\t".join(processed) + "\n")

with open("stats.json", "w+") as stats:
    json.dump({
        "total": replacements,
        "replacements": replacement_dict
    }, stats)

