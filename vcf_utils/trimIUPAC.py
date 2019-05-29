import sys
import re
import os

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


with open(inpath) as inputfp, open(outpath,'w') as outputfp:
    extrahead_elems = []


    indRef = None   # header.index("REF")
    indAlt = None   # header.index("ALT")

    p = re.compile("W|K|Y|R|S|M")
    for line in inputfp:

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
            processed[indRef]=p.sub("N", processed[indRef])
            processed[indAlt]=p.sub("N", processed[indAlt])

        outputfp.write("\t".join(processed) + "\n")
