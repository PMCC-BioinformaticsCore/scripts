#!/usr/bin/python3
# Based on github repo:
# https://github.com/J35P312/MegaFusion
import argparse
import json


def retrieve_required_entry(data, tag, content):
    if type(data["required"][tag]) is dict:
        info = content[data["required"][tag]["column"]].split(
            data["required"][tag]["delimiter"]
        )[data["required"][tag]["element"]]
    else:
        info = content[data["required"][tag]]
    return info


version = "0.0.2"
parser = argparse.ArgumentParser(
    """MegaFusion.py-{}: convert gene fusion tab file to SV vcf""".format(version)
)
parser.add_argument(
    "--fusion",
    required=True,
    type=str,
    help="path to fusion tab file",
)
parser.add_argument(
    "--json",
    required=True,
    type=str,
    help="path to a config json file",
)
parser.add_argument(
    "--sample",
    type=str,
    default="Bob",
    help="Sample name (default=Bob)",
)
parser.add_argument(
    "--contig",
    required=True,
    type=str,
    help="path to a contig file to add to vcf",
)
parser.add_argument(
    "--tool_version",
    required=True,
    type=str,
    help="version of the fusion tool",
)
args = parser.parse_args()
args.version = version

with open(args.json) as json_file:
    data = json.load(json_file)

# start to print headers
print("##fileformat=VCFv4.3")
print("##source={}_v{}".format(data["source"], args.tool_version))
print('##ALT=<ID=BND,Description="Break end">')
print(
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural \
variant">'
)
print(
    '##INFO=<ID=ROWID,Number=1,Type=String,Description="Breakend ID of the \
breakend at reference position">'
)
print(
    '##INFO=<ID=MATEID,Number=1,Type=String,Description="Breakend ID of the \
mate breakend">'
)
print(
    '##INFO=<ID=EVENT,Number=1,Type=String,Description="Event identifier for \
a rearrangement">'
)
print(
    '##INFO=<ID=VARTYPE,Number=1,Type=String,Description="Type of \
variant/arrangement">'
)
print(
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the \
variant described in this record">'
)
print('##INFO=<ID=GENEA,Number=1,Type=String,Description="Gene A">')
print('##INFO=<ID=GENEB,Number=1,Type=String,Description="Gene B">')
print(
    '##INFO=<ID=ORIENTATION,Number=1,Type=String,Description="Orientation of \
the fusion">'
)
print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')


# define the INFO and FORMAT columns
required_INFO = "SVTYPE=BND;ROWID={};MATEID={};EVENT=AF{};VARTYPE=fusion;GENEA={};GENEB={};ORIENTATION={},{}"
required_format = "GT"
vcf_columns = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
required_sample = "0/1"


# add custom columns to INFO or FORMAT
for entry in sorted(data["custom"].keys()):
    try:
        print(
            '##{}=<ID={},Number={},Type={},Description="{}">'.format(
                data["custom"][entry]["entry"],
                entry,
                data["custom"][entry]["number"],
                data["custom"][entry]["type"],
                data["custom"][entry]["description"],
            )
        )
        if "FORMAT" == data["custom"][entry]["entry"]:
            required_format += ":{}".format(entry)
    except:
        print("invalid config")
        quit()

# check if the json file has the required columns
checklist = [
    "chromosomeA",
    "chromosomeB",
    "posA",
    "posB",
    "strandA",
    "strandB",
    "GeneA",
    "GeneB",
]
for item in data["required"].keys():
    if item in checklist:
        pass
    else:
        print("missing required entry: {}".format(item))
        quit()

# filter the variant if specified
filter_string = '##FILTER=<ID={},Description="Not up to snuff">'
filter_set = set([])
if data["filter"]["filter"]:
    for line in open(args.fusion):
        if line.startswith(data["header"]):
            continue
        content = line.strip().split(data["delimiter"])
        if data["filter"]["pass"] == content[data["filter"]["column"]]:
            continue

        if not "delimiter" in data["filter"]:
            if not content[data["filter"]["column"]] in filter_set:
                print(filter_string.format(content[data["filter"]["column"]]))
            filter_set.add(content[data["filter"]["column"]])
        else:
            for entry in content[data["filter"]["column"]].split(
                data["filter"]["delimiter"]
            ):
                if not entry in filter_set:
                    print(filter_string.format(entry))
            filter_set.add(entry)

# print contigs to the header
with open(args.contig, "r") as contig_file:
    for line in contig_file:
        print(line.strip("\n"))

# print the last line of header
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(args.sample))


i = 1
for line in open(args.fusion):
    if line.startswith(data["header"]):
        continue

    content = line.strip().split(data["delimiter"])

    chrA = retrieve_required_entry(data, "chromosomeA", content)
    chrB = retrieve_required_entry(data, "chromosomeB", content)
    posA = retrieve_required_entry(data, "posA", content)
    posB = retrieve_required_entry(data, "posB", content)
    geneA = retrieve_required_entry(data, "GeneA", content)
    geneB = retrieve_required_entry(data, "GeneB", content)
    strandA = retrieve_required_entry(data, "strandA", content)
    strandB = retrieve_required_entry(data, "strandB", content)

    # get ID
    IDA = "{}_Fusion_{}A".format(data["source"], i)
    IDB = "{}_Fusion_{}B".format(data["source"], i)

    # get ALT
    if strandA == strandB == "+":
        altA = "N[{}:{}[".format(chrB, posB)
        altB = "]{}:{}]N".format(chrA, posA)
    elif strandA == strandB == "-":
        altA = "]{}:{}]N".format(chrB, posB)
        altB = "N[{}:{}[".format(chrA, posA)
    elif strandA == "+" and strandB == "-":
        altA = "N]{}:{}]".format(chrB, posB)
        altB = "N]{}:{}]".format(chrA, posA)
    elif strandA == "-" and strandB == "+":
        altA = "[{}:{}[N".format(chrB, posB)
        altB = "[{}:{}[N".format(chrA, posA)
    else:
        altA, altB = ".", "."

    # get INFO
    INFOA = required_INFO.format(IDA, IDB, i, geneA, geneB, strandA, strandB)
    INFOB = required_INFO.format(IDB, IDA, i, geneA, geneB, strandA, strandB)

    # get FORMAT
    FORMAT = required_sample

    # get custom INFO/FORMAT
    for entry in data["custom"]:
        if data["custom"][entry]["entry"] == "INFO":
            if "none" in data["custom"][entry]:
                if (
                    data["custom"][entry]["none"]
                    == content[data["custom"][entry]["column"]]
                ):
                    continue
            if "remove" in data["custom"][entry]:
                for r in data["custom"][entry]["remove"]:
                    content[data["custom"][entry]["column"]] = content[
                        data["custom"][entry]["column"]
                    ].replace(r, "")

            INFOA += ";{}={}".format(entry, content[data["custom"][entry]["column"]])
            INFOB += ";{}={}".format(entry, content[data["custom"][entry]["column"]])

    for entry in sorted(data["custom"].keys()):
        if data["custom"][entry]["entry"] == "FORMAT":
            if "none" in data["custom"][entry]:
                if (
                    data["custom"][entry]["none"]
                    == content[data["custom"][entry]["column"]]
                ):
                    content[data["custom"][entry]["column"]] = "."
            if "remove" in data["custom"][entry]:
                for r in data["custom"][entry]["remove"]:
                    entry = entry.replace(r, "")
            FORMAT += ":{}".format(content[data["custom"][entry]["column"]])

    # get QUAL
    qual = "."

    # get FILTER (and apply filter if specified)
    filt = "PASS"
    if data["filter"]["filter"]:
        if content[data["filter"]["column"]] == data["filter"]["pass"]:
            pass
        else:
            if not "delimiter" in data["filter"]:
                qual = content[data["filter"]["column"]]
            else:
                qual = content[data["filter"]["column"]].replace(
                    data["filter"]["delimiter"], ","
                )

    print(
        vcf_columns.format(
            chrA, posA, ".", "N", altA, qual, filt, INFOA, required_format, FORMAT
        )
    )
    print(
        vcf_columns.format(
            chrB, posB, ".", "N", altB, qual, filt, INFOB, required_format, FORMAT
        )
    )
    i += 1
