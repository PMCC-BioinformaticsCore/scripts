"""
Split on line, and write a file for each line as *.bed, where * is the first col
"""
from os import makedirs
from os.path import expanduser, isdir

## OUTPUT LOCATION
output_dir = "~/reference/hg38beds/"




HG38_BEDS = """1	1	248956422
2	1	242193529
3	1	198295559
4	1	190214555
5	1	181538259
6	1	170805979
7	1	159345973
8	1	145138636
9	1	138394717
10	1	133797422
11	1	135086622
12	1	133275309
13	1	114364328
14	1	107043718
15	1	101991189
16	1	90338345
17	1	83257441
18	1	80373285
19	1	58617616
20	1	64444167
21	1	46709983
22	1	50818468
X	1	156040895
Y	1	57227415
M	1	16569"""

output_dir = expanduser(output_dir)
if not isdir(output_dir):
    makedirs(output_dir)

files = []

for l in HG38_BEDS.splitlines():
    fn = output_dir + l.split("\t")[0] + ".bed"

    print("Writing %s ... " % fn, end="\t")

    with open(fn, "w+") as mf:
        mf.write(l)

    print("Done")
