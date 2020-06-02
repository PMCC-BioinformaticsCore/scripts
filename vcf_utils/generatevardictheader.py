"""
Authors:
    - Michael Franklin (author of this converted script) [2020-06-02]
    - Jiaan Yu (author of original script)

This is integrated into Janis as "GenerateVardictHeaderLines"

Converted from script:

    out="header_file.txt"
    #Additional line for vardict source
    echo -e "##source=vardict" > $out 
    # Parse fasta dictionary
    while read line; do
    if [[ "$line" == @SQ* ]]; then
        #echo $line
        IFS=$'\t' read -ra array <<< "$line"
        chrom="${array[1]}"; 
        chrom="${chrom/SN:/}"
        length="${array[2]}"; 
        length="${length/LN:/}"
        echo -e "##contig=<ID=${chrom},length=${length}>" >> $out
        #echo $chrom $length
    fi
    done < $1
"""

def generate_vardictheader_from_fasta_dict(fasta_dict: str, output_filename: str="out.txt"):
    with open(output_filename, "w+") as out, open(fasta_dict) as inp:
        out.write("##source=vardict\n")
        for line in inp:
            if not line.startswith("@SQ"):
                continue
            pieces = line.split("\t")
            chrom = pieces[1].replace("SN:", "")
            length = pieces[2].replace("LN:", "")

            out.write(f"##contig=<ID={chrom},length={length}>\n")
