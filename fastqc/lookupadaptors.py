def lookup_adaptors(fastqc_datafile: str, cutadapt_adaptors_lookup: str):
    import mmap, re, csv
    from sys import stderr
    from io import StringIO

    def get_overrepresented_text(fp):
        adapt_section_query = (
            br"(?s)>>Overrepresented sequences\t\S+\n(.*?)>>END_MODULE"
        )

        overrepresented_sequences_match = re.search(adapt_section_query, fp)
        if overrepresented_sequences_match is None:
            raise Exception(
                f"Couldn't find query ('{adapt_section_query.decode('utf8')}') in {fastqc_datafile}"
            )

        return overrepresented_sequences_match.groups()[0].decode("utf8")

    def get_adaptor_ids_from_unparsed_overrpresented_table(tbl):
        rd = csv.reader(StringIO(tbl), delimiter="\t", quotechar='"')
        next(rd)  # discard headers
        return set(r[0] for r in rd)

    def get_cutadapt_map():
        cutadapt_map = {}
        with open(cutadapt_adaptors_lookup) as fp:
            for row in fp:
                st = row.strip()
                if not st or st.startswith("#"):
                    continue
                split = [f for f in st.split("\t") if len(f)]

                if len(split) != 2:
                    print(
                        f"Skipping cutadapt line '{st}' as irregular elements ({len(split)})",
                        file=stderr,
                    )
                    continue

                cutadapt_map[split[1]] = split[0]
        return cutadapt_map

    with open(fastqc_datafile) as f, mmap.mmap(
        f.fileno(), 0, access=mmap.ACCESS_READ
    ) as s:
        text = get_overrepresented_text(s)
        adaptor_ids = get_adaptor_ids_from_unparsed_overrpresented_table(text)

        adaptor_sequences = []

        if adaptor_ids:
            cutadapt_map = get_cutadapt_map()
            for aid in adaptor_ids:
                if aid in cutadapt_map:
                    print(
                        f"Identified sequence '{aid}' as '{cutadapt_map.get(aid)}' in lookup",
                        file=stderr,
                    )
                    adaptor_sequences.append(aid)

                else:
                    print(
                        f"Couldn't find a corresponding sequence for '{aid}' in lookup map",
                        file=stderr,
                    )

        return {"adaptor_sequences": adaptor_sequences}
