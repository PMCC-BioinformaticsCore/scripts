# VCF Utils

In this README, we'll document some of the tools present in this repository.

## Dockerized

The tools in this repository have been [dockerized](https://hub.docker.com/r/michaelfranklin/pmacutil) in:

```bash
docker pull michaelfranklin/pmacutil
```

Autobuilds have been configured, this can be triggered by tagging a commit with the format `pmacutil-0.x.x` and pushing to Github. This will automatically create a tag `0.x.x`. For example:
```bash
git tag -a "pmacutil-0.0.1" -m "<Message for release>"
```

The tools have been added to the path, so you can simply call the tool from the command line.

## Tools

### Combine VCFs
File: `combine_vcf.py`

Usage:
```
usage: combine_vcf.py [-h] -i I --columns COLUMNS -o O --type
                      {germline,somatic} [--regions REGIONS] [--normal NORMAL]
                      [--tumor TUMOR] [--priority PRIORITY [PRIORITY ...]]

Extracts and combines the information from germline / somatic vcfs into one

required arguments:
  -i I                  input vcfs, the priority of the vcfs will be based on
                        the order of the input. This parameter can be
                        specified more than once
  --columns COLUMNS     Columns to keep. This parameter can be specified more
                        than once
  -o O                  output vcf (unsorted)
  --type {germline,somatic}
                        must be either germline or somatic
  --regions REGIONS     Region file containing all the variants, used as
                        samtools mpileup
  --normal NORMAL       Sample ID of germline, or tumor sample ID  of
                        somatic vcf
  --tumor TUMOR         tumor sample ID, required if inputs are somatic vcfs
  --priority PRIORITY [PRIORITY ...]
                        The priority of the callers, must match with the
                        callers in the source header

optional arguments:
  -h, --help            show this help message and exit
  --priority PRIORITY [PRIORITY ...]
                        The priority of the callers, must match with the
                        callers in the source header. Default: sort by alphabetical 
                        order
  --columns COLUMNS
                        Columns to be extracted, seperated by comma
 ```


### Add Bam Stats
File: `add_bam_stats.py`

Usage:
```
usage: add_bam_stats.py [-h] -i I -o O --type TYPE [--mpileup MPILEUP]
                        [--normal_mpileup NORMAL_MPILEUP]
                        [--tumor_mpileup TUMOR_MPILEUP]
                        [--normal_id NORMAL_ID] [--tumor_id TUMOR_ID]

Get stats from bam file and write to vcf

required arguments:
  -i I                  input vcf
  -o O                  output vcf
  --type TYPE           must be either germline or somatic
  --mpileup MPILEUP     mpileup file extracted from bam file
  --normal_mpileup NORMAL_MPILEUP
                        mpileup file extracted from the normal sample bam,
                        required if input is somatic vcf
  --tumor_mpileup TUMOR_MPILEUP
                        mpileup file extracted from the tumor sample, required
                        if input is somatic vcf
  --normal_id NORMAL_ID
                        Normal sample id, required if input is somatic vcf
  --tumor_id TUMOR_ID   Tumor sample id, required if input is somatic vcf

optional arguments:
  -h, --help            show this help message and exit

```


## List of files:

Scripts:
- combine_vcf.py: Combine multiple vcfs (called from the same sample) into one. Variant callers are marked with "Identified=" in INFO. Offer options to select columns (INFO/FORMAT) to include.
- add_bam_stats.py: Add variants statistics from bam file (mpileup)

## Example usage


Supporting classes:
- `normalisedvcf.py`: For parsing vcfs
- `variant.py`: For parsing variants
- `vcfheader.py`: For parsing headers

Usage examples:

- Germline VCFs

```bash
python3 combine_vcf.py -i vcf1 -i vcf2 -i vcf3 --columns AD,DP,AF,GT -o combined.vcf --type germline --regions regions.tsv
samtools mpileup -A -B -Q 0 -d 10000 -l regions.tsv -f reference.fa sample.bam > regions.mpileup
bcftools sort combined.vcf -o combined.sorted.vcf
python3 add_bam_stats.py -i combined.sorted.vcf -o combined.sorted.addbamstats.vcf --type germline -mpileup sample.mpileup
```

- For somatic vcfs

```bash
python combine_vcf.py -i vcf1 -i vcf2 -i vcf3 --columns AD,DP,GT -o combined.vcf --type somatic --regions regions.tsv --normal_id normal --tumor_id tumor
samtools mpileup -A -B -Q 0 -d 10000 -l regions.tsv -f reference.fa sample_normal.bam > normal.mpileup
samtools mpileup -A -B -Q 0 -d 10000 -l regions.tsv -f reference.fa sample_tumor.bam > tumor.mpileup
bcftools sort combined.vcf -o combined.sorted.vcf python3 add_bam_stats.py -i combined.sorted.vcf -o combined.sorted.addbamstats.vcf --type somatic --normal_id normal --tumor_id tumor --normal_mpileup normal.mpileup --tumor_mpileup tumor.mpileup
```
