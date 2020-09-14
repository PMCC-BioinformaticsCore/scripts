FROM python:3.7
# Run this in the root directory to get the correct build context for vcf_utils/
COPY vcf_utils/ /app/vcf_utils/
COPY fastqc/ /app/fastqc_utils/
COPY performance/ /app/performance/
COPY rna_seq_utils/ /app/rna_seq_utils/
ENV PATH="/app/vcf_utils/:/app/fastqc_utils/:/app/performance/:/app/rna_seq_utils/:${PATH}"
