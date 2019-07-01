FROM python:3.7
# Run this in the root directory to get the correct build context for vcf_utils/
COPY vcf_utils/ /app/vcf_utils/
ENV PATH="/app/vcf_utils/:${PATH}"