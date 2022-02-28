# Use R base image and install python 3.7
FROM rocker/r-ver:4.0.2

RUN apt-get update && \
    apt-get -y install zlib1g-dev \
    libxml2 \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    python3.7 python3-pip python3-setuptools

RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/BiocManager/BiocManager_1.30.10.tar.gz",repos = NULL, type = "source")' && \
    R -e 'BiocManager::install("OmicCircos", version="3.12")' && \
    R -e 'BiocManager::install("dplyr", version="3.12")' && \
    R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/matrixStats/matrixStats_0.60.0.tar.gz", repos=NULL, type="source")' && \
    R -e 'BiocManager::install("VariantAnnotation", version="3.12")' && \
    R -e 'BiocManager::install("StructuralVariantAnnotation", version="3.12")' && \
    R -e 'BiocManager::install("rtracklayer", version="3.12")'

COPY circos_plot/ /app/circos_plot/
COPY vcf_utils/ /app/vcf_utils/
COPY fastqc/ /app/fastqc_utils/
COPY performance/ /app/performance/
COPY rna_seq_utils/ /app/rna_seq_utils/

ENV PATH="/app/vcf_utils/:/app/fastqc_utils/:/app/performance/:/app/rna_seq_utils/:/app/circos_plot/:${PATH}"
