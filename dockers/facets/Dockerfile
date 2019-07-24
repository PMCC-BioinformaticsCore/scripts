FROM ubuntu:16.04

# Based on the Dockerfile papaemmelab/docker-facets
# This file contains 

MAINTAINER Michael Franklin <michael.franklin@petermac.org>

ENV VERSION_PCTGCDATA=v0.2.0
ENV VERSION_FACETS=v0.5.14
ENV HTSLIB_URL=https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2

LABEL \
  version="${VERSION_FACETS}" \
  description="Algorithm to implement Fraction and Copy number Estimate from Tumor/normal Sequencing."

RUN \
    apt-get update && \
    apt-get -y --fix-missing install \
        build-essential \
        bzip2 \
        libbz2-dev \
        libcurl4-gnutls-dev \
        liblzma-dev \
        libxml2-dev \
        r-base-core \
        zlib1g-dev \
        git \
        libssl-dev

RUN \
    # install devtools and argparse
    Rscript -e 'install.packages("argparse", repos="http://cran.us.r-project.org")' && \
    Rscript -e 'install.packages("devtools", repos="http://cran.us.r-project.org")' && \
    \
    # install specific version of pctGCdata
    git clone https://github.com/mskcc/pctGCdata /pctGCdata && \
    cd /pctGCdata && git checkout tags/${VERSION_PCTGCDATA} && \
    Rscript -e 'devtools::install_local("/pctGCdata", build_vignettes=FALSE)' && \
    \
    # install specific version of facets
    git clone https://github.com/mskcc/facets /facets && \
    cd /facets && git checkout tags/${VERSION_FACETS} && \
    Rscript -e 'devtools::install_local("/facets", build_vignettes=FALSE)'

RUN \
    cd /tmp && \
    curl -sSL -o tmp.tar.gz --retry 10 ${HTSLIB_URL} && \
    tar -xjf tmp.tar.gz && \
    cd htslib-1.9 && ./configure && make && make install && \
    cd /facets/inst/extcode && \
    g++ \
        -std=c++11 \
        -I/usr/local/include/htslib snp-pileup.cpp \
        -L/usr/local/lib \
        -lhts \
        -Wl,-rpath=/usr/local/lib \
        -o /bin/snp-pileup

COPY build/run_facets.R /facets/run_pmac_facets.R
VOLUME [ "/gpfs", "/ifs", "/data" ]
ENTRYPOINT [ "Rscript", "/facets/run_pmac_facets.R" ]
