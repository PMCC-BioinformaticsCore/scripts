FROM ubuntu:16.04

MAINTAINER Michael Franklin <michael.franklin@petermac.org>

RUN apt-get update -qq \
  && apt-get install -qq bzip2 gcc g++ make zlib1g-dev wget libncurses5-dev liblzma-dev libbz2-dev

ENV BCFTOOLS_VERSION 1.7

LABEL \
  version="${BCFTOOLS_VERSION}" \
  description="BCFTools image for use in Workflows"
  
RUN cd /opt/ \
	&& wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
	&& tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
	&& rm -rf bcftoolsbcftools-${BCFTOOLS_VERSION}.tar.bz2  \
	&& cd bcftools-${BCFTOOLS_VERSION}/ \
	&& make && make install

ENV PATH="/opt/bcftools-${BCFTOOLS_VERSION}/:${PATH}"