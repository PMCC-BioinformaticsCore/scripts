FROM ubuntu:16.04

MAINTAINER Michael Franklin <michael.franklin@petermac.org>

RUN apt-get update -qq \
  && apt-get install -qq bzip2 gcc g++ make zlib1g-dev wget libncurses5-dev liblzma-dev libbz2-dev

ENV SAMTOOLS_VERSION 1.9

LABEL \
  version="${SAMTOOLS_VERSION}" \
  description="Samtools image for use in Workflows"
  
RUN cd /opt/ \
	&& wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2  \
	&& cd samtools-${SAMTOOLS_VERSION}/ \
	&& make && make install

ENV PATH="/opt/samtools-${SAMTOOLS_VERSION}/:${PATH}"