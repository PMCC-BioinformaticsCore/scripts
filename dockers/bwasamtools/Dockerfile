# Merge BWA and SAMTOOLS together to allow us to pipe bwa mem to samtools view

FROM ubuntu:16.04

RUN apt-get update -qq \
  && apt-get install -qq bzip2 gcc g++ make zlib1g-dev wget libncurses5-dev liblzma-dev libbz2-dev

ENV BWA_VERSION 0.7.17
ENV SAMTOOLS_VERSION 1.9


RUN cd /opt/ \
	&& wget https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2 \
	&& tar -xjf bwa-${BWA_VERSION}.tar.bz2 \
	&& rm -f bwa-${BWA_VERSION}.tar.bz2 \
	&& cd /opt/bwa-${BWA_VERSION}/ \
	&& make

RUN cd /opt/ \
	&& wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2  \
	&& cd samtools-${SAMTOOLS_VERSION}/ \
	&& make && make install

ENV PATH="/opt/bwa-${BWA_VERSION}/:/opt/samtools-${SAMTOOLS_VERSION}/:${PATH}"