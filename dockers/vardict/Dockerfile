# VarDictJava - https://github.com/AstraZeneca-NGS/VarDictJava/

FROM ubuntu:16.04

RUN apt-get update && apt-get install -y \
    curl \
    default-jre \
    r-base

LABEL name="VarDictJava"
LABEL maintainer="michael.franklin@petermac.org"
LABEL build_date="2019-10-14"
LABEL version=1.7.0

WORKDIR /opt

ENV VARDICT_VERSION=1.7.0 BCFTOOLS_VERSION=1.9

ENV VARDICT_INSTALL_DIR /opt/vardict/

RUN curl -LO "https://github.com/AstraZeneca-NGS/VarDictJava/releases/download/${VARDICT_VERSION}/VarDict-${VARDICT_VERSION}.tar"
RUN mkdir ${VARDICT_INSTALL_DIR} \
    && tar -xvf VarDict-${VARDICT_VERSION}.tar \
    && mv VarDict-${VARDICT_VERSION}/* vardict/ \
    && rm VarDict-${VARDICT_VERSION}.tar
  
RUN cd /opt/ \
	&& curl -LO https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
	&& tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
	&& rm -rf bcftoolsbcftools-${BCFTOOLS_VERSION}.tar.bz2  \
	&& cd bcftools-${BCFTOOLS_VERSION}/ \
	&& make && make install

ENV PATH="${VARDICT_INSTALL_DIR}bin:/opt/samtools-${SAMTOOLS_VERSION}/:${PATH}"

CMD ["VarDict"]