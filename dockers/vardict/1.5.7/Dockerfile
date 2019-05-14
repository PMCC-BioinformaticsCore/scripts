# VarDictJava - https://github.com/AstraZeneca-NGS/VarDictJava/

FROM ubuntu:16.04

RUN apt-get update && apt-get install -y \
    curl \
    default-jre \
    r-base

LABEL name="VarDictJava"
LABEL maintainer="michael.franklin@petermac.org"
LABEL build_date="2019-04-05"
LABEL version=1.5.7

WORKDIR /opt

ENV VARDICT_VERSION 1.5.7
ENV VARDICT_INSTALL_DIR /opt/vardict/

RUN curl -LO https://github.com/AstraZeneca-NGS/VarDictJava/releases/download/v${VARDICT_VERSION}/VarDict-${VARDICT_VERSION}.tar

RUN mkdir ${VARDICT_INSTALL_DIR} \
    && tar -xvf VarDict-${VARDICT_VERSION}.tar \
    && mv VarDict-${VARDICT_VERSION}/* vardict/ \
    && rm VarDict-${VARDICT_VERSION}.tar

ENV PATH="${VARDICT_INSTALL_DIR}bin:${PATH}"

CMD ["VarDict"]