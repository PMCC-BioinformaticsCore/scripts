# Based on https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/installation.md

FROM ubuntu:16.04
MAINTAINER Michael Franklin <michael.franklin@petermac.org>

ENV STRELKA_VERSION 2.9.9

LABEL \
  version="2.9.9" \
  description="Strelka image for use in Workflows"

RUN apt-get update -qq && \
	apt-get install -y -qq \
		bzip2 \
		gcc \
		g++ \
		make \
		python \
		wget \
		zlib1g-dev

ENV STRELKA_INSTALL_DIR /opt/strelka/

RUN wget https://github.com/Illumina/strelka/releases/download/v${STRELKA_VERSION}/strelka-${STRELKA_VERSION}.release_src.tar.bz2 \
  && tar -xjf strelka-${STRELKA_VERSION}.release_src.tar.bz2 \
  && rm -rf strelka-${STRELKA_VERSION}.release_src.tar.bz2 \
  && mkdir -p ${STRELKA_INSTALL_DIR}/build \
  && mv strelka-${STRELKA_VERSION}.release_src ${STRELKA_INSTALL_DIR} \
  && cd ${STRELKA_INSTALL_DIR} \
  && strelka-${STRELKA_VERSION}.release_src/configure --jobs=4 --prefix=${STRELKA_INSTALL_DIR} \
  && make install

ENV PATH="${STRELKA_INSTALL_DIR}/bin:${PATH}"