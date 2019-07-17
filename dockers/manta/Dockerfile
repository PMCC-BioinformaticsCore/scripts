# Based on https://github.com/Illumina/manta/blob/master/docs/userGuide/installation.md


FROM ubuntu:16.04
MAINTAINER Michael Franklin <michael.franklin@petermac.org>

ENV MANTA_VERSION 1.5.1

LABEL \
  version="1.5.1" \
  description="Manta image for use in Workflows"

RUN apt-get update -qq \
  && apt-get install -qq bzip2 gcc g++ make python zlib1g-dev wget

ENV MANTA_INSTALL_DIR /opt/manta/

RUN wget https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.release_src.tar.bz2
RUN tar -xjf manta-${MANTA_VERSION}.release_src.tar.bz2 \
  && rm -rf manta-${MANTA_VERSION}.release_src.tar.bz2 \
  && mkdir -p $MANTA_INSTALL_DIR/build \
  && mv manta-${MANTA_VERSION}.release_src $MANTA_INSTALL_DIR \
  && cd $MANTA_INSTALL_DIR \
  && manta-${MANTA_VERSION}.release_src/configure --jobs=4 --prefix=$MANTA_INSTALL_DIR \
  && make install

ENV PATH="${MANTA_INSTALL_DIR}bin:${PATH}"
