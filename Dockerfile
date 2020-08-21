FROM ubuntu:16.04
MAINTAINER Yuichi Shiraishi <friend1ws@gmail.com> 

ARG MINIMAP2_VERSION=2.17
ARG MINIMAP2_RELEASE=minimap2-${MINIMAP2_VERSION}_x64-linux

RUN apt-get update && apt-get install -y \
    git \
    wget \
    bzip2 \
    make \
    gcc \
    zlib1g-dev \
    python \
    python-pip \
    libbz2-dev \
    liblzma-dev

RUN wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 && \
    tar jxvf htslib-1.7.tar.bz2 && \
    cd htslib-1.7 && \
    make && \
    make install

RUN wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/${MINIMAP2_RELEASE}.tar.bz2 && \
    tar xf ${MINIMAP2_RELEASE}.tar.bz2 && \
    cp ${MINIMAP2_RELEASE}/minimap2 /usr/local/bin && \
    rm -rf ${MINIMAP2_RELEASE}

RUN pip install --upgrade pip
RUN pip install --upgrade setuptools

RUN pip install pysam==0.13
RUN pip install annot-utils==0.2.0
RUN wget https://github.com/Genomon-Project/fusionfusion/archive/v0.4.1.tar.gz && \
    tar -zxvf v0.4.1.tar.gz && \
    cd fusionfusion-0.4.1 && \
    python setup.py build install

RUN apt-get update && apt-get install -y \
    libkrb5-3 \
    libpng12-0

