FROM ubuntu:16.04
MAINTAINER Yuichi Shiraishi <friend1ws@gmail.com> 


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
    liblzma-dev \
    bwa

RUN wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 && \
    tar jxvf htslib-1.7.tar.bz2 && \
    cd htslib-1.7 && \
    make && \
    make install

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

