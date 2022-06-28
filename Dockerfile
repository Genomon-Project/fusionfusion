FROM ubuntu:22.04
MAINTAINER Yuichi Shiraishi <friend1ws@gmail.com> 


RUN apt-get update && apt-get install -y \
    git \
    wget \
    bzip2 \
    make \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libkrb5-3 \
    libpng16-16 \
    python3 \
    python3-dev \
    python3-pip
    
RUN wget -q https://github.com/samtools/htslib/releases/download/1.15/htslib-1.15.tar.bz2 && \
    tar jxvf htslib-1.15.tar.bz2 && \
    cd htslib-1.15 && \
    make && \
    make install

RUN pip3 install pysam==0.19.1
RUN pip3 install annot-utils==0.3.1
RUN pip3 install chimera_utils==0.6.0

RUN cd  /usr/local/bin && \
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
    chmod a+x /usr/local/bin/blat

RUN git clone https://github.com/Genomon-Project/fusionfusion.git && \
    cd fusionfusion && \
    python3 setup.py build install
