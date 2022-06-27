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
    python2 \
    python2-dev \
    python-pip
    
RUN wget https://github.com/samtools/htslib/releases/download/1.15/htslib-1.15.tar.bz2 && \
    tar jxvf htslib-1.15.tar.bz2 && \
    cd htslib-1.15 && \
    make && \
    make install

#RUN pip install --upgrade pip
#RUN pip2 install --upgrade setuptools

RUN pip2 install pysam==0.19.1
RUN pip2 install annot-utils==0.3.1
RUN pip2 install chimera_utils==0.6.0

RUN cd  /usr/local/bin && \
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
    chmod a+x /usr/local/bin/blat

RUN wget -q https://github.com/friend1ws/fusion_utils/archive/v0.2.0.tar.gz && \
    tar xzvf v0.2.0.tar.gz && \
    cd fusion_utils-0.2.0 && \
    python2 setup.py install
    
RUN wget -q https://github.com/aokad/fusionfusion/archive/refs/tags/v0.5.2b1.tar.gz && \
    tar xzvf v0.5.2b1.tar.gz && \
    cd fusionfusion-0.5.2b1 && \
    python2 setup.py install && \
    echo "success!"
