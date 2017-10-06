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
    python-pip

RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 && \
    tar jxvf htslib-1.3.2.tar.bz2 && \
    cd htslib-1.3.2 && \
    make && \
    make install

RUN pip install --upgrade pip
RUN pip install --upgrade setuptools

RUN wget https://github.com/friend1ws/annot_utils/archive/v0.1.0.tar.gz && \
    tar xzvf v0.1.0.tar.gz && \
    cd annot_utils-0.1.0/resource && \
    bash prep_data.sh && \
    cd .. && \
    python setup.py build install



RUN apt-get update && apt-get install -y \
    libkrb5-3 \
    libpng12-0

RUN cd  /usr/local/bin && \
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
RUN chmod a+x /usr/local/bin/blat

RUN pip install --upgrade pip
RUN pip install pysam

RUN wget https://github.com/Genomon-Project/fusionfusion/archive/v0.3.0.tar.gz && \
    tar xzvf v0.3.0.tar.gz && \
    cd fusionfusion-0.3.0 && \
    python setup.py build install

