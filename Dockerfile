FROM friend1ws/annot_utils 
MAINTAINER Yuichi Shiraishi <friend1ws@gmail.com> 

WORKDIR /usr/local/bin

RUN apt-get update && apt-get install -y \
    libkrb5-3 \
    libpng12-0

RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
RUN chmod a+x /usr/local/bin/blat

RUN pip install --upgrade pip
RUN pip install pysam

RUN wget https://github.com/Genomon-Project/fusionfusion/archive/v0.3.0.tar.gz && \
    tar xzvf v0.3.0.tar.gz && \
    cd fusionfusion-0.3.0 && \
    python setup.py build install

WORKDIR /data

ENTRYPOINT ["/usr/local/bin/fusionfusion"]
CMD ["--help"]

