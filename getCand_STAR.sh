#! /bin/sh

INPUT=$1

REF=/home/yshira/common/ref/hg19_all/hg19.all.fasta
BLAT_ALL_REF=/home/yshira/common/ref/hg19_all/hg19.all.2bit
BLAT_OOC=/home/yshira/common/ref/hg19_all/11.ooc

:<<_COMMENT_OUT_

echo "python getJuncInfo_STAR.py ${INPUT} | sort -k1,1 -k2,2n -k4,4 -k 5,5n - > junction_STAR.txt"
python getJuncInfo_STAR.py ${INPUT} | sort -k1,1 -k2,2n -k4,4 -k 5,5n - > junction_STAR.txt

echo "python summarizeJunc_ms.py junction_STAR.txt > junction_STAR.summarized.txt"
python summarizeJunc_ms.py junction_STAR.txt > junction_STAR.summarized.txt

echo "python filterJunc_ms.py junction_STAR.summarized.txt 3 0.8 50 > junction_STAR.summarized.filt.txt"
python filterJunc_ms.py junction_STAR.summarized.txt 3 0.8 50 > junction_STAR.summarized.filt.txt

_COMMENT_OUT_
echo "python getSplicing.py junction_STAR.summarized.filt.txt ${REF} > junction_STAR.summarized.filt2.txt"
python getSplicing.py junction_STAR.summarized.filt.txt ${REF} > junction_STAR.summarized.filt2.txt 

echo "python makeJuncSeqPairFa.py junction_STAR.summarized.filt2.txt > junction_STAR.summarized.contig.fa"
python makeJuncSeqPairFa.py junction_STAR.summarized.filt2.txt > junction_STAR.summarized.contig.fa 


echo "blat -stepSize=5 -repMatch=2253 -ooc=${BLAT_OOC} ${BLAT_ALL_REF} junction_STAR.summarized.contig.fa junction_STAR.summarized.contig.psl"
blat -stepSize=5 -repMatch=2253 -ooc=${BLAT_OOC} ${BLAT_ALL_REF} junction_STAR.summarized.contig.fa junction_STAR.summarized.contig.psl


