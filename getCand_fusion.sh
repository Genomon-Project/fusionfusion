#! /bin/sh

INPUT=$1
OUTPUTDIR=$2
TOOL=$3

REF=/home/yshira/common/ref/hg19_all/hg19.all.fasta
BLAT_ALL_REF=/home/yshira/common/ref/hg19_all/hg19.all.2bit
BLAT_OOC=/home/yshira/common/ref/hg19_all/11.ooc

if [ ! -d ${OUTPUTDIR} ]
then
    mkdir -p ${OUTPUTDIR}
fi

# :<<_COMMENT_OUT_

if [ ${TOOL} = "STAR" ]
then
    echo "python getJuncInfo_STAR.py ${INPUT} | sort -k1,1 -k2,2n -k4,4 -k 5,5n - > ${OUTPUTDIR}/junction.txt" 
    python getJuncInfo_STAR.py ${INPUT} | sort -k1,1 -k2,2n -k4,4 -k 5,5n - > ${OUTPUTDIR}/junction.txt 
elif [ ${TOOL} = "ms2" ]
then
    echo "python getJuncInfo_ms2.py ${INPUT} | sort -k1,1 -k2,2n -k4,4 -k 5,5n - > ${OUTPUTDIR}/junction.txt" 
    python getJuncInfo_ms2.py ${INPUT} | sort -k1,1 -k2,2n -k4,4 -k 5,5n - > ${OUTPUTDIR}/junction.txt 
else
    echo "The specified alignment tool should be STAR or ms2"
    exit
fi


echo "python summarizeJunc.py ${OUTPUTDIR}/junction.txt > ${OUTPUTDIR}/junction.summarized.txt"
python summarizeJunc.py ${OUTPUTDIR}/junction.txt > ${OUTPUTDIR}/junction.summarized.txt 

echo "python filterJunc.py ${OUTPUTDIR}/junction.summarized.txt 3 0.8 50 > ${OUTPUTDIR}/junction.summarized.filt1.txt"
python filterJunc.py ${OUTPUTDIR}/junction.summarized.txt 3 0.8 50 > ${OUTPUTDIR}/junction.summarized.filt1.txt 

echo "python getSplicing.py ${OUTPUTDIR}/junction.summarized.filt1.txt ${REF} > ${OUTPUTDIR}/junction.summarized.filt2.txt"
python getSplicing.py ${OUTPUTDIR}/junction.summarized.filt1.txt ${REF} > ${OUTPUTDIR}/junction.summarized.filt2.txt 

echo "python makeJuncSeqPairFa.py ${OUTPUTDIR}/junction.summarized.filt2.txt > ${OUTPUTDIR}/junction.summarized.contig.fa"
python makeJuncSeqPairFa.py ${OUTPUTDIR}/junction.summarized.filt2.txt > ${OUTPUTDIR}/junction.summarized.contig.fa 

echo "blat -stepSize=5 -repMatch=2253 -ooc=${BLAT_OOC} ${BLAT_ALL_REF} ${OUTPUTDIR}/junction.summarized.contig.fa ${OUTPUTDIR}/junction.summarized.contig.psl"
blat -stepSize=5 -repMatch=2253 -ooc=${BLAT_OOC} ${BLAT_ALL_REF} ${OUTPUTDIR}/junction.summarized.contig.fa ${OUTPUTDIR}/junction.summarized.contig.psl  

# _COMMENT_OUT_

echo "python checkMatching.py ${OUTPUTDIR}/junction.summarized.contig.psl > ${OUTPUTDIR}/junction.summarized.contig.check.txt"
python checkMatching.py ${OUTPUTDIR}/junction.summarized.contig.psl > ${OUTPUTDIR}/junction.summarized.contig.check.txt

echo "python filterMatchCheck.py ${OUTPUTDIR}/junction.summarized.filt2.txt ${OUTPUTDIR}/junction.summarized.contig.check.txt 3 > ${OUTPUTDIR}/junction.summarized.filt3.txt"
python filterMatchCheck.py ${OUTPUTDIR}/junction.summarized.filt2.txt ${OUTPUTDIR}/junction.summarized.contig.check.txt 3 > ${OUTPUTDIR}/junction.summarized.filt3.txt 

