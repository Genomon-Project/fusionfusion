#! /usr/bin/env python

from __future__ import print_function

import sys
import re
import pysam
import collections
# import config
from .config import *
from . import cigar_utils

# for mapsplice2
ReFus_ms2 = re.compile('FUS_(\d+)_(\d+)\(([\-\+])([\-\+])\)')

# for TopHat2
ReFus_th2 = re.compile('XF:Z:(\d) ([^ \t\n\r\f\v,]+)\-([^ \t\n\r\f\v,]+) (\d+) ([\dMmNnIiDd]+[Mm])(\d+F)([\dMmNnIiDd]+[Mm])')

# for STAR
cigarSRe_right = re.compile('(\d+)S$')
cigarSRe_left = re.compile('^(\d+)S')


def extractFusionReads_ms2(inputFilePath, outputFilePath):

    inBamFile = pysam.AlignmentFile(inputFilePath, "r") 

    fusionReadID = {}
    for read in inBamFile.fetch():

        for item in read.tags:
            if item[0] == "ZF" and ReFus_ms2.search(item[1]) is not None:
                fusionReadID[read.qname] = 1

    
    # it seem that reset can be used only for indexed bam files...
    # inBamFile.reset() # 

    inBamFile.close()

    inBamFile = pysam.AlignmentFile(inputFilePath, "r")
    outBamFile = pysam.AlignmentFile(outputFilePath, "w", template = inBamFile)
    for read in inBamFile.fetch():
        
        if read.qname in fusionReadID:
            outBamFile.write(read)


    inBamFile.close()
    outBamFile.close()

    
def extractFusionReads_th2(inputFilePath, outputFilePath):

    inBamFile = pysam.AlignmentFile(inputFilePath, "r")

    fusionReadID = {}
    for read in inBamFile.fetch():

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip if either of the read pair is unmapped
        if flags[2] == "1" or flags[3] == "1": continue

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue

        for item in read.tags:
            if item[0] == "XF":
                fusionReadID[read.qname] = 1


    inBamFile.close()

    inBamFile = pysam.AlignmentFile(inputFilePath, "r")
    outBamFile = pysam.AlignmentFile(outputFilePath, "w", template = inBamFile)
    for read in inBamFile.fetch():

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue


        if read.qname in fusionReadID:
            outBamFile.write(read)


    inBamFile.close()
    outBamFile.close()



def getFusInfo_ms2(tempID, tempLine, fusInfo):

    # abnormal_insert_size = config.param_conf.getint("parse_condition", "abnormal_insert_size")
    abnormal_insert_size = param_conf.abnormal_insert_size

    # check the fusion validity
    ufusInfo = list(set(fusInfo))
    for fus in ufusInfo:
        if fus == 0: continue

        chr_primary, pos_primary, dir_primary, chr_pair, pos_pair, dir_pair, chr_chimera, pos_chimera, dir_chimera = "*", "*", "*", "*", "*", "*", "*", "*", "*"
        mq_primary, cover_primary, mq_pair, cover_pair, mq_chimera, cover_chimera = "*", "*", "*", "*", "*", "*"
        # index showing the order of fusion 
        # 1: FUS_(breakPos_primary)_(breakPos_chimera)_(breakDir_primary)_(breakDir_chimera),
        # -1: FUS_(breakPos_chimera)_(breakPos_primary)_(breakDir_chimera)_(breakDir_primary)
        fusOrder = 0
        
        for i in range(0, len(tempLine)):
             
            FF = tempLine[i].split('\t')
            flags = format(int(FF[1]), "#014b")[:1:-1]
            if flags[8] == "1" and fusInfo[i] == fus:
                chr_chimera = FF[2]
                pos_chimera = FF[3]
                dir_chimera = ("+" if flags[4] != "1" else "-")
                mq_chimera = FF[4]
                cover_chimera = cigar_utils.getCoverRegion(FF[2], FF[3], FF[5])
                if fusOrder == 0: fusOrder = -1
            elif flags[8] != "1" and fusInfo[i] == fus:
                chr_primary = FF[2]
                pos_primary = FF[3]
                dir_primary = ("+" if flags[4] != "1" else "-")
                mq_primary = FF[4]
                cover_primary = cigar_utils.getCoverRegion(FF[2], FF[3], FF[5])
                if fusOrder == 0: fusOrder = 1
            elif flags[8] != "1" and flags[2] != "1" and (fusInfo[i] != fus or chr_primary != "*"):
                chr_pair = FF[2]
                pos_pair = FF[3]
                dir_pair = ("+" if flags[4] != "1" else "-")
                mq_pair = FF[4]
                cover_pair = cigar_utils.getCoverRegion(FF[2], FF[3], FF[5])

        if chr_primary == "*": continue
        fusSplit = fus.split(',')
        if fusOrder == 1:
            breakPos_primary, breakPos_chimera, breakDir_primary, breakDir_chimera = fusSplit[0], fusSplit[1], fusSplit[2], fusSplit[3]
            breakDir_chimera = ("+" if breakDir_chimera == "-" else "-") 
        else:
            breakPos_primary, breakPos_chimera, breakDir_primary, breakDir_chimera = fusSplit[1], fusSplit[0], fusSplit[3], fusSplit[2]
            breakDir_primary = ("+" if breakDir_primary == "-" else "-")

        pairPos = 0
        if pos_pair != "*":
            if breakDir_primary == "+" and dir_pair == "+" and int(breakPos_primary) - abnormal_insert_size <= int(pos_pair) <= int(breakPos_primary): pairPos = 1
            if breakDir_primary == "-" and dir_pair == "-" and int(breakPos_primary) <= int(pos_pair) <= int(breakPos_primary) + abnormal_insert_size: pairPos = 1
            if breakDir_chimera == "+" and dir_pair == "+" and int(breakPos_chimera) - abnormal_insert_size <= int(pos_pair) <= int(breakPos_chimera): pairPos = 2
            if breakDir_chimera == "-" and dir_pair == "-" and int(breakPos_chimera) <= int(pos_pair) <= int(breakPos_chimera) + abnormal_insert_size: pairPos = 2

        if chr_primary < chr_chimera or chr_primary == chr_chimera and breakPos_primary <= breakPos_chimera:
            return '\t'.join([chr_primary, breakPos_primary, breakDir_primary, chr_chimera, breakPos_chimera, breakDir_chimera, "---", tempID, \
                             mq_primary, cover_primary, dir_primary, mq_pair, cover_pair, dir_pair, \
                             mq_chimera, cover_chimera, dir_chimera, str(pairPos), "1"] )
        else:
            return '\t'.join([chr_chimera, breakPos_chimera, breakDir_chimera, chr_primary, breakPos_primary, breakDir_primary, "---", tempID, \
                             mq_primary, cover_primary, dir_primary, mq_pair, cover_pair, dir_pair, \
                             mq_chimera, cover_chimera, dir_chimera, str(pairPos), "2"] )



def getFusInfo_th2(tempID, tempLine, fusInfo, SAFlag):

    # abnormal_insert_size = config.param_conf.getint("parse_condition", "abnormal_insert_size")
    abnormal_insert_size = param_conf.abnormal_insert_size

    # check the fusion validity
    ufusInfo = list(set(fusInfo))
    for fus in ufusInfo:
        if fus == '*': continue

        chr_primary, pos_primary, dir_primary, chr_pair, pos_pair, dir_pair, chr_chimera, pos_chimera, dir_chimera = "*", "*", "*", "*", "*", "*", "*", "*", "*"
        mq_primary, cover_primary, mq_pair, cover_pair, mq_chimera, cover_chimera = "*", "*", "*", "*", "*", "*"

        for i in range(0, len(tempLine)):

            FF = tempLine[i].split('\t')
            if FF[0] == "HWI-ST1165:73:AD0JHGACXX:7:1102:18080:166569":
                pass

            flags = format(int(FF[1]), "#014b")[:1:-1]
            if fusInfo[i] == fus:
                if str(SAFlag[i]) == "1":
                    if chr_primary == '*':
                        chr_primary = FF[2]
                        pos_primary = FF[3]
                        dir_primary = ("+" if flags[4] != "1" else "-")
                        mq_primary = FF[4]
                        cover_primary = cigar_utils.getCoverRegion(FF[2], FF[3], FF[5])
                    else:
                        chr_pair = FF[2]
                        pos_pair = FF[3]
                        dir_pair = ("+" if flags[4] != "1" else "-")
                        mq_pair = FF[4]
                        cover_pair = cigar_utils.getCoverRegion(FF[2], FF[3], FF[5])
                else:
                    chr_chimera = FF[2]
                    pos_chimera = FF[3]
                    dir_chimera = ("+" if flags[4] != "1" else "-")
                    mq_chimera = FF[4]
                    cover_chimera = cigar_utils.getCoverRegion(FF[2], FF[3], FF[5])
            else:
                chr_pair = FF[2]
                pos_pair = FF[3]
                dir_pair = ("+" if flags[4] != "1" else "-")
                mq_pair = FF[4]
                cover_pair = cigar_utils.getCoverRegion(FF[2], FF[3], FF[5])


        if chr_primary == "*": continue
        fusSplit = fus.split(',')
        cigar_primary = fusSplit[3]
        cigar_chimera = fusSplit[5] 
        
        breakDir_primary = ("-" if cigar_primary.islower() else "+")
        breakDir_chimera = ("+" if cigar_chimera.islower() else "-")
        breakPos_primary = str(pos_primary if cigar_primary.islower() else cigar_utils.getEndPos(pos_primary, cigar_primary.upper()))
        breakPos_chimera = str(cigar_utils.getEndPos(pos_chimera, cigar_chimera.upper()) if cigar_chimera.islower() else pos_chimera)

        pairPos = 0
        if pos_pair != "*":
            if chr_primary == chr_pair and breakDir_primary == "+" and dir_pair == "+" and int(breakPos_primary) - abnormal_insert_size <= int(pos_pair) <= int(breakPos_primary): pairPos = 1
            if chr_primary == chr_pair and breakDir_primary == "-" and dir_pair == "-" and int(breakPos_primary) <= int(pos_pair) <= int(breakPos_primary) + abnormal_insert_size: pairPos = 1
            if chr_chimera == chr_pair and breakDir_chimera == "+" and dir_pair == "+" and int(breakPos_chimera) - abnormal_insert_size <= int(pos_pair) <= int(breakPos_chimera): pairPos = 2
            if chr_chimera == chr_pair and breakDir_chimera == "-" and dir_pair == "-" and int(breakPos_chimera) <= int(pos_pair) <= int(breakPos_chimera) + abnormal_insert_size: pairPos = 2

        if pairPos == 0: continue

        if chr_primary < chr_chimera or chr_primary == chr_chimera and breakPos_primary <= breakPos_chimera:
            return '\t'.join([chr_primary, breakPos_primary, breakDir_primary, chr_chimera, breakPos_chimera, breakDir_chimera, "---", tempID, \
                             mq_primary, cover_primary, dir_primary, mq_pair, cover_pair, dir_pair, \
                             mq_chimera, cover_chimera, dir_chimera, str(pairPos), "1"] )
        else:
            return '\t'.join([chr_chimera, breakPos_chimera, breakDir_chimera, chr_primary, breakPos_primary, breakDir_primary, "---", tempID, \
                             mq_primary, cover_primary, dir_primary, mq_pair, cover_pair, dir_pair, \
                             mq_chimera, cover_chimera, dir_chimera, str(pairPos), "2"] )


def parseJuncInfo_ms2(inputFilePath, outputFilePath):

    """
    script for collecting short reads supporting fusion candidates in MapSplice2 sam file
    """

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    tempID = ""
    fusFlag = [] 
    fusInfo = []
    tempLine = []

    for line in hIN:
        if line[0] == "@": continue
        line = line.rstrip('\n')
        F = line.split('\t')

        if tempID != F[0]:
            if fusInfo.count(0) != len(tempLine):
                tempFusInfo = getFusInfo_ms2(tempID, tempLine, fusInfo)
                if tempFusInfo is not None:
                    print(tempFusInfo, file = hOUT)

            tempID = F[0]
            fusFlag = []
            fusInfo = [] 
            tempLine = [] 
    
        tempLine.append(line)
        mFus = ReFus_ms2.search('\t'.join(F[11:]))
        if mFus is not None:
            fusFlag.append(1)
            fusInfo.append(','.join([mFus.group(1), mFus.group(2), mFus.group(3), mFus.group(4)]))
        else:
            fusFlag.append(0)
            fusInfo.append(0)

    hIN.close()


    if fusInfo.count(0) != len(tempLine):
        tempFusInfo = getFusInfo_ms2(tempID, tempLine, fusInfo)
        if tempFusInfo is not None:
            print(tempFusInfo, file = hOUT)

    hOUT.close()


def parseJuncInfo_th2(inputFilePath, outputFilePath):

    """
    script for collecting short reads supporting fusion candidates in TopHat2 sam file
    """

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    tempID = ""
    fusInfo = []
    SAFlag = []
    tempLine = []

    for line in hIN:
        if line[0] == "@": continue
        line = line.rstrip('\n')
        F = line.split('\t')

        if tempID != F[0]:
            if fusInfo.count('*') != len(tempLine):
                tempFusInfo = getFusInfo_th2(tempID, tempLine, fusInfo, SAFlag)
                if tempFusInfo is not None:
                    print(tempFusInfo, file = hOUT)

            tempID = F[0]
            fusInfo = []
            SAFlag = []
            tempLine = []

        tempLine.append(line)
        mFus = ReFus_th2.search('\t'.join(F[11:]))
        if mFus is not None:
            fusInfo.append(','.join([mFus.group(2), mFus.group(3), mFus.group(4), mFus.group(5), mFus.group(6), mFus.group(7)]))
            SAFlag.append(mFus.group(1))
        else:
            fusInfo.append('*')
            SAFlag.append('*')

    hIN.close()


    if fusInfo.count('*') != len(tempLine):
        tempFusInfo = getFusInfo_th2(tempID, tempLine, fusInfo, SAFlag)
        if tempFusInfo is not None:
            print(tempFusInfo, file = hOUT)


    hOUT.close()


def getFusInfo_STAR(juncLine, source=None):

    # abnormal_insert_size = config.param_conf.getint("parse_condition", "abnormal_insert_size")
    # min_major_clip_size = config.param_conf.getint("parse_condition", "min_major_clip_size")
    abnormal_insert_size = param_conf.abnormal_insert_size
    min_major_clip_size = param_conf.min_major_clipping_size
 
    """
    function for organizing and print junction information
    """

    # be reminded necessary variables
    chr_primary, pos_primary, dir_primary, chr_pair, pos_pair, dir_pair, chr_SA, pos_SA, dir_SA = "", "", "", "", "", "", "", "", ""
    mq_primary, coverRegion_primary, mq_pair, coverRegion_pair, mq_SA, coverRegion_SA = "", "", "", "", "", ""
    juncChr_primary, juncPos_primary, juncDir_primary, juncChr_SA, juncPos_SA, juncDir_SA = "", "", "", "", "", "" 

    # collect about the supplementary alingment information
    right_clipping_SA = 0
    left_clipping_SA = 0
    flags_SA = ""
    for line in juncLine:
        F = line.split('\t')
        flags = format(int(F[1]), "#014b")[:1:-1] 
        if flags[8] != "1": continue

        if chr_SA != "":
            print("Multiple supplementary alignment at:" + '\n' +  '\n'.join(juncLine), file = sys.stderr)

        chr_SA = F[2]
        juncChr_SA = F[2]
        pos_SA = int(F[3])
        dir_SA = ("-" if flags[4] == "1" else "+")
        mq_SA = F[4]
        coverRegion_SA = cigar_utils.getCoverRegion(F[2], F[3], F[5])
        endPos_SA = cigar_utils.getEndPos(pos_SA, F[5])

        flags_SA = flags
        if flags_SA[6] == flags_SA[7]: print("The supplementary Read is both first and second reads at:" + '\n' + '\n'.join(juncLine), file = sys.stderr)

        tempMatch = cigarSRe_right.search(F[5])
        if tempMatch is not None: right_clipping_SA = int(tempMatch.group(1))

        tempMatch = cigarSRe_left.search(F[5])
        if tempMatch is not None: left_clipping_SA = int(tempMatch.group(1))

    # collect information about primary junction read and its pair read
    right_clipping_primary = 0
    left_clipping_primary = 0
    cigar_primary = ""
    readLength_primary = 0
    readID_primary = ""
    for line in juncLine:
        F = line.split('\t')
        flags = format(int(F[1]), "#014b")[:1:-1]
        if flags[8] == "1": continue

        # primary junction read
        if flags[6] == flags_SA[6] and flags[7] == flags_SA[7]:
            chr_primary = F[2]
            juncChr_primary = F[2]
            pos_primary = int(F[3])
            dir_primary = ("-" if flags[4] == "1" else "+")
            mq_primary = F[4]
            coverRegion_primary = cigar_utils.getCoverRegion(F[2], F[3], F[5])
            readLength_primary = len(F[9])
            endPos_primary = cigar_utils.getEndPos(pos_primary, F[5])
            readID_primary = "{qname}/{read}{suffix}".format(
                qname=F[0],
                read=(1 if flags[6] == "1" else 2),
                suffix=("@" + source if source else "")
            )

            tempMatch = cigarSRe_right.search(F[5])
            if tempMatch is not None: right_clipping_primary = int(tempMatch.group(1))
        
            tempMatch = cigarSRe_left.search(F[5])
            if tempMatch is not None: left_clipping_primary = int(tempMatch.group(1))

        elif flags[6] == flags_SA[7] and flags[7] == flags_SA[6]:
            chr_pair = F[2]
            juncChr_pair = F[2]
            pos_pair = int(F[3])
            dir_pair = ("-" if flags[4] == "1" else "+")
            mq_pair = F[4]
            coverRegion_pair = cigar_utils.getCoverRegion(F[2], F[3], F[5])
        else:
            print("The following read is both first and second reads at:" + '\n' + line, file = sys.stderr)


    if right_clipping_primary >= min_major_clip_size:

        clipLen_primary = right_clipping_primary

        juncChr_primary = chr_primary
        juncPos_primary = endPos_primary
        juncDir_primary = "+"

        expected_clipLen_SA = readLength_primary - clipLen_primary
        expected_clipDir_SA = ("-" if dir_primary== dir_SA else "+")

        validFlag = 0
        juncDir_SA = ""
        juncPos_SA = ""

        # the pair read is aligned at the same chromosome with the primary read
        if dir_primary == "-" and dir_pair == "+" and chr_primary == chr_pair and 0 <= pos_primary - pos_pair < abnormal_insert_size:

            if dir_SA == "+" and expected_clipDir_SA == "+" and right_clipping_SA > 0:
                clipLen_SA = right_clipping_SA
                juncDir_SA = "+"
                juncPos_SA = endPos_SA
                if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA - (expected_clipLen_SA - clipLen_SA)
                juncType = 1
                validFlag = 1

            if dir_SA == "-" and expected_clipDir_SA == "-" and left_clipping_SA > 0:
                clipLen_SA = left_clipping_SA
                juncDir_SA = "-"
                juncPos_SA = int(pos_SA)
                if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA + (expected_clipLen_SA - clipLen_SA)
                juncType = 1
                validFlag = 1

        # when the supplementary read is aligned on the same chromosome with the paired read
        if dir_primary == "+" and chr_SA == chr_pair:

            if dir_SA == "-" and dir_pair == "+" and 0 <= pos_SA - pos_pair < abnormal_insert_size and expected_clipDir_SA == "+" and right_clipping_SA > 0:
                clipLen_SA = right_clipping_SA
                juncDir_SA = "+"
                juncPos_SA = endPos_SA
                if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA - (expected_clipLen_SA - clipLen_SA)
                juncType = 2
                validFlag = 1

            if dir_SA == "+" and dir_pair == "-" and 0 <= pos_pair - pos_SA < abnormal_insert_size and expected_clipDir_SA == "-" and left_clipping_SA > 0:
                clipLen_SA = left_clipping_SA
                juncDir_SA = "-"
                juncPos_SA = pos_SA
                if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA + (expected_clipLen_SA - clipLen_SA)
                juncType = 2
                validFlag = 1

        if validFlag == 1:

            juncSurplus = "---"
            if clipLen_SA > expected_clipLen_SA:
                surPlus_start = readLength_primary - clipLen_primary
                surPlus_end = surPlus_start + clipLen_SA - expected_clipLen_SA
                juncSurplus = read.seq[surPlus_start:surPlus_end]

            # reorder by the chromosome position and print
            if juncChr_primary < juncChr_SA or juncChr_primary == juncChr_SA and juncPos_primary <= juncPos_SA:
                return '\t'.join([juncChr_primary, str(juncPos_primary), juncDir_primary, juncChr_SA, str(juncPos_SA), juncDir_SA, \
                                 juncSurplus, readID_primary, mq_primary, coverRegion_primary, dir_primary, \
                                 mq_pair, coverRegion_pair, dir_pair, mq_SA, coverRegion_SA, dir_SA, str(juncType) , "1"])

            else:
                return '\t'.join([juncChr_SA, str(juncPos_SA), juncDir_SA, juncChr_primary, str(juncPos_primary), juncDir_primary, \
                                 juncSurplus, readID_primary, mq_primary, coverRegion_primary, dir_primary, \
                                 mq_pair, coverRegion_pair, dir_pair, mq_SA, coverRegion_SA, dir_SA, str(juncType) , "2"])


    if left_clipping_primary >= min_major_clip_size:

        clipLen_primary = left_clipping_primary
        juncChr_primary = chr_primary
        juncPos_primary = int(pos_primary)
        juncDir_primary = "-"
        juncChr_SA = chr_SA

        expected_clipLen_SA = readLength_primary - clipLen_primary
        expected_clipDir_SA = ("+" if dir_primary== dir_SA else "-")

        validFlag = 0
        juncDir_SA = ""
        juncPos_SA = ""
        # the pair read is aligned at the same chromosome with the primary read
        if dir_primary == "+" and dir_pair == "-" and chr_primary == chr_pair and 0 <= pos_pair - pos_primary < abnormal_insert_size:

            if dir_SA == "+" and expected_clipDir_SA == "+" and right_clipping_SA > 0:
                clipLen_SA = right_clipping_SA
                juncDir_SA = "+"
                juncPos_SA = endPos_SA
                if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA - (expected_clipLen_SA - clipLen_SA)
                juncType = 1
                validFlag = 1

            if dir_SA == "-" and expected_clipDir_SA == "-" and left_clipping_SA > 0:
                clipLen_SA = left_clipping_SA
                juncDir_SA = "-"
                juncPos_SA = pos_SA
                if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA + (expected_clipLen_SA - clipLen_SA)
                juncType = 1
                validFlag = 1

        if dir_primary == "-" and chr_SA == chr_pair:

            if dir_SA == "-" and dir_pair == "+" and 0 <= pos_SA - pos_pair < abnormal_insert_size and expected_clipDir_SA == "+" and right_clipping_SA > 0:
                clipLen_SA = right_clipping_SA
                juncDir_SA = "+"
                juncPos_SA = endPos_SA 
                if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA - (expected_clipLen_SA - clipLen_SA)
                juncType = 2
                validFlag = 1

            if dir_SA == "+" and dir_pair == "-" and 0 <= pos_pair - pos_SA < abnormal_insert_size and expected_clipDir_SA == "-" and left_clipping_SA:
                clipLen_SA = left_clipping_SA
                juncDir_SA = "-"
                juncPos_SA = pos_SA
                if clipLen_SA < expected_clipLen_SA: juncPos_SA = juncPos_SA + (expected_clipLen_SA - clipLen_SA)
                juncType = 2
                validFlag = 1


        if validFlag == 1:

            juncSurplus = "---"
            if clipLen_SA > expected_clipLen_SA:
                surPlus_end = clipLen_primary # this is right
                surPlus_start = surPlus_end - (clipLen_SA - expected_clipLen_SA)
                juncSurplus = read.seq[surPlus_start:surPlus_end]

            # reorder by the chromosome position and print
            if juncChr_primary < juncChr_SA or juncChr_primary == juncChr_SA and juncPos_primary <= juncPos_SA:
                return '\t'.join([juncChr_primary, str(juncPos_primary), juncDir_primary, juncChr_SA, str(juncPos_SA), juncDir_SA, \
                                 juncSurplus, readID_primary, mq_primary, coverRegion_primary, dir_primary, \
                                 mq_pair, coverRegion_pair, dir_pair, mq_SA, coverRegion_SA, dir_SA, str(juncType) , "1"])
                                 
            else:                
                return '\t'.join([juncChr_SA, str(juncPos_SA), juncDir_SA, juncChr_primary, str(juncPos_primary), juncDir_primary, \
                                 juncSurplus, readID_primary, mq_primary, coverRegion_primary, dir_primary, \
                                 mq_pair, coverRegion_pair, dir_pair, mq_SA, coverRegion_SA, dir_SA, str(juncType) , "2"])



def parseJuncInfo_STAR(inputFilePath, outputFilePath, source=None):

     
    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    tempID = ""
    tempLine = []
    for line in hIN:
        line = line.rstrip('\n')
        F = line.split('\t')
        if F[0][0] == "@": continue

        if tempID != F[0]:
            if tempID != "" and len(tempLine) == 3:
                tempFusInfo = getFusInfo_STAR(tempLine, source)
                if tempFusInfo is not None:
                    print(tempFusInfo, file = hOUT)

            tempID = F[0]
            tempLine = []

        tempLine.append(line)

    hIN.close()


    if tempID != "" and len(tempLine) == 3:
        tempFusInfo = getFusInfo_STAR(tempLine, source)
        if tempFusInfo is not None:
            print(tempFusInfo, file = hOUT)


    hOUT.close()



def clusterJuncInfo(inputFilePath, outputFilePath):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    # I believe this is appropriate..
    check_margin_size = 30

    cluster_junc = {}
    cluster_inseq, cluster_id, cluster_pair_pos, cluster_primary_pos = {}, {}, {}, {}
    cluster_MQ_primary, cluser_cover_primary, cluster_dir_primary = {}, {}, {}
    cluster_MQ_pair, cluster_cover_pair, cluster_dir_pair = {}, {}, {}
    cluster_MQ_SA, cluster_cover_SA, cluster_dir_SA = {}, {}, {}

    for line in hIN:
 
        F = line.rstrip('\n').split('\t')
        match = 0
        delList = []

        for key in sorted(cluster_junc):
            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, inseq_size = key.split('\t')

            # the investigated key is sufficiently far from the current line in the input file and no additional line to merge is expected. therefore flush the key and information
            if F[0] != tchr1 or int(F[1]) > int(tpos1) + check_margin_size:

                # obtain the most frequent junction
                junc_counter = collections.Counter(cluster_junc[key])
                best_junc = junc_counter.most_common(1)[0][0]
                btchr1, btpos1, btdir1, btchr2, btpos2, btdir2, btinseq = best_junc.split('\t') 

                print('\t'.join([btchr1, btpos1, btdir1, btchr2, btpos2, btdir2, btinseq, \
                      ';'.join(cluster_id[key]), ';'.join(cluster_MQ_primary[key]), ';'.join(cluser_cover_primary[key]), ';'.join(cluster_dir_primary[key]), \
                      ';'.join(cluster_MQ_pair[key]), ';'.join(cluster_cover_pair[key]), ';'.join(cluster_dir_pair[key]), \
                      ';'.join(cluster_MQ_SA[key]), ';'.join(cluster_cover_SA[key]), ';'.join(cluster_dir_SA[key]), \
                      ';'.join(cluster_pair_pos[key]), ';'.join(cluster_primary_pos[key])]), file = hOUT)


                # add to the deletion list (later the key will be removed from the dictionaries)
                delList.append(key)
                continue

            else:

                # check whether the investigated key and the current line should be merged or not
                if F[0] == tchr1 and F[3] == tchr2 and F[2] == tdir1 and F[5] == tdir2:

                    flag = 0
                    temp_seq_size = 0 if F[6] == "---" else len(F[6])
                    # detailed check on the junction position considering inserted sequences
                    if F[2] == "+":
                        expected_diff_size = (int(F[1]) - int(tpos1)) + (temp_seq_size - int(inseq_size))
                        if (F[5] == "+" and int(F[4]) == int(tpos2) - expected_diff_size) or (F[5] == "-" and int(F[4]) == int(tpos2) + expected_diff_size):
                            flag = 1
                    else:
                        expected_diff_size = (int(F[1]) - int(tpos1)) + (int(inseq_size) - temp_seq_size)
                        if (F[5] == "+" and int(F[4]) == int(tpos2) + expected_diff_size) or (F[5] == "-" and int(F[4]) == int(tpos2) - expected_diff_size):
                            flag = 1

                    if flag == 1:

                        match = 1
                        cluster_inseq[key].append(F[6])
                        cluster_id[key].append(F[7])
                        cluster_MQ_primary[key].append(F[8])
                        cluser_cover_primary[key].append(F[9])
                        cluster_dir_primary[key].append(F[10])
                        cluster_MQ_pair[key].append(F[11])
                        cluster_cover_pair[key].append(F[12])
                        cluster_dir_pair[key].append(F[13])
                        cluster_MQ_SA[key].append(F[14])
                        cluster_cover_SA[key].append(F[15])
                        cluster_dir_SA[key].append(F[16])
                        cluster_pair_pos[key].append(F[17])
                        cluster_primary_pos[key].append(F[18])

                        # whether to check whether the inserted sequence should be reverse-complemented?
                        cluster_junc[key].append('\t'.join(F[0:7]))

        for key in delList:
            del cluster_junc[key]
            del cluster_inseq[key]
            del cluster_id[key]
            del cluster_MQ_primary[key]
            del cluser_cover_primary[key]
            del cluster_dir_primary[key]
            del cluster_MQ_pair[key]
            del cluster_cover_pair[key]
            del cluster_dir_pair[key]
            del cluster_MQ_SA[key]
            del cluster_cover_SA[key]
            del cluster_dir_SA[key]
            del cluster_pair_pos[key]
            del cluster_primary_pos[key]

        # if the current line in the input file does not match any of the pooled keys
        if match == 0:
            temp_seq_size = 0 if F[6] == "---" else len(F[6])
            key = '\t'.join(F[0:6] + [str(temp_seq_size)])
            # whether to check whether the inserted sequence should be reverse-complemented?
            cluster_junc[key] = ['\t'.join(F[0:7])]

            cluster_inseq[key], cluster_id[key] = [F[6]], [F[7]]
            cluster_MQ_primary[key], cluser_cover_primary[key], cluster_dir_primary[key] = [F[8]], [F[9]], [F[10]]
            cluster_MQ_pair[key], cluster_cover_pair[key], cluster_dir_pair[key] = [F[11]], [F[12]], [F[13]]
            cluster_MQ_SA[key], cluster_cover_SA[key], cluster_dir_SA[key] = [F[14]], [F[15]], [F[16]]
            cluster_pair_pos[key], cluster_primary_pos[key] = [F[17]], [F[18]]

    hIN.close()

    for key in sorted(cluster_junc):

        # obtain the most frequent junction
        junc_counter = collections.Counter(cluster_junc[key])
        best_junc = junc_counter.most_common(1)[0][0]
        btchr1, btpos1, btdir1, btchr2, btpos2, btdir2, btinseq = best_junc.split('\t')

        print('\t'.join([btchr1, btpos1, btdir1, btchr2, btpos2, btdir2, btinseq, \
              ';'.join(cluster_id[key]), ';'.join(cluster_MQ_primary[key]), ';'.join(cluser_cover_primary[key]), ';'.join(cluster_dir_primary[key]), \
              ';'.join(cluster_MQ_pair[key]), ';'.join(cluster_cover_pair[key]), ';'.join(cluster_dir_pair[key]), \
              ';'.join(cluster_MQ_SA[key]), ';'.join(cluster_cover_SA[key]), ';'.join(cluster_dir_SA[key]), \
              ';'.join(cluster_pair_pos[key]), ';'.join(cluster_primary_pos[key])]), file = hOUT)


    hOUT.close()


