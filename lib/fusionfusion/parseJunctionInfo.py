#! /usr/bin/env python

import sys
import re
import cigar_utils
import pysam

import config

# for mapsplice2
ReFus_ms2 = re.compile('FUS_(\d+)_(\d+)\(([\-\+])([\-\+])\)')

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

    

def getFusInfo_ms2(tempID, tempLine, fusInfo):

    abnormal_insert_size = config.param_conf.getint("parse_condition", "abnormal_insert_size")

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
                    print >> hOUT, tempFusInfo

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
            print >> hOUT, tempFusInfo


    hOUT.close()



def getFusInfo_STAR(juncLine):

    abnormal_insert_size = config.param_conf.getint("parse_condition", "abnormal_insert_size")
    min_major_clip_size = config.param_conf.getint("parse_condition", "min_major_clip_size")
 
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
            print >> sys.stderr, "Multiple supplementary alignment at:" + '\n' +  '\n'.join(juncLine)

        chr_SA = F[2]
        juncChr_SA = F[2]
        pos_SA = int(F[3])
        dir_SA = ("-" if flags[4] == "1" else "+")
        mq_SA = F[4]
        coverRegion_SA = cigar_utils.getCoverRegion(F[2], F[3], F[5])
        endPos_SA = cigar_utils.getEndPos(pos_SA, F[5])

        flags_SA = flags
        if flags_SA[6] == flags_SA[7]: print >> sys.stderr, "The supplementary Read is both first and second reads at:" + '\n' + '\n'.join(juncLine)

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
            readID_primary = F[0] + ("/1" if flags[6] == "1" else "/2")

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
            print >> sys.stderr, "The following read is both first and second reads at:" + '\n' + line


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



def parseJuncInfo_STAR(inputFilePath, outputFilePath):

    abnormal_insert_size = config.param_conf.getint("parse_condition", "abnormal_insert_size")
 
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
                tempFusInfo = getFusInfo_STAR(tempLine)
                if tempFusInfo is not None:
                    print >> hOUT, tempFusInfo

            tempID = F[0]
            tempLine = []

        tempLine.append(line)

    hIN.close()


    if tempID != "" and len(tempLine) == 3:
        tempFusInfo = getFusInfo_STAR(tempLine)
        if tempFusInfo is not None:
            print >> hOUT, tempFusInfo


    hOUT.close()



def clusterJuncInfo(inputFilePath, outputFilePath):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    tempKey = ""
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        key = '\t'.join(F[0:6])

        if tempKey != key:
            if tempKey != "":
                print >> hOUT,'\t'.join([tempKey, ';'.join(tempInseq), ';'.join(tempIDs),  ';'.join(tempMQ_primary), ';'.join(tempCover_primary), ';'.join(tempDir_primary), \
                                ';'.join(tempMQ_pair), ';'.join(tempCover_pair), ';'.join(tempDir_pair), \
                                 ';'.join(tempMQ_SA), ';'.join(tempCover_SA), ';'.join(tempDir_SA), ';'.join(tempPairPos), ';'.join(tempPrimaryPos)])
            

            tempKey = key
            tempIDs = []
            tempInseq = []
            tempMQ_primary = []
            tempCover_primary = []
            tempDir_primary = []
            tempMQ_pair = []
            tempCover_pair = []
            tempDir_pair = []
            tempMQ_SA = []
            tempCover_SA = []
            tempDir_SA = []
            tempPairPos = []
            tempPrimaryPos = []

        tempInseq.append(F[6])
        tempIDs.append(F[7])
        tempMQ_primary.append(F[8])
        tempCover_primary.append(F[9])
        tempDir_primary.append(F[10])
        tempMQ_pair.append(F[11])
        tempCover_pair.append(F[12])
        tempDir_pair.append(F[13])
        tempMQ_SA.append(F[14])
        tempCover_SA.append(F[15])
        tempDir_SA.append(F[16])
        tempPairPos.append(F[17])
        tempPrimaryPos.append(F[18])

    hIN.close()

    if tempKey != "":
        print >> hOUT, '\t'.join([tempKey, ';'.join(tempInseq), ';'.join(tempIDs),  ';'.join(tempMQ_primary), ';'.join(tempCover_primary), ';'.join(tempDir_primary), \
                        ';'.join(tempMQ_pair), ';'.join(tempCover_pair), ';'.join(tempDir_pair), \
                         ';'.join(tempMQ_SA), ';'.join(tempCover_SA), ';'.join(tempDir_SA), ';'.join(tempPairPos), ';'.join(tempPrimaryPos)])

    hOUT.close()


