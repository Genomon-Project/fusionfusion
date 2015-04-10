#! /usr/local/bin/python

import sys, re, myCigar

inputFile = sys.argv[1]
cigarSRe_right = re.compile('(\d+)S$')
cigarSRe_left = re.compile('^(\d+)S')

abnormal_insert_size = 500000
min_major_clip_size = 15 
# max_minor_clip_size = 15 

# function for organizing and print junction information
def printJuncInfo(juncLine):

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
        coverRegion_SA = myCigar.getCoverRegion(F[2], F[3], F[5])
        endPos_SA = myCigar.getEndPos(pos_SA, F[5])

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
            coverRegion_primary = myCigar.getCoverRegion(F[2], F[3], F[5])
            readLength_primary = len(F[9])
            endPos_primary = myCigar.getEndPos(pos_primary, F[5])
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
            coverRegion_pair = myCigar.getCoverRegion(F[2], F[3], F[5])
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
                print '\t'.join([juncChr_primary, str(juncPos_primary), juncDir_primary, juncChr_SA, str(juncPos_SA), juncDir_SA, \
                                 juncSurplus, readID_primary, mq_primary, coverRegion_primary, dir_primary, \
                                 mq_pair, coverRegion_pair, dir_pair, mq_SA, coverRegion_SA, dir_SA, str(juncType) , "1"])

            else:
                print '\t'.join([juncChr_SA, str(juncPos_SA), juncDir_SA, juncChr_primary, str(juncPos_primary), juncDir_primary, \
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
                print '\t'.join([juncChr_primary, str(juncPos_primary), juncDir_primary, juncChr_SA, str(juncPos_SA), juncDir_SA, \
                                 juncSurplus, readID_primary, mq_primary, coverRegion_primary, dir_primary, \
                                 mq_pair, coverRegion_pair, dir_pair, mq_SA, coverRegion_SA, dir_SA, str(juncType) , "1"])
                                 
            else:                
                print '\t'.join([juncChr_SA, str(juncPos_SA), juncDir_SA, juncChr_primary, str(juncPos_primary), juncDir_primary, \
                                 juncSurplus, readID_primary, mq_primary, coverRegion_primary, dir_primary, \
                                 mq_pair, coverRegion_pair, dir_pair, mq_SA, coverRegion_SA, dir_SA, str(juncType) , "2"])


hIN = open(inputFile, 'r')

tempID = ""
tempLine = []
for line in hIN:
    line = line.rstrip('\n')
    F = line.split('\t')
    if F[0][0] == "@": continue

    if tempID != F[0]:
        if tempID != "" and len(tempLine) == 3:
            printJuncInfo(tempLine)

        tempID = F[0]
        tempLine = []

    tempLine.append(line)

hIN.close()


if tempID != "" and len(tempLine) == 3:
    printJuncInfo(tempLine)

