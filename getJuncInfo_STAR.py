#! /usr/local/bin/python

import sys, re, myCigar

inputFile = sys.argv[1]
cigarSRe_right = re.compile('(\d+)S$')
cigarSRe_left = re.compile('^(\d+)S')

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
        F = line.rstrip('\t')
        flags = format(int(F[1]), "#014b")[:1:-1] 
        if flags[8] != "1": continue

        if chr_SA != "":
            print >> sys.stderr, "Multiple supplementary alignment at:" + '\n' +  '\n'.join(juncLine)

        chr_SA = F[2]
        juncChr_SA = F[2]
        pos_SA = F[3]
        dir_SA = ("-" if flags[4] == "1" else "+")
        mq_SA = F[4]
        coverRegion_SA = myCigar.getCoverRegion(F[2], F[3], F[5])

        flags_SA = flags
        if flags_SA[6] == flags_SA[7]: print >> sys.stderr, "The supplementary Read is both first and second reads at:" + '\n' + '\n'.join(juncLine)

        tempMatch = cigarSRe_right.search(F[5])
        if tempMatch is not None: right_clipping_SA = int(tempMatch.group(1))

        tempMatch = cigarSRe_left.search(F[5])
        if tempMatch is not None: left_clipping_SA = int(tempMatch.group(1))


    # collect information about primary junction read and its pair read
    right_clipping_primary = 0
    left_clipping_primary = 0
    for line in juncLine:
        F = line.rstrip('\t')
        flags = format(int(F[1]), "#014b")[:1:-1]
        if flags[8] == "1": continue

        # primary junction read
        if flags[6] == flags_SA[6] and flags[7] == flags_SA[7]:
            chr_primary = F[2]
            juncChr_primary = F[2]
            pos_primary = F[3]
            dir_primary = ("-" if flags[4] == "1" else "+")
            mq_primary = F[4]
            coverRegion_primary = myCigar.getCoverRegion(F[2], F[3], F[5])

            tempMatch = cigarSRe_right.search(F[5])
            if tempMatch is not None: right_clipping_primary = int(tempMatch.group(1))
        
            tempMatch = cigarSRe_left.search(F[5])
            if tempMatch is not None: left_clipping_primary = int(tempMatch.group(1))

        elif flags[6] == flags_SA[7] and flags[7] == flags_SA[6]:
            chr_pair = F[2]
            juncChr_pair = F[2]
            pos_pair = F[3]
            dir_pair = ("-" if flags[4] == "1" else "+")
            mq_pair = F[4]
            coverRegion_pair = myCigar.getCoverRegion(F[2], F[3], F[5])
        else:
            print >> sys.stderr, "The following read is both first and second reads at:" + '\n' + line


    if right_clipping_SA > min_clipping_size:

        juncDir_SA = "+"
        juncPos_SA = myCigar.getEndPos(F[3], F[5])


hIN = open(inputFile, 'r')

tempID = ""
tempLine = []
for line in hIN:
    line = line.rstrip('\n')
    F = line.split('\n')
    if tempID != F[0]:
        if tempID != "" and len(tempLine) == 3:
            printJuncInfo(tempLine)

        tempID = F[0]
        tempLine = []

    tempLine.append(line)

hIN.close()


if tempID != "" and len(tempLine) == 3:
    printJuncInfo(tempLine)

