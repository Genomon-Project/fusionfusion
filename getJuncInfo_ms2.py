#! /usr/local/bin/python

"""
    script for collecting short reads supporting fusion candidates in MapSplice2 sam file
"""

import sys, re, getCoverRegion

inputFile = sys.argv[1]

hIN = open(inputFile, 'r')


def printFusInfo(tempLine, fusInfo):

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
                cover_chimera = getCoverRegion.getCoverRegion(FF[2], FF[3], FF[5])
                if fusOrder == 0: fusOrder = -1
            elif flags[8] != "1" and fusInfo[i] == fus:
                chr_primary = FF[2]
                pos_primary = FF[3]
                dir_primary = ("+" if flags[4] != "1" else "-")
                mq_primary = FF[4]
                cover_primary = getCoverRegion.getCoverRegion(FF[2], FF[3], FF[5])
                if fusOrder == 0: fusOrder = 1
            elif flags[8] != "1" and flags[2] != "1" and (fusInfo[i] != fus or chr_primary != "*"):
                chr_pair = FF[2]
                pos_pair = FF[3]
                dir_pair = ("+" if flags[4] != "1" else "-")
                mq_pair = FF[4]
                cover_pair = getCoverRegion.getCoverRegion(FF[2], FF[3], FF[5])

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
            if breakDir_primary == "+" and dir_pair == "+" and int(breakPos_primary) - 500000 <= int(pos_pair) <= int(breakPos_primary): pairPos = 1
            if breakDir_primary == "-" and dir_pair == "-" and int(breakPos_primary) <= int(pos_pair) <= int(breakPos_primary) + 500000: pairPos = 1
            if breakDir_chimera == "+" and dir_pair == "+" and int(breakPos_chimera) - 500000 <= int(pos_pair) <= int(breakPos_chimera): pairPos = 2
            if breakDir_chimera == "-" and dir_pair == "-" and int(breakPos_chimera) <= int(pos_pair) <= int(breakPos_chimera) + 500000: pairPos = 2

        if chr_primary < chr_chimera or chr_primary == chr_chimera and breakPos_primary <= breakPos_chimera:
            print '\t'.join([chr_primary, breakPos_primary, breakDir_primary, chr_chimera, breakPos_chimera, breakDir_chimera, tempID, \
                             mq_primary, cover_primary, dir_primary, mq_pair, cover_pair, dir_pair, \
                             mq_chimera, cover_chimera, dir_chimera, str(pairPos), "1"] )
        else:
            print '\t'.join([chr_chimera, breakPos_chimera, breakDir_chimera, chr_primary, breakPos_primary, breakDir_primary, tempID, \
                             mq_primary, cover_primary, dir_primary, mq_pair, cover_pair, dir_pair, \
                             mq_chimera, cover_chimera, dir_chimera, str(pairPos), "2"] )



ReFus = re.compile('FUS_(\d+)_(\d+)\(([\-\+])([\-\+])\)')

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
            printFusInfo(tempLine, fusInfo)

        tempID = F[0]
        fusFlag = []
        fusInfo = [] 
        tempLine = [] 
    
    tempLine.append(line)
    mFus = ReFus.search('\t'.join(F[11:]))
    if mFus is not None:
        fusFlag.append(1)
        fusInfo.append(','.join([mFus.group(1), mFus.group(2), mFus.group(3), mFus.group(4)]))
    else:
        fusFlag.append(0)
        fusInfo.append(0)

hIN.close()


if fusInfo.count(0) != len(tempLine):
    printFusInfo(tempLine, fusInfo)


