#! /usr/local/bin/python

"""
    script for collecting short reads supporting fusion candidates in MapSplice2 sam file
"""

import sys, re

inputFile = sys.argv[1]

hIN = open(inputFile, 'r')

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

            # check the fusion validity
            ufusInfo = list(set(fusInfo))
            for fus in ufusInfo:
                if fus == 0: continue

                chr_primary, pos_primary, dir_primary, chr_pair, pos_pair, dir_pair, chr_chimera, pos_chimera, dir_chimera = "*", "*", "*", "*", "*", "*", "*", "*", "*"
                fusOrder = 0
                for i in range(0, len(tempLine)):

                    FF = tempLine[i].split('\t')
                    flags = format(int(FF[1]), "#014b")[:1:-1]
                    if flags[8] == "1" and fusInfo[i] == fus:
                        chr_chimera = FF[2]
                        pos_chimera = FF[3]
                        dir_chimera = ("+" if flags[4] != "1" else "-")
                        if fusOrder == 0: fusOrder = -1
                    elif flags[8] != "1" and fusInfo[i] == fus:
                        chr_primary = FF[2]
                        pos_primary = FF[3]
                        dir_primary = ("+" if flags[4] != "1" else "-")
                        if fusOrder == 0: fusOrder = 1
                    elif flags[8] != "1" and (fusInfo[i] != fus or chr_primary != "*"):
                        chr_pair = FF[2]
                        pos_pair = FF[3]
                        dir_pair = ("+" if flags[4] != "1" else "-")

                    fusSplit = fus.split(',')
                    if fusOrder == 1:
                        breakPos_primary, breakPos_chimera, breakDir_primary, breakDir_chimera = fusSplit[0], fusSplit[1], fusSplit[2], fusSplit[3]
                    else:
                        breakPos_primary, breakPos_chimera, breakDir_primary, breakDir_chimera = fusSplit[1], fusSplit[0], fusSplit[3], fusSplit[2]
                        breakDir_primary = ("+" if breakDir_primary == "-" else "-")
                        breakDir_chimera = ("+" if breakDir_chimera == "-" else "-")

                    validFlag = 0
                    if pos_pair != "*":
                        if breakDir_primary == "+" and dir_pair == "+" and int(breakPos_primary) - 500000 <= int(pos_pair) <= int(breakPos_primary): validFlag = 1
                        if breakDir_primary == "-" and dir_pair == "-" and int(breakPos_primary) <= int(pos_pair) <= int(breakPos_primary) + 500000: validFlag = 1
                        if breakDir_chimera == "+" and dir_pair == "-" and int(breakPos_chimera) <= int(pos_pair) <= int(breakPos_chimera) + 500000: validFlag = 1
                        if breakDir_chimera == "-" and dir_pair == "+" and int(breakPos_chimera) - 500000 <= int(pos_pair) <= int(breakPos_chimera): validFlag = 1

                print '\t'.join([breakPos_primary, breakPos_chimera, breakDir_primary, breakDir_chimera, tempID, chr_primary, pos_primary, chr_chimera, pos_chimera, chr_pair, pos_pair, dir_pair, str(validFlag)] )

                # 133729451,23632600,-,-  1       HWI-ST1165:87:AC0RYJACXX:2:1101:3462:76598      chr9    133729451       -       chr22   23631767        +       chr22   23632550        -


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

