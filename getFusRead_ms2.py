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
tempLine = []
for line in hIN:
    if line[0] == "@": continue
    line = line.rstrip('\n')
    F = line.split('\t')

    if tempID != F[0]:

        if fusFlag.count(0) != len(tempLine):

            print fusInfo
            # if tempID == "HWI-ST1165:87:AC0RYJACXX:2:1101:10138:158119": 
            #     print tempID
            # check the fusion validity
            chr_primary, pos_primary, dir_primary, chr_pair, pos_pair, dir_pair, chr_chimera, pos_chimera, dir_chimera = "*", "*", "*", "*", "*", "*", "*", "*", "*"
            for i in range(0, len(tempLine)):

                FF = tempLine[i].split('\t')
                flags = format(int(FF[1]), "#014b")[:1:-1]
                if flags[8] == "1" and fusFlag[i] != 0:
                    chr_chimera = FF[2]
                    pos_chimera = FF[3]
                    dir_chimera = ("+" if FF[4] != "1" else "-")
                elif flags[8] != "1" and fusFlag[i] != 0:
                    chr_primary = FF[2]
                    pos_primary = FF[3]
                    dir_primary = ("+" if FF[4] != "1" else "-")
                elif flags[8] != "1" and (fusFlag[i] == 0 or chr_primary != ""):
                    chr_pair = FF[2]
                    pos_pair = FF[3]
                    dir_pair = ("+" if FF[4] != "1" else "-")

            print '\t'.join([tempID, chr_primary, pos_primary, dir_primary, chr_pair, pos_pair, dir_pair, chr_chimera, pos_chimera, dir_chimera]) 


        tempID = F[0]
        fusFlag = []
        fusInfo = [] 
        tempLine = [] 
    
    tempLine.append(line)
    mFus = ReFus.search('\t'.join(F[11:]))
    if mFus is not None:
        fusFlag.append(1)
        fusInfo.append(mFus.group(0))
    else:
        fusFlag.append(0)

hIN.close()

