#! /usr/bin/env python

import re

reCigar = re.compile(r'(\d+)([DIMNS])')

# get the region covered by the alignment
def getCoverRegion(chr, pos, cigar):

    tempStart = int(pos)
    tempSize = 0 
    coverRegion = []
    for m in reCigar.finditer(cigar):
        if m.group(2) == "N":
            coverRegion.append(chr + ":" + str(tempStart) + "-" + str(tempStart + tempSize - 1))
            tempStart = tempStart + tempSize + int(m.group(1))
            tempSize = 0 
        elif m.group(2) in ["M", "D"]:
            tempSize = int(tempSize) + int(m.group(1))

    coverRegion.append(chr + ":" + str(tempStart) + "-" + str(tempStart + tempSize - 1))

    return ','.join(coverRegion)


# get the end positon
def getEndPos(pos, cigar):

    tempPos = int(pos) - 1
    for m in reCigar.finditer(cigar):
        if m.group(2) in ["M", "N", "D"]:
            tempPos = tempPos + int(m.group(1))

    return tempPos 


