#! /usr/local/bin/python

import pysam

def getSeq(reference, region_tuple):
    region_tuple_sorted = sorted(region_tuple, key=lambda x: x[1])

    seq = ""    
    for reg in region_tuple_sorted:
        for item in pysam.faidx(reference, reg[0] + ":" + str(reg[1]) + "-" + str(reg[2])):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n')

    return(seq)


if __name__ == "__main__":
    import sys, pysam
    print getSeq(sys.argv[1], [(sys.argv[2], sys.argv[3], sys.argv[4])])

