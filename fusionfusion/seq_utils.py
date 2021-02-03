#! /usr/bin/env python

import pysam

def getSeq(reference, region_tuple):
    region_tuple_sorted = sorted(region_tuple, key=lambda x: x[1])

    seq = ""    
    for reg in region_tuple_sorted:
        for line in pysam.faidx(reference, reg[0] + ":" + str(reg[1]) + "-" + str(reg[2]), split_lines=True):
            if line[0] == ">": continue
            seq += line

    return(seq)

def reverseComplement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return("".join(complement.get(base, base) for base in reversed(seq)))


