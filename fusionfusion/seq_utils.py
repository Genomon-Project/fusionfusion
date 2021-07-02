#! /usr/bin/env python

import pysam

def getSeq(reference, region_tuple):
    region_tuple_sorted = sorted(region_tuple, key=lambda x: x[1])

    seq = ""    
    for reg in region_tuple_sorted:
        lines = []
        try:
            lines = pysam.faidx(reference, reg[0] + ":" + str(reg[1]) + "-" + str(reg[2]), split_lines=True)
        except pysam.utils.SamtoolsError:
            # With latest versions of pysam, faidx raises an error if it gets fed with
            # a contradictory range (such as `chr:start-end` where start > end).
            # This try-except statement is for handling such cases.
            pass
        for line in lines:
            if line[0] == ">": continue
            seq += line

    return(seq)

def reverseComplement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return("".join(complement.get(base, base) for base in reversed(seq)))


