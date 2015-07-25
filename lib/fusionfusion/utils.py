#!/usr/bin/env python

import sys, os, subprocess, pysam


def sortBedpe(inputFile, outputFile):

    hOUT = open(outputFile, "w")
    subprocess.call(["sort", "-k1,1", "-k2,2n", "-k4,4", "-k5,5n", inputFile], stdout = hOUT)
    hOUT.close()

def make_directory(inputDir):
    """
    make input directory if it does not exist.
    """
    if not os.path.exists(inputDir):
        os.makedirs(inputDir)

 
def getSeq(reference, region_tuple):
    region_tuple_sorted = sorted(region_tuple, key=lambda x: x[1])

    seq = ""
    for reg in region_tuple_sorted:
        for item in pysam.faidx(reference, reg[0] + ":" + str(reg[1]) + "-" + str(reg[2])):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n')

    return(seq)

def reverseComplement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return("".join(complement.get(base, base) for base in reversed(seq)))

