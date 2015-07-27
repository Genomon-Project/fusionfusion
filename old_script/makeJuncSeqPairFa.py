#! /usr/local/bin/python

import sys

inputFile = sys.argv[1]

hIN = open(inputFile, 'r')

for line in hIN:
    F = line.rstrip('\n').split('\t')

    print '>' + F[0] + ":" + F[2] + F[1] + "-" + F[3] + ":" + F[5] + F[4] + "_contig" + '1'
    print F[10]

    print '>' + F[0] + ":" + F[2] + F[1] + "-" + F[3] + ":" + F[5] + F[4] + "_contig" + '2'
    print F[11]

hIN.close()

