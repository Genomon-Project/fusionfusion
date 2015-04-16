#! /usr/local/bin/python

import sys, re

inputFile = sys.argv[1]
checkMatchFile = sys.argv[2]
readNumThres = sys.argv[3]

ReContig = re.compile(r'(\w+):([\+\-])(\d+)\-(\w+):([\+\-])(\d+)_contig([12])')

hIN = open(checkMatchFile, 'r')
key2match1 = {}
key2match2 = {}
for line in hIN:
    F = line.rstrip('\n').split('\t')
    keyMatch = ReContig.match(F[0])

    chr1, strand1, pos1, chr2, strand2, pos2, contigNum = keyMatch.group(1), keyMatch.group(2), keyMatch.group(3), keyMatch.group(4), keyMatch.group(5), keyMatch.group(6), keyMatch.group(7)

    if contigNum == "1":
        key2match1['\t'.join([chr1, pos1, strand1, chr2, pos2, strand2])] = '\t'.join(F[1:4])

    if contigNum == "2":
        key2match2['\t'.join([chr1, pos1, strand1, chr2, pos2, strand2])] = '\t'.join(F[1:4])

hIN.close()


hIN = open(inputFile, 'r')
for line in hIN:
    F = line.rstrip('\n').split('\t')

    # temporary treatment
    if F[0] == "chrM" or F[3] == "chrM": continue

    key = '\t'.join(F[0:6])

    match1 = key2match1[key]
    match2 = key2match2[key]

    matches1 = match1.split('\t')
    matches2 = match2.split('\t')

    if matches1[0] == "---" or matches2[0] == "---": continue
    if matches1[1] != "---" or matches2[1] != "---": continue
    if float(matches1[2]) > 0 or float(matches2[2]) > 0: continue
    if int(F[7]) < int(readNumThres): continue

    print key + '\t' + F[6] + '\t' + F[7] + '\t' +  match1 + '\t' + match2

hIN.close()

