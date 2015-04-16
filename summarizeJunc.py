#! /usr/local/bin/python

import sys

inputFile = sys.argv[1]

hIN = open(inputFile, 'r')

tempKey = ""
for line in hIN:
    F = line.rstrip('\n').split('\t')
    key = '\t'.join(F[0:6])

    if tempKey != key:
        if tempKey != "":
            print '\t'.join([tempKey, ';'.join(tempInseq), ';'.join(tempIDs),  ';'.join(tempMQ_primary), ';'.join(tempCover_primary), ';'.join(tempDir_primary), \
                             ';'.join(tempMQ_pair), ';'.join(tempCover_pair), ';'.join(tempDir_pair), \
                             ';'.join(tempMQ_SA), ';'.join(tempCover_SA), ';'.join(tempDir_SA), ';'.join(tempPairPos), ';'.join(tempPrimaryPos)])
            

        tempKey = key
        tempIDs = []
        tempInseq = []
        tempMQ_primary = []
        tempCover_primary = []
        tempDir_primary = []
        tempMQ_pair = []
        tempCover_pair = []
        tempDir_pair = []
        tempMQ_SA = []
        tempCover_SA = []
        tempDir_SA = []
        tempPairPos = []
        tempPrimaryPos = []

    tempInseq.append(F[6])
    tempIDs.append(F[7])
    tempMQ_primary.append(F[8])
    tempCover_primary.append(F[9])
    tempDir_primary.append(F[10])
    tempMQ_pair.append(F[11])
    tempCover_pair.append(F[12])
    tempDir_pair.append(F[13])
    tempMQ_SA.append(F[14])
    tempCover_SA.append(F[15])
    tempDir_SA.append(F[16])
    tempPairPos.append(F[17])
    tempPrimaryPos.append(F[18])

hIN.close()

print '\t'.join([tempKey, ';'.join(tempInseq), ';'.join(tempIDs),  ';'.join(tempMQ_primary), ';'.join(tempCover_primary), ';'.join(tempDir_primary), \
                 ';'.join(tempMQ_pair), ';'.join(tempCover_pair), ';'.join(tempDir_pair), \
                 ';'.join(tempMQ_SA), ';'.join(tempCover_SA), ';'.join(tempDir_SA), ';'.join(tempPairPos), ';'.join(tempPrimaryPos)])



