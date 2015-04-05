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
            print '\t'.join([tempKey, ';'.join(tempIDs),  ';'.join(tempMQ_primary), ';'.join(tempCover_primary), ';'.join(tempDir_primary), \
                             ';'.join(tempMQ_pair), ';'.join(tempCover_pair), ';'.join(tempDir_pair), \
                             ';'.join(tempMQ_SA), ';'.join(tempCover_SA), ';'.join(tempDir_SA), ';'.join(tempPairPos), ';'.join(tempPrimaryPos)])
            

        tempKey = key
        tempIDs = []
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

    tempIDs.append(F[6])
    tempMQ_primary.append(F[7])
    tempCover_primary.append(F[8])
    tempDir_primary.append(F[9])
    tempMQ_pair.append(F[10])
    tempCover_pair.append(F[11])
    tempDir_pair.append(F[12])
    tempMQ_SA.append(F[13])
    tempCover_SA.append(F[14])
    tempDir_SA.append(F[15])
    tempPairPos.append(F[16])
    tempPrimaryPos.append(F[17])

hIN.close()

print '\t'.join([tempKey, ';'.join(tempIDs),  ';'.join(tempMQ_primary), ';'.join(tempCover_primary), ';'.join(tempDir_primary), \
                 ';'.join(tempMQ_pair), ';'.join(tempCover_pair), ';'.join(tempDir_pair), \
                 ';'.join(tempMQ_SA), ';'.join(tempCover_SA), ';'.join(tempDir_SA), ';'.join(tempPairPos), ';'.join(tempPrimaryPos)])



