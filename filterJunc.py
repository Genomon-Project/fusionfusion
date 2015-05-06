#! /usr/local/bin/python

import sys, regions

inputFile = sys.argv[1]
minReadPairNum = int(sys.argv[2])
minValidThres = float(sys.argv[3])
minCoveredBase = int(sys.argv[4])
minSize = 1000

hIN = open(inputFile, 'r')

for line in hIN:
    F = line.rstrip('\n').split('\t')

    if F[0] == F[3] and abs(int(F[1]) - int(F[4])) < minSize: continue

    if F[1] == "5313724":
        pass

    # filter by the number of supporting read pairs
    # IDs = F[7].split(';')
    # uIDs = list(set(IDs))
    # if len(uIDs) < minReadPairNum: continue
    coveredRegion_primary = F[9].split(';')
    coveredRegion_pair = F[12].split(';')
    coveredRegion_SA = F[15].split(';')
    coveredRegion_meta = [coveredRegion_primary[i] + ';' + coveredRegion_pair[i] + ';' + coveredRegion_SA[i] for i in range(0, len(coveredRegion_primary))]
    uniqueCoverdRegion_meta = list(set(coveredRegion_meta))
    if len(uniqueCoverdRegion_meta) < minReadPairNum: continue

    # filter by the ratio of valid read pairs
    pairPos = F[17].split(';')
    if float(pairPos.count("0")) / float(len(pairPos)) > 1 - minValidThres: continue


    # check for the covered region
    coveredRegion_primary = F[9].split(';')
    coveredRegion_pair = F[12].split(';')
    coveredRegion_SA = F[15].split(';')
    pairPos = F[17].split(';')
    primaryPos = F[18].split(';')

    region1 =  regions.Regions()
    region2 =  regions.Regions()

    for i in range(0, len(coveredRegion_primary)):
        if primaryPos[i] == "1" and pairPos != "0":
            for reg in coveredRegion_primary[i].split(','):
                region1.addMerge(reg)
        elif primaryPos[i] == "2" and pairPos != "0":
            for reg in coveredRegion_primary[i].split(','):
                region2.addMerge(reg)

    for i in range(0, len(coveredRegion_SA)):
        if primaryPos[i] == "1" and pairPos != "0":
            for reg in coveredRegion_SA[i].split(','):
                region2.addMerge(reg)
        elif primaryPos[i] == "2" and pairPos != "0":
            for reg in coveredRegion_SA[i].split(','):
                region1.addMerge(reg)

    for i in range(0, len(coveredRegion_pair)):
        if (primaryPos[i] == "1" and pairPos[i] == "1") or (primaryPos[i] == "2" and pairPos[i] == "2"):
            for reg in coveredRegion_pair[i].split(','):
                region1.addMerge(reg)
        elif (primaryPos[i] == "1" and pairPos[i] == "2") or (primaryPos[i] == "2" and pairPos[i] == "1"):
            for reg in coveredRegion_pair[i].split(','):
                region2.addMerge(reg)

    region1.reduceMerge()
    region2.reduceMerge()

    if region1.regionSize() < minCoveredBase or region2.regionSize() < minCoveredBase: continue

    print '\t'.join(F)

hIN.close()
        
