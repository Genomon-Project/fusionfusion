#! /usr/local/bin/python

import sys, tabix

inputFile = sys.argv[1]
geneFile = sys.argv[2]

hIN = open(inputFile, 'r')
gene_tb = tabix.open(geneFile)

for line in hIN:
    F = line.rstrip('\n').split('\t')
 
    ##########
    # check gene annotation for the side 1  
    tabixErrorFlag = 0
    try:
        records = gene_tb.query(F[0], int(F[1]) - 1, int(F[1]))
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    gene1 = [];
    if tabixErrorFlag == 0:
        for record in records:
            gene1.append(record[3])

    if len(gene1) == 0: gene1.append("---")
    gene1 = list(set(gene1))
    ##########

    ##########
    # check gene annotation for the side 2
    tabixErrorFlag = 0
    try:
        records = gene_tb.query(F[3], int(F[4]) - 1, int(F[4]))
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    gene2 = [];
    if tabixErrorFlag == 0:
        for record in records:
            gene2.append(record[3])

    if len(gene2) == 0: gene2.append("---")
    gene2 = list(set(gene2))
    ##########

    sameGeneFlag = 0
    for g1 in gene1:
        for g2 in gene2:
            if g1 == g2 and g1 != "---": sameGeneFlag = 1

    # if sameGeneFlag == 1: continue

    print '\t'.join(F[0:8]) + '\t' + ';'.join(gene1) + '\t' + ';'.join(gene2) 
 


hIN.close()



