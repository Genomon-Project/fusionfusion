#! /usr/bin/env python

import sys, pysam

def filterAndAnnotation(inputFilePath, outputFilePath, Params):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    gene_bed = Params["gene_bed"] 
    filter_same_gene = Params["filter_same_gene"]

    gene_tb = pysam.TabixFile(gene_bed)

    for line in hIN:

        F = line.rstrip('\n').split('\t')
     
        ##########
        # check gene annotation for the side 1  
        tabixErrorFlag = 0
        try:
            records = gene_tb.fetch(F[0], int(F[1]) - 1, int(F[1]) + 1)
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorFlag = 1

        gene1 = [];
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                gene1.append(record[3])

        if len(gene1) == 0: gene1.append("---")
        gene1 = list(set(gene1))
        ##########

        ##########
        # check gene annotation for the side 2
        tabixErrorFlag = 0
        try:
            records = gene_tb.fetch(F[3], int(F[4]) - 1, int(F[4]) + 1) 
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorFlag = 1

        gene2 = [];
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                gene2.append(record[3])

        if len(gene2) == 0: gene2.append("---")
        gene2 = list(set(gene2))
        ##########

        sameGeneFlag = 0
        for g1 in gene1:
            for g2 in gene2:
                if g1 == g2 and g1 != "---": sameGeneFlag = 1

        if filter_same_gene == True and sameGeneFlag == 1: continue

        print >> hOUT, '\t'.join(F[0:8]) + '\t' + ';'.join(gene1) + '\t' + ';'.join(gene2) 
     


    hIN.close()
    hOUT.close()


