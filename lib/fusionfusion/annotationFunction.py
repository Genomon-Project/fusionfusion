#! /usr/bin/env python

import sys, pysam, os

import config
junction_margin = 5

def filterAndAnnotation(inputFilePath, outputFilePath):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    gene_bed = config.param_conf.get("annotation", "gene_bed")
    exon_bed = config.param_conf.get("annotation", "exon_bed")
    filter_same_gene = config.param_conf.getboolean("filter_condition", "filter_same_gene")

    gene_tb = pysam.TabixFile(gene_bed)
    exon_tb = pysam.TabixFile(exon_bed)

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
        # check exon-intron junction annotation for the side 1  
        tabixErrorFlag = 0
        try:
            records = exon_tb.fetch(F[0], int(F[1]) - junction_margin, int(F[1]) + junction_margin)
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorFlag = 1
        
        junction1 = [] 
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                if abs(int(F[1]) - int(record[1])) < junction_margin:
                    if record[5] == "+": junction1.append(record[3] + ".start")
                    if record[5] == "-": junction1.append(record[3] + ".end")
                if abs(int(F[1]) - int(record[2])) < junction_margin:
                    if record[5] == "+": junction1.append(record[3] + ".end")
                    if record[5] == "-": junction1.append(record[3] + ".start") 

        if len(junction1) == 0: junction1.append("---")
        junction1 = list(set(junction1))
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

        ##########
        # check exon-intron junction annotation for the side 2  
        tabixErrorFlag = 0
        try:
            records = exon_tb.fetch(F[3], int(F[4]) - junction_margin, int(F[4]) + junction_margin)
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorFlag = 1

        junction2 = []
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                if abs(int(F[4]) - int(record[1])) < junction_margin:
                    if record[5] == "+": junction2.append(record[3] + ".start")
                    if record[5] == "-": junction2.append(record[3] + ".end")
                if abs(int(F[4]) - int(record[2])) < junction_margin:
                    if record[5] == "+": junction2.append(record[3] + ".end")
                    if record[5] == "-": junction2.append(record[3] + ".start")

        if len(junction2) == 0: junction2.append("---")
        junction2 = list(set(junction2))
        ##########


        sameGeneFlag = 0
        for g1 in gene1:
            for g2 in gene2:
                if g1 == g2 and g1 != "---": sameGeneFlag = 1

        if filter_same_gene == True and sameGeneFlag == 1: continue

        print >> hOUT, '\t'.join(F[0:8]) + '\t' + ';'.join(gene1) + '\t' + ';'.join(junction1) + '\t' + ';'.join(gene2) + '\t' + ';'.join(junction2) + '\t' + \
                         F[11] + '\t' + F[12] + '\t' + F[16] + '\t' + F[17]

    hIN.close()
    hOUT.close()


def merge_fusion_result(input_dir, output_file_path):

    fus2count_ms2 = {}
    fus2count_star = {}
    fus2count_th2 = {}

    label_ms2, label_star, label_th2 = 0, 0, 0

    if os.path.exists(input_dir + "/ms2.fusion.result.txt"):
        label_ms2 = 1
        with open(input_dir + "/ms2.fusion.result.txt") as hIN:
            for line in hIN:
                F = line.rstrip('\n').split('\t')
                key = '\t'.join(F[0:7] + F[8:12])
                fus2count_ms2[key] = F[7]

    if os.path.exists(input_dir + "/star.fusion.result.txt"): 
        label_star = 1
        with open(input_dir + "/star.fusion.result.txt") as hIN:
            for line in hIN:
                F = line.rstrip('\n').split('\t')
                key = '\t'.join(F[0:7] + F[8:12])
                fus2count_star[key] = F[7]

    if os.path.exists(input_dir + "/th2.fusion.result.txt"):
        label_th2 = 1 
        with open(input_dir + "/th2.fusion.result.txt") as hIN:
            for line in hIN:
                F = line.rstrip('\n').split('\t')
                key = '\t'.join(F[0:7] + F[8:12])
                fus2count_th2[key] = F[7]

    fus_keys = list(set(fus2count_star.keys() + fus2count_ms2.keys() + fus2count_th2.keys()))
    hOUT = open(output_file_path, 'w')

    for fus in fus_keys:
        count_star = fus2count_star[fus] if fus in fus2count_star else "---"
        count_ms2 = fus2count_ms2[fus] if fus in fus2count_ms2 else "---"
        count_th2 = fus2count_th2[fus] if fus in fus2count_th2 else "---"

        print_line = fus
        if label_ms2 == 1: print_line = print_line + '\t' + count_ms2
        if label_star == 1: print_line = print_line + '\t' + count_star
        if label_th2 == 1: print_line = print_line + '\t' + count_th2

        print >> hOUT, print_line

    hOUT.close()
