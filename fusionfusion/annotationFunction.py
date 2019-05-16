#! /usr/bin/env python

from __future__ import print_function
import sys, pysam, os, subprocess
import annot_utils.gene
import annot_utils.exon

# import config
from .config import *

junction_margin = 5


def get_gene_info(chr, pos, ref_gene_tb, ens_gene_tb):

    # check gene annotation for refGene 
    tabixErrorFlag = 0
    try:
        records = ref_gene_tb.fetch(chr, int(pos) - 1, int(pos) + 1)
    except Exception as inst:
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1
        
    gene = [];
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            gene.append(record[3])

    # if any gene cannot be found in refGene, then search ensGene
    if len(gene) == 0:
        tabixErrorFlag = 0
        try:
            records = ens_gene_tb.fetch(chr, int(pos) - 1, int(pos) + 1)
        except Exception as inst:
            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorFlag = 1
           
        """ 
        # for ensGene, just the longest gene is shown
        temp_length = 0
        temp_gene = ""
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                if int(record[4]) > temp_length:
                     temp_gene = record[3]
    
            if temp_gene != "": gene.append(temp_gene)
        """

        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                gene.append(record[3])

 
    if len(gene) == 0: gene.append("---")

    return list(set(gene))


def get_junc_info(chr, pos, ref_exon_tb, ens_exon_tb, junction_margin):

    # check exon annotation for refGene 
    tabixErrorFlag = 0
    try:
        records = ref_exon_tb.fetch(chr, int(pos) - junction_margin, int(pos) + junction_margin)
    except Exception as inst:
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1
        
    junction = []
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if abs(int(pos) - int(record[1])) < junction_margin:
                if record[5] == "+": junction.append(record[3] + ".start")
                if record[5] == "-": junction.append(record[3] + ".end")
            if abs(int(pos) - int(record[2])) < junction_margin:
                if record[5] == "+": junction.append(record[3] + ".end")
                if record[5] == "-": junction.append(record[3] + ".start")

    # if any exon-intron junction cannot be found in refGene, then search ensGene
    if len(junction) == 0: 
        tabixErrorFlag = 0
        try:
            records = ens_exon_tb.fetch(chr, int(pos) - junction_margin, int(pos) + junction_margin)
        except Exception as inst:
            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorFlag = 1

        """             
        # for ensGene, just the longest gene is shown
        temp_length = 0
        temp_junc = ""
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                # if int(record[4]) > temp_length:
                if True:
                    if abs(int(pos) - int(record[1])) < junction_margin:
                        if record[5] == "+": temp_junc = record[3] + ".start"
                        if record[5] == "-": temp_junc = record[3] + ".end"
                    if abs(int(pos) - int(record[2])) < junction_margin: 
                        if record[5] == "+": temp_junc = record[3] + ".end"
                        if record[5] == "-": temp_junc = record[3] + ".start"
                
            if temp_junc != "": junction.append(temp_junc)
        """

        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                if abs(int(pos) - int(record[1])) < junction_margin:
                    if record[5] == "+": junction.append(record[3] + ".start")
                    if record[5] == "-": junction.append(record[3] + ".end")
                if abs(int(pos) - int(record[2])) < junction_margin: 
                    if record[5] == "+": junction.append(record[3] + ".end")
                    if record[5] == "-": junction.append(record[3] + ".start")
    
                
    if len(junction) == 0: junction.append("---")
    
    return list(set(junction))



def filterAndAnnotation(inputFilePath, outputFilePath, genome_id, is_grc):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    # annotation_dir = config.param_conf.get("annotation", "annotation_dir")
    # filter_same_gene = config.param_conf.getboolean("filter_condition", "filter_same_gene")
    # annotation_dir = param_conf.resource_dir
    filter_same_gene = param_conf.filter_same_gene

    """
    # old procedure
    # ref_gene_bed = annotation_dir + "/refGene.bed.gz"
    ref_exon_bed = annotation_dir + "/refExon.bed.gz"
    ens_gene_bed = annotation_dir + "/ensGene.bed.gz"
    ens_exon_bed = annotation_dir + "/ensExon.bed.gz"
    grch2ucsc_file = annotation_dir + "/grch2ucsc.txt"

    # relationship between CRCh and UCSC chromosome names
    grch2ucsc = {}
    with open(grch2ucsc_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            grch2ucsc[F[0]] = F[1]

    ref_gene_tb = pysam.TabixFile(ref_gene_bed)
    ref_exon_tb = pysam.TabixFile(ref_exon_bed)
    ens_gene_tb = pysam.TabixFile(ens_gene_bed)
    ens_exon_tb = pysam.TabixFile(ens_exon_bed)
    """

    annot_utils.gene.make_gene_info(outputFilePath + ".tmp.refGene.bed.gz", "refseq", genome_id, is_grc, False) 
    annot_utils.gene.make_gene_info(outputFilePath + ".tmp.ensGene.bed.gz", "gencode", genome_id, is_grc, False)
    annot_utils.exon.make_exon_info(outputFilePath + ".tmp.refExon.bed.gz", "refseq", genome_id, is_grc, False)
    annot_utils.exon.make_exon_info(outputFilePath + ".tmp.ensExon.bed.gz", "gencode", genome_id, is_grc, False)
 
    ref_gene_tb = pysam.TabixFile(outputFilePath + ".tmp.refGene.bed.gz")
    ens_gene_tb = pysam.TabixFile(outputFilePath + ".tmp.ensGene.bed.gz")
    ref_exon_tb = pysam.TabixFile(outputFilePath + ".tmp.refExon.bed.gz")
    ens_exon_tb = pysam.TabixFile(outputFilePath + ".tmp.ensExon.bed.gz")

    for line in hIN:

        F = line.rstrip('\n').split('\t')
    
        # check gene annotation for the side 1
        gene1 = get_gene_info(F[0], F[1], ref_gene_tb, ens_gene_tb)
        
        # check gene annotation for the side 2
        gene2 = get_gene_info(F[3], F[4], ref_gene_tb, ens_gene_tb)

        # check exon-intron junction annotation for the side 1 
        junction1 = get_junc_info(F[0], F[1], ref_exon_tb, ens_exon_tb, junction_margin) 

        # check exon-intron junction annotation for the side 2 
        junction2 = get_junc_info(F[3], F[4], ref_exon_tb, ens_exon_tb, junction_margin)

        sameGeneFlag = 0
        for g1 in gene1:
            for g2 in gene2:
                if g1 == g2 and g1 != "---": sameGeneFlag = 1

        if filter_same_gene == True and sameGeneFlag == 1: continue

        print('\t'.join(F[0:8]) + '\t' + ';'.join(gene1) + '\t' + ';'.join(junction1) + '\t' + ';'.join(gene2) + '\t' + ';'.join(junction2) + '\t' + \
              F[11] + '\t' + F[12] + '\t' + F[16] + '\t' + F[17], file = hOUT)

    hIN.close()
    hOUT.close()


    subprocess.check_call(["rm", "-rf", outputFilePath + ".tmp.refGene.bed.gz"])
    subprocess.check_call(["rm", "-rf", outputFilePath + ".tmp.ensGene.bed.gz"])
    subprocess.check_call(["rm", "-rf", outputFilePath + ".tmp.refExon.bed.gz"])
    subprocess.check_call(["rm", "-rf", outputFilePath + ".tmp.ensExon.bed.gz"])

    subprocess.check_call(["rm", "-rf", outputFilePath + ".tmp.refGene.bed.gz.tbi"])
    subprocess.check_call(["rm", "-rf", outputFilePath + ".tmp.ensGene.bed.gz.tbi"])
    subprocess.check_call(["rm", "-rf", outputFilePath + ".tmp.refExon.bed.gz.tbi"])
    subprocess.check_call(["rm", "-rf", outputFilePath + ".tmp.ensExon.bed.gz.tbi"])



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

    fus_keys = list(set(list(fus2count_star) + list(fus2count_ms2) + list(fus2count_th2)))
    hOUT = open(output_file_path + ".unsorted.tmp", 'w')

    for fus in fus_keys:
        count_star = fus2count_star[fus] if fus in fus2count_star else "---"
        count_ms2 = fus2count_ms2[fus] if fus in fus2count_ms2 else "---"
        count_th2 = fus2count_th2[fus] if fus in fus2count_th2 else "---"

        print_line = fus
        if label_ms2 == 1: print_line = print_line + '\t' + count_ms2
        if label_star == 1: print_line = print_line + '\t' + count_star
        if label_th2 == 1: print_line = print_line + '\t' + count_th2

        print(print_line, file = hOUT)

    hOUT.close()

    hOUT = open(output_file_path, 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k4,4", "-k5,5n", output_file_path + '.unsorted.tmp'], stdout = hOUT)
    hOUT.close()

    subprocess.check_call(["rm", "-rf", output_file_path + ".unsorted.tmp"])
 
