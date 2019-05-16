#! /usr/bin/env python

from __future__ import print_function

import re
import pysam

# import config
from .config import *
from . import regions, seq_utils, region_utils

regRe = re.compile(r'([^ \t\n\r\f\v,]+):(\d+)\-(\d+)')
ReContig = re.compile(r'([^ \t\n\r\f\v,]+):([\+\-])(\d+)\-([^ \t\n\r\f\v,]+):([\+\-])(\d+)_contig([12])')

def filterCoverRegion(inputFilePath, outputFilePath):

    # min_read_pair_num = config.param_conf.getint("filter_condition", "min_read_pair_num")
    # min_valid_read_pair_ratio = config.param_conf.getfloat("filter_condition", "min_valid_read_pair_ratio")
    # min_cover_size = config.param_conf.getint("filter_condition", "min_cover_size")
    # min_chimeric_size = config.param_conf.getint("filter_condition", "min_chimeric_size")
    # anchor_size_thres = config.param_conf.getint("filter_condition", "anchor_size_thres")
    min_read_pair_num = param_conf.min_read_pair_num
    min_valid_read_pair_ratio = param_conf.min_valid_read_pair_ratio
    min_cover_size = param_conf.min_cover_size
    min_chimeric_size = param_conf.min_chimeric_size
    anchor_size_thres = param_conf.anchor_size_thres


    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    for line in hIN:
        F = line.rstrip('\n').split('\t')

        if F[0] == F[3] and abs(int(F[1]) - int(F[4])) < min_chimeric_size: continue

        coveredRegion_primary = F[9].split(';')
        coveredRegion_pair = F[12].split(';')
        coveredRegion_SA = F[15].split(';')

        # check the number of unique read pairs
        coveredRegion_meta = [coveredRegion_primary[i] + ';' + coveredRegion_pair[i] + ';' + coveredRegion_SA[i] for i in range(0, len(coveredRegion_primary))]
        uniqueCoverdRegion_meta = list(set(coveredRegion_meta))
        if len(uniqueCoverdRegion_meta) < min_read_pair_num: continue

        # check the maximum anchor size
        coverRegionSize_primary = list(map(region_utils.getCoverSize, coveredRegion_primary))
        coverRegionSize_SA = list(map(region_utils.getCoverSize, coveredRegion_SA))
        anchor_size = [min(coverRegionSize_primary[i], coverRegionSize_SA[i]) for i in range(len(coverRegionSize_primary))]
        if max(anchor_size) < anchor_size_thres: continue

            
        # filter by the ratio of valid read pairs
        pairPos = F[17].split(';')
        if float(pairPos.count("0")) / float(len(pairPos)) > 1 - min_valid_read_pair_ratio: continue

        # check for the covered region
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

        if region1.regionSize() < min_cover_size or region2.regionSize() < min_cover_size: continue

        print('\t'.join(F), file = hOUT)

    hIN.close()
    hOUT.close()



def filterPoolControl(input_file, output_file, control_file):

    control_db = pysam.TabixFile(control_file)

    hout = open(output_file, 'w') 

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
 
            # skip if the junction is included in the control file
            tabixErrorFlag = 0
            if control_file is not None:
                try:
                    records = control_db.fetch(F[0], int(F[1]) - 5, int(F[1]) + 5)
                except Exception as inst:
                    # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    # tabixErrorMsg = str(inst.args)
                    tabixErrorFlag = 1

            control_flag = 0;
            if tabixErrorFlag == 0:
                 for record_line in records:
                    record = record_line.split('\t')
                    if F[0] == record[0] and F[1] == record[1] and F[2] == record[2] and \
                       F[3] == record[3] and F[4] == record[4] and F[5] == record[5]:
                        control_flag = 1

            if control_flag == 0:
                print('\t'.join(F), file = hout)

    hout.close()


def extractSplicingPattern(inputFilePath, outputFilePath):

    # reference_genome = config.param_conf.get("alignment", "reference_genome")
    reference_genome = param_conf.reference_genome

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    for line in hIN:
        F = line.rstrip('\n').split('\t')

        # check for the covered region
        coveredRegion_primary = F[9].split(';')
        coveredRegion_pair = F[12].split(';')
        coveredRegion_SA = F[15].split(';')
        pairPos = F[17].split(';')
        primaryPos = F[18].split(';')


        ####################
        splice2count1 = {}
        splice2count2 = {}

        # get the splicing position for the primary reads 
        for i in range(0, len(coveredRegion_primary)):
            regs = coveredRegion_primary[i].split(',')
            if len(regs) == 1: continue

            chr_reg = ""
            start_reg = []
            end_reg = [] 
            for reg in regs:
                mReg = regRe.match(reg)
                chr_reg = mReg.group(1)
                start_reg.append(int(mReg.group(2)))
                end_reg.append(int(mReg.group(3)))
            
            start_reg.sort()
            end_reg.sort()
        
            for j in range(0, len(start_reg) - 1):
                regKey = (chr_reg, int(end_reg[j]), int(start_reg[j + 1]))
                if primaryPos[i] == "1" and pairPos != "0":
                    splice2count1[regKey] = (splice2count1[regKey] + 1 if regKey in splice2count1 else 1)
                elif primaryPos[i] == "2" and pairPos != "0":
                    splice2count2[regKey] = (splice2count2[regKey] + 1 if regKey in splice2count2 else 1)


        # get the splicing position for the supplementary reads
        for i in range(0, len(coveredRegion_SA)):
            regs = coveredRegion_SA[i].split(',')
            if len(regs) == 1: continue

            chr_reg = ""
            start_reg = []
            end_reg = []
            for reg in regs:
                mReg = regRe.match(reg)
                chr_reg = mReg.group(1)
                start_reg.append(int(mReg.group(2)))
                end_reg.append(int(mReg.group(3)))

            start_reg.sort()
            end_reg.sort()

            for j in range(0, len(start_reg) - 1):
                regKey = (chr_reg, int(end_reg[j]), int(start_reg[j + 1]))
                if primaryPos[i] == "1" and pairPos != "0":
                    splice2count2[regKey] = (splice2count2[regKey] + 1 if regKey in splice2count2 else 1)
                elif primaryPos[i] == "2" and pairPos != "0":
                    splice2count1[regKey] = (splice2count1[regKey] + 1 if regKey in splice2count1 else 1)


        # get the splicing position for the pair reads
        for i in range(0, len(coveredRegion_pair)):
            regs = coveredRegion_pair[i].split(',')
            if len(regs) == 1: continue

            chr_reg = ""
            start_reg = []
            end_reg = []
            for reg in regs:
                mReg = regRe.match(reg)
                chr_reg = mReg.group(1)
                start_reg.append(int(mReg.group(2)))
                end_reg.append(int(mReg.group(3)))

            start_reg.sort()
            end_reg.sort()

            for j in range(0, len(start_reg) - 1):
                regKey = (chr_reg, int(end_reg[j]), int(start_reg[j + 1]))
                if (primaryPos[i] == "1" and pairPos[i] == "1") or (primaryPos[i] == "2" and pairPos[i] == "2"):
                    splice2count1[regKey] = (splice2count1[regKey] + 1 if regKey in splice2count1 else 1)
                elif (primaryPos[i] == "1" and pairPos[i] == "2") or (primaryPos[i] == "2" and pairPos[i] == "1"):
                    splice2count2[regKey] = (splice2count2[regKey] + 1 if regKey in splice2count2 else 1)
        ####################

        ####################
        splice2count1_ref = {}
        splice2count2_ref = {}
        for key in splice2count1:
            splice2count1_ref[key] = 0

        for key in splice2count2:
            splice2count2_ref[key] = 0

        
        # check about each splicing position for the primary reads 
        for i in range(0, len(coveredRegion_primary)):
            for reg in coveredRegion_primary[i].split(','):
                mReg = regRe.match(reg)
                chr_reg, start_reg, end_reg  = mReg.group(1), int(mReg.group(2)), int(mReg.group(3))
     
                if primaryPos[i] == "1" and pairPos != "0":
                    for key in splice2count1:
                        chr_key, start_key, end_key  = key[0], key[1], key[2]

                        # if the given region cover arround the spliced sites
                        if F[2] == "+" and start_reg <= end_key - 5 and end_key + 5 <= end_reg: splice2count1_ref[key] = splice2count1_ref[key] + 1
                        if F[2] == "-" and start_reg <= start_key - 5 and start_key + 5 <= end_reg: splice2count1_ref[key] = splice2count1_ref[key] + 1

                if primaryPos[i] == "2" and pairPos != "0":
                    for key in splice2count2:
                        chr_key, start_key, end_key  = key[0], key[1], key[2]

                        # if the given region cover arround the spliced sites
                        if F[2] == "+" and start_reg <= end_key - 5 and end_key + 5 <= end_reg: splice2count2_ref[key] = splice2count2_ref[key] + 1
                        if F[2] == "-" and start_reg <= start_key - 5 and start_key + 5 <= end_reg: splice2count2_ref[key] = splice2count2_ref[key] + 1


        # check about each splicing position for the supplementary reads 
        for i in range(0, len(coveredRegion_SA)):
            for reg in coveredRegion_SA[i].split(','):
                mReg = regRe.match(reg)
                chr_reg, start_reg, end_reg  = mReg.group(1), int(mReg.group(2)), int(mReg.group(3))

                if primaryPos[i] == "1" and pairPos != "0":
                    for key in splice2count2:
                        chr_key, start_key, end_key  = key[0], key[1], key[2]

                        # if the given region cover arround the spliced sites
                        if F[2] == "+" and start_reg <= end_key - 5 and end_key + 5 <= end_reg: splice2count2_ref[key] = splice2count2_ref[key] + 1
                        if F[2] == "-" and start_reg <= start_key - 5 and start_key + 5 <= end_reg: splice2count2_ref[key] = splice2count2_ref[key] + 1

                if primaryPos[i] == "2" and pairPos != "0":
                    for key in splice2count1:
                        chr_key, start_key, end_key  = key[0], key[1], key[2]

                        # if the given region cover arround the spliced sites
                        if F[2] == "+" and start_reg <= end_key - 5 and end_key + 5 <= end_reg: splice2count1_ref[key] = splice2count1_ref[key] + 1
                        if F[2] == "-" and start_reg <= start_key - 5 and start_key + 5 <= end_reg: splice2count1_ref[key] = splice2count1_ref[key] + 1


        # check about each splicing position for the pair reads 
        for i in range(0, len(coveredRegion_pair)):
            for reg in coveredRegion_pair[i].split(','):
                if reg == "*": continue
                mReg = regRe.match(reg)
                chr_reg, start_reg, end_reg  = mReg.group(1), int(mReg.group(2)), int(mReg.group(3))

                if (primaryPos[i] == "1" and pairPos[i] == "1") or (primaryPos[i] == "2" and pairPos[i] == "2"):
                    for key in splice2count1:
                        chr_key, start_key, end_key  = key[0], key[1], key[2]

                        # if the given region cover arround the spliced sites
                        if F[2] == "+" and start_reg <= end_key - 5 and end_key + 5 <= end_reg: splice2count1_ref[key] = splice2count1_ref[key] + 1
                        if F[2] == "-" and start_reg <= start_key - 5 and start_key + 5 <= end_reg: splice2count1_ref[key] = splice2count1_ref[key] + 1

                if (primaryPos[i] == "1" and pairPos[i] == "2") or (primaryPos[i] == "2" and pairPos[i] == "1"):
                    for key in splice2count2:
                        chr_key, start_key, end_key  = key[0], key[1], key[2]

                        # if the given region cover arround the spliced sites
                        if F[2] == "+" and start_reg <= end_key - 5 and end_key + 5 <= end_reg: splice2count2_ref[key] = splice2count2_ref[key] + 1
                        if F[2] == "-" and start_reg <= start_key - 5 and start_key + 5 <= end_reg: splice2count2_ref[key] = splice2count2_ref[key] + 1
        ####################


        ####################
        # check whether we should adopt each splice junction
        spliceJunction1 = []
        spliceJunction2 = []
        for key in splice2count1:
            if splice2count1[key] > splice2count1_ref[key]: spliceJunction1.append(key)

        for key in splice2count2:
            if splice2count2[key] > splice2count2_ref[key]: spliceJunction2.append(key)
        ####################


        ####################
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
        ####################


        ####################
        contig1 = []
        contig2 = []

        tempPos = int(F[1])
        if F[2] == "+":
            for key in sorted(spliceJunction1, key=lambda x: x[2], reverse=True):
                if tempPos < key[2]: continue

                minStart = float("inf")
                for reg in region1.regionVec:
                    mReg = regRe.match(reg)
                    chr_reg, start_reg, end_reg  = mReg.group(1), int(mReg.group(2)), int(mReg.group(3))
                    if tempPos <= end_reg and start_reg < minStart:
                        minStart = start_reg

                if minStart != float("inf"):
                    contig1.append((chr_reg, max(key[2], minStart), tempPos))
                    if minStart > key[2]: break 
                    tempPos = key[1]

            minStart = float("inf")
            for reg in region1.regionVec:
                mReg = regRe.match(reg)
                chr_reg, start_reg, end_reg  = mReg.group(1), int(mReg.group(2)), int(mReg.group(3))
                if tempPos <= end_reg and start_reg < minStart:
                    minStart = start_reg

            if minStart != float("inf"):
                contig1.append((chr_reg, minStart, tempPos))

        if F[2] == "-":
            for key in sorted(spliceJunction1, key=lambda x: x[1]):
                if tempPos > key[1]: continue

                maxEnd = - float("inf")            
                for reg in region1.regionVec:
                    mReg = regRe.match(reg)
                    chr_reg, start_reg, end_reg  = mReg.group(1), int(mReg.group(2)), int(mReg.group(3))
                    if start_reg <= tempPos and end_reg > maxEnd:
                        maxEnd = end_reg
       
                if maxEnd != - float("inf"):
                    contig1.append((chr_reg, tempPos, min(key[1], maxEnd)))
                    if maxEnd < key[1]: break 
                    tempPos = key[2]

            maxEnd = - float("inf")
            for reg in region1.regionVec:
                mReg = regRe.match(reg)
                chr_reg, start_reg, end_reg  = mReg.group(1), int(mReg.group(2)), int(mReg.group(3))
                if start_reg <= tempPos and end_reg > maxEnd:
                    maxEnd = end_reg 

            if maxEnd != - float("inf"):
                contig1.append((chr_reg, tempPos, maxEnd))


        tempPos = int(F[4])
        if F[5] == "+":
            for key in sorted(spliceJunction2, key=lambda x: x[2], reverse=True):
                if tempPos < key[2]: continue

                minStart = float("inf")
                for reg in region2.regionVec:
                    mReg = regRe.match(reg)
                    chr_reg, start_reg, end_reg  = mReg.group(1), int(mReg.group(2)), int(mReg.group(3))
                    if tempPos <= end_reg and start_reg < minStart:
                        minStart = start_reg

                if minStart != float("inf"):
                    contig2.append((chr_reg, max(key[2], minStart), tempPos))
                    if minStart > key[2]: break
                    tempPos = key[1]

            minStart = float("inf")
            for reg in region2.regionVec:
                mReg = regRe.match(reg)
                chr_reg, start_reg, end_reg  = mReg.group(1), int(mReg.group(2)), int(mReg.group(3))
                if tempPos <= end_reg and start_reg < minStart:
                    minStart = start_reg

            if minStart != float("inf"):
                contig2.append((chr_reg, minStart, tempPos))

        if F[5] == "-":
            for key in sorted(spliceJunction2, key=lambda x: x[1]):
                if tempPos > key[1]: continue

                maxEnd = - float("inf")
                for reg in region2.regionVec:
                    mReg = regRe.match(reg)
                    chr_reg, start_reg, end_reg  = mReg.group(1), int(mReg.group(2)), int(mReg.group(3))
                    if start_reg <= tempPos and end_reg > maxEnd:
                        maxEnd = end_reg
      
                if maxEnd != - float("inf"):
                    contig2.append((chr_reg, tempPos, min(key[1], maxEnd)))
                    if maxEnd < key[1]: break
                    tempPos = key[2]

            maxEnd = - float("inf")
            for reg in region2.regionVec:
                mReg = regRe.match(reg)
                chr_reg, start_reg, end_reg  = mReg.group(1), int(mReg.group(2)), int(mReg.group(3))
                if start_reg <= tempPos and end_reg > maxEnd:
                    maxEnd = end_reg

            if maxEnd != - float("inf"):
                contig2.append((chr_reg, tempPos, maxEnd))
        ####################

        contigSeq1 = ""
        contigSeq2 = ""

        inSeq = ';'.join(list(set(F[6].split(';'))))

        contigSeq1 = seq_utils.getSeq(reference_genome, contig1)
        contigSeq2 = seq_utils.getSeq(reference_genome, contig2)
        if F[2] == "+": contigSeq1 = seq_utils.reverseComplement(contigSeq1)
        if F[5] == "+": contigSeq2 = seq_utils.reverseComplement(contigSeq2)

        print('\t'.join(['\t'.join(F[0:6]), str(inSeq), str(len(F[8].split(';'))), str(contig1), str(contig2), contigSeq1, contigSeq2]), file = hOUT)


    hIN.close()
    hOUT.close()



def makeJucSeqPairFa(inputFilePath, outputFilePath):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    for line in hIN:
        F = line.rstrip('\n').split('\t')

        print('>' + F[0] + ":" + F[2] + F[1] + "-" + F[3] + ":" + F[5] + F[4] + "_contig" + '1', file = hOUT)
        print(F[10], file = hOUT)

        print('>' + F[0] + ":" + F[2] + F[1] + "-" + F[3] + ":" + F[5] + F[4] + "_contig" + '2', file = hOUT)
        print(F[11], file = hOUT)

    hIN.close()
    hOUT.close()



def checkMatching(inputFilePath, outputFilePath):

    # min_allowed_contig_match_diff = config.param_conf.getint("filter_condition", "min_allowed_contig_match_diff")
    # check_contig_size_other_breakpoint = config.param_conf.getint("filter_condition", "check_contig_size_other_breakpoint")
    min_allowed_contig_match_diff = param_conf.min_allowed_contig_match_diff
    check_contig_size_other_breakpoint = param_conf.check_contig_size_other_breakpoint

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    
    tempID = ""
    for line in hIN:
        F = line.rstrip('\n').split('\t')

        if re.match(r'^\d', F[0]) is None: continue

        if F[9] != tempID:
            if tempID != "":
                otherMatch = []
                for site, value in sorted(site2Match.items(), key = lambda x: x[1], reverse=True):
                    if site == targetAln: continue
                    if value >= targetScore - min_allowed_contig_match_diff: otherMatch.append(site)
                    if len(otherMatch) >= 10: break
                otherMatch_str = ("---" if len(otherMatch) == 0 else ';'.join(otherMatch))
                print(tempID + '\t' + targetAln + '\t' + otherMatch_str + '\t' + str(crossMatch) + '\t' + str(targetScore) + '\t' + str(targetSize), file = hOUT)

            tempID = F[9]
            targetAln = "---"
            targetScore= 0
            site2Match = {}
            crossMatch = 0
            targetSize = "---"

        site2Match[F[13] + ':' + str(int(F[15]) + 1) + '-' + F[16]] = int(F[0])

        chr1, strand1, pos1, chr2, strand2, pos2, contigNum = "", "", "", "", "", "", 0
        contigMatch = ReContig.match(F[9])
        if contigMatch is None:
            print("the format of the fasta name is inconsistent at", file = sys.stderr)
            print('\t'.join(F), file = sys.stderr)

        contigNum = contigMatch.group(7)
        if contigNum == "1":
            chr1, strand1, pos1, chr2, strand2, pos2 = contigMatch.group(1), contigMatch.group(2), int(contigMatch.group(3)), contigMatch.group(4), contigMatch.group(5), int(contigMatch.group(6))
        elif contigNum == "2":
            chr1, strand1, pos1, chr2, strand2, pos2 = contigMatch.group(4), contigMatch.group(5), int(contigMatch.group(6)), contigMatch.group(1), contigMatch.group(2), int(contigMatch.group(3))

        if chr1 == F[13]:
            if (strand1 == "+" and F[8] == "-" and abs(pos1 - int(F[16])) < 5) or (strand1 == "-" and F[8] == "+" and abs(pos1 - 1 - int(F[15])) < 5):

                targetAln = F[13] + ':' + str(int(F[15]) + 1) + '-' + F[16] 
                targetScore= int(F[0])
                targetSize = int(F[10])

        if chr2 == F[13] and (pos2 >= int(F[15]) - check_contig_size_other_breakpoint and pos2 <= int(F[16]) + check_contig_size_other_breakpoint):
            crossMatch = float(F[0]) / float(F[10])


    hIN.close()

    # last procedure    
    if tempID != "":
        otherMatch = []
        for site, value in sorted(site2Match.items(), key = lambda x: x[1], reverse=True):
            if site == targetAln: continue
            if value >= targetScore - min_allowed_contig_match_diff: otherMatch.append(site)
            if len(otherMatch) >= 10: break
        otherMatch_str = ("---" if len(otherMatch) == 0 else ';'.join(otherMatch))
        print(tempID + '\t' + targetAln + '\t' + otherMatch_str + '\t' + str(crossMatch) + '\t' + str(targetScore) + '\t' + str(targetSize), file = hOUT)

    hOUT.close()



def filterContigCheck(inputFilePath, outputFilePath, checkMatchFile):

    # min_cover_size = config.param_conf.getint("filter_condition", "min_cover_size")
    min_cover_size = param_conf.min_cover_size

    hIN = open(checkMatchFile, 'r')
    key2match1 = {}
    key2match2 = {}
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        keyMatch = ReContig.match(F[0])

        chr1, strand1, pos1, chr2, strand2, pos2, contigNum = keyMatch.group(1), keyMatch.group(2), keyMatch.group(3), keyMatch.group(4), keyMatch.group(5), keyMatch.group(6), keyMatch.group(7)

        if contigNum == "1":
            key2match1['\t'.join([chr1, pos1, strand1, chr2, pos2, strand2])] = '\t'.join(F[1:6])

        if contigNum == "2":
            key2match2['\t'.join([chr1, pos1, strand1, chr2, pos2, strand2])] = '\t'.join(F[1:6])

    hIN.close()


    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    for line in hIN:
        F = line.rstrip('\n').split('\t')

        key = '\t'.join(F[0:6])

        match1 = (key2match1[key] if key in key2match1 else "---\t---\t---\t---\t---")
        match2 = (key2match2[key] if key in key2match2 else "---\t---\t---\t---\t---") 

        matches1 = match1.split('\t')
        matches2 = match2.split('\t')

        if matches1[0] == "---" or matches2[0] == "---": continue
        if matches1[1] != "---" or matches2[1] != "---": continue
        if float(matches1[2]) > 0 or float(matches2[2]) > 0: continue
        # if int(matches1[3]) < min_cover_size or int(matches2[3]) < min_cover_size: continue

        print(key + '\t' + F[6] + '\t' + F[7] + '\t' +  match1 + '\t' + match2, file = hOUT)

    hIN.close()
    hOUT.close()

