#! /usr/bin/env python

import re
import regions, seq_utils

regRe = re.compile(r'(\w+):(\d+)\-(\d+)')


def filterCoverRegion(inputFilePath, outputFilePath, Params):

    min_read_pair_num = Params["min_read_pair_num"]
    min_valid_read_pair_ratio = Params["min_valid_read_pair_ratio"]
    min_cover_size = Params["min_cover_size"]
    min_chimeric_size = Params["min_chimeric_size"]

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    for line in hIN:
        F = line.rstrip('\n').split('\t')

        if F[0] == F[3] and abs(int(F[1]) - int(F[4])) < min_chimeric_size: continue

        if F[1] == "4489051":
            pass

        # filter by the number of supporting read pairs
        # IDs = F[7].split(';')
        # uIDs = list(set(IDs))
        # if len(uIDs) < min_read_pair_num: continue
        coveredRegion_primary = F[9].split(';')
        coveredRegion_pair = F[12].split(';')
        coveredRegion_SA = F[15].split(';')
        coveredRegion_meta = [coveredRegion_primary[i] + ';' + coveredRegion_pair[i] + ';' + coveredRegion_SA[i] for i in range(0, len(coveredRegion_primary))]
        uniqueCoverdRegion_meta = list(set(coveredRegion_meta))
        if len(uniqueCoverdRegion_meta) < min_read_pair_num: continue

        # filter by the ratio of valid read pairs
        pairPos = F[17].split(';')
        if float(pairPos.count("0")) / float(len(pairPos)) > 1 - min_valid_read_pair_ratio: continue

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

        if region1.regionSize() < min_cover_size or region2.regionSize() < min_cover_size: continue

        print >> hOUT, '\t'.join(F)

    hIN.close()
    hOUT.close()


def extractSplicingPattern(inputFilePath, outputFilePath, Params):

    reference_genome = Params['reference_genome']    
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

        print >> hOUT, '\t'.join(['\t'.join(F[0:6]), str(inSeq), str(len(F[8].split(';'))), str(contig1), str(contig2), contigSeq1, contigSeq2])


    hIN.close()
    hOUT.close()



def makeJucSeqPairFa(inputFilePath, outputFilePath, Params):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    for line in hIN:
        F = line.rstrip('\n').split('\t')

        print >> hOUT, '>' + F[0] + ":" + F[2] + F[1] + "-" + F[3] + ":" + F[5] + F[4] + "_contig" + '1'
        print >> hOUT, F[10]

        print >> hOUT, '>' + F[0] + ":" + F[2] + F[1] + "-" + F[3] + ":" + F[5] + F[4] + "_contig" + '2'
        print >> hOUT, F[11]

    hIN.close()
    hOUT.close()


