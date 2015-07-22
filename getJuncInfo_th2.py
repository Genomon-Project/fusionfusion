#! /usr/local/bin/python

import sys, re, myCigar

inputFile = sys.argv[1]

abnormal_insert_size = 500000
fusionCIGAR = re.compile('([\dMmNnIiDd]+[Mm])(\d+F)([\dMmNnIiDd]+[Mm])')


hIN = open(inputFile, 'r')
for line in hIN:
    F = line.rstrip('\n').split('\t')

    fusionFlag = 0
    XFInfo = ""
    for i in range(11, len(F)):
        if F[i].startswith("XF:Z:1"):
            fusionFlag = 1
            XFInfo = F[i] 
            break

    if fusionFlag == 0: continue

    # be reminded necessary variables
    chr_primary, pos_primary, dir_primary, chr_pair, pos_pair, dir_pair, chr_SA, pos_SA, dir_SA = "", "", "", "", "", "", "", "", ""
    mq_primary, coverRegion_primary, mq_pair, coverRegion_pair, mq_SA, coverRegion_SA = "", "", "", "", "", ""
    juncChr_primary, juncPos_primary, juncDir_primary, juncChr_SA, juncPos_SA, juncDir_SA = "", "", "", "", "", "" 

    # about the samflag
    flags = format(int(F[1]), "#014b")[:1:-1]
    if flags[8] == "1" or flags[11] == "1": continue # pass non-primary alignment
    if flags[10] == "1": continue # pass duplicate reads
    if flags[2] == "1" or flags[3] == "1": continue # pass unless both the pair are aligned

    # search for XP tag
    XPInfo = ""
    for i in range(11, len(F)):
        if F[i].startswith("XP:Z:"):
            XPInfo = F[i]
            break


    # for the primary alignment
    chr_primary = F[2]
    pos_primary = F[3]
    dir_primary = ("-" if flags[4] == "1" else "+")
    mq_primary = F[4]
    cigar_primary = F[5]
    coverRegion_primary = myCigar.getCoverRegion(chr_primary, pos_primary, cigar_primary)

    # for the pair
    chr_pair = F[6]
    pos_pair = F[7]
    dir_pair = ("-" if flags[5] == "1" else "+")
    mq_pair = "0" # since we cannot get the mapping quality of pair reads directly, this is just a dummy
    if len(XPInfo.split(" ")) < 3:
        pass
    cigar_pair = XPInfo.split(" ")[2]
    coverRegion_pair = myCigar.getCoverRegion(chr_pair, pos_pair, cigar_pair)


    # organize fusion infomation
    FF = XFInfo.split(" ")
    juncChr_primary, juncChr_SA = FF[1].split("-")
    juncPos_primary = FF[2]

    FCIGARMatch = fusionCIGAR.match(FF[3])
    juncPos_SA = FCIGARMatch.group(2).rstrip('F')
    cigar_primary = FCIGARMatch.group(1)
    cigar_SA = FCIGARMatch.group(3)
    juncDir_primary = ("-" if cigar_primary.islower() else "+")
    juncDir_SA = ("+" if cigar_SA.islower() else "+")
    


    # the pair read is aligned at the same chromosome with the primary read
    validFlag = 0
    if chr_pair == juncChr_primary:
        if juncDir_primary == "+" and dir_primary == "-" and dir_pair == "+" and 0 <= juncPos_primary - pos_pair < abnormal_insert_size:
            juncType = 1
            validFlag = 1
        if juncDir_primary == "-" and dir_primary == "+" and dir_pair == "-" and 0 <= pos_pair - juncPos_primary < abnormal_insert_size:
            juncType = 1
            validFlag = 1

    if chr_pair == juncChr_SA:
        if juncDir_SA == "+" and dir_SA == "-" and dir_pair == "+" and 0 <= juncPos_SA - pos_pair < abnormal_insert_size:
            juncType = 1
            validFlag = 1
        if juncDir_SA == "-" and dir_SA == "+" and dir_pair == "-" and 0 <= pos_pair - juncPos_SA < abnormal_insert_size:
            juncType = 2 
            validFlag = 1


    if validFlag == 1:

        # TopHat2 seems not to consider inserted short bases
        juncSurplus = "---"

        # reorder by the chromosome position and print
        if juncChr_primary < juncChr_SA or juncChr_primary == juncChr_SA and juncPos_primary <= juncPos_SA:
            print '\t'.join([juncChr_primary, str(juncPos_primary), juncDir_primary, juncChr_SA, str(juncPos_SA), juncDir_SA, \
                             juncSurplus, readID_primary, mq_primary, coverRegion_primary, dir_primary, \
                             mq_pair, coverRegion_pair, dir_pair, mq_SA, coverRegion_SA, dir_SA, str(juncType) , "1"])

        else:
            print '\t'.join([juncChr_SA, str(juncPos_SA), juncDir_SA, juncChr_primary, str(juncPos_primary), juncDir_primary, \
                             juncSurplus, readID_primary, mq_primary, coverRegion_primary, dir_primary, \
                             mq_pair, coverRegion_pair, dir_pair, mq_SA, coverRegion_SA, dir_SA, str(juncType) , "2"])


hIN.close()


