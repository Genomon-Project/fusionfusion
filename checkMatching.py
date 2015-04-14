#! /usr/local/bin/python

import sys, re

inputFile = sys.argv[1]

hIN = open(inputFile, 'r')
ReContig = re.compile(r'(\w+):([\+\-])(\d+)\-(\w+):([\+\-])(\d+)_contig([12])')

tempID = ""
for line in hIN:
    F = line.rstrip('\n').split('\t')

    if re.match(r'^\d', F[0]) is None: continue

    if F[9] != tempID:
        if tempID != "":
            otherMatch = []
            for site, value in sorted(site2Match.items(), key = lambda x: x[1], reverse=True):
                if site == targetAln: continue
                if value >= targetScore - 3: otherMatch.append(site)
                if len(otherMatch) >= 10: break
            otherMatch_str = ("---" if len(otherMatch) == 0 else ';'.join(otherMatch))
            print tempID + '\t' + targetAln + '\t' + otherMatch_str + '\t' + str(crossMatch)

        tempID = F[9]
        targetAln = "---"
        targetScore= 0
        site2Match = {}
        crossMatch = 0

    site2Match[F[13] + ':' + str(int(F[15]) + 1) + '-' + F[16]] = int(F[0])

    chr1, strand1, pos1, chr2, strand2, pos2, contigNum = "", "", "", "", "", "", 0
    contigMatch = ReContig.match(F[9])
    if contigMatch is None:
        print >> sys.stderr, "the format of the fasta name is inconsistent at"
        print >> sys.stderr, '\t'.join(F)

    contigNum = contigMatch.group(7)
    if contigNum == "1":
        chr1, strand1, pos1, chr2, strand2, pos2 = contigMatch.group(1), contigMatch.group(2), int(contigMatch.group(3)), contigMatch.group(4), contigMatch.group(5), int(contigMatch.group(6))
    elif contigNum == "2":
        chr1, strand1, pos1, chr2, strand2, pos2 = contigMatch.group(4), contigMatch.group(5), int(contigMatch.group(6)), contigMatch.group(1), contigMatch.group(2), int(contigMatch.group(3))

    if chr1 == F[13]:
        if (strand1 == "+" and F[8] == "-" and abs(pos1 - int(F[16])) < 5) or (strand1 == "-" and F[8] == "+" and abs(pos1 - 1 - int(F[15]))) < 5:
            targetAln = F[13] + ':' + str(int(F[15]) + 1) + '-' + F[16] 
            targetScore= int(F[0])

    if chr2 == F[13] and (pos2 >= int(F[15]) - 1000 and pos2 <= int(F[16]) + 1000):
        crossMatch = float(F[0]) / float(F[10])


hIN.close()


otherMatch = []
for site, value in sorted(site2Match.items(), key = lambda x: x[1], reverse=True):
    if site == targetAln: continue
    if value >= targetScore - 3: otherMatch.append(site)
    if len(otherMatch) >= 10: break
otherMatch_str = ("---" if len(otherMatch) == 0 else ';'.join(otherMatch))
print tempID + '\t' + targetAln + '\t' + otherMatch_str + '\t' + str(crossMatch)


