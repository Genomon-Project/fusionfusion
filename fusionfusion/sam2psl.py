class InvalidSAMError(Exception):
    def __init__(self, message):
        super(Exception, self).__init__(message)

class SAMRecord:
    BAM_FPAIRED = 1
    BAM_FPROPER_PAIR = 2
    BAM_FUNMAP = 4
    BAM_FMUNMAP = 8
    BAM_FREVERSE = 16
    BAM_FMREVERSE = 32
    BAM_FREAD1 = 64
    BAM_FREAD2 = 128
    BAM_FSECONDARY = 256
    BAM_FQCFAIL = 512
    BAM_FDUP = 1024

    def __init__(self, fields):
        self.QNAME = fields[0]
        self.FLAG = int(fields[1])
        self.RNAME = fields[2]
        self.POS = int(fields[3])
        self.MAPQ = int(fields[4])
        self.CIGAR = fields[5]
        self.RNEXT = fields[6]
        self.PNEXT = int(fields[7])
        self.TLEN = int(fields[8])
        self.SEQ = fields[9]
        self.QUAL = fields[10]
        self.AS = -len(self.SEQ)
        self.NM = len(self.SEQ)
        self.MD = ''

        optionalField = 0;
        for field in fields[11:]:
            if field.startswith('AS:'):
                self.AS = int(field[field.rfind(':')+1:])
                optionalField += 1
            elif field.startswith('NM:'):
                self.NM = int(field[field.rfind(':')+1:])
                optionalField += 1
            elif field.startswith('MD:'):
                self.MD = field[field.rfind(':')+1:]
                optionalField += 1
            if optionalField == 3:
                break

class CigarOp:
    def __init__(self, op, len):
        self.op = op
        self.len = len
        self.qPos = 0
        self.tPos = 0

class PosCigar:
    BAM_CIGAR_OPS = set('MIDNSHP=XB')

    def __init__(self):
        self.tid = -1
        self.base = -1
        self.pos = -1
        self.qStart = -1
        self.tStart = -1
        self.qEnd = -1
        self.tEnd = -1
        self.qual = -1
        self.anchor = -1
        self.iclip = -1
        self.ops = []
        self.l_qseq = -1

    @staticmethod
    def resolve_cigar_pos(pos, cigar):
        m = PosCigar()

        if cigar == '*':
            return m

        m.base = 0
        m.pos = pos
        m.tStart = pos
        m.l_qseq = 0
        digits = ''
        for c in cigar:
            if '0' <= c <= '9':
                digits += c
            elif c not in PosCigar.BAM_CIGAR_OPS:
                raise InvalidSAMError('unknown cigar op: ' + c)
            else:
                m.ops.append(CigarOp(c, int(digits)))
                digits = ''

        ns = 0
        end = 0
        for k, e in enumerate(m.ops):
            op = e.op
            l = e.len
            e.qPos = end
            if op in set('MIS=X'):
                end += l
            if op == 'S' and l > ns:
                ns = l
                m.iclip = k
            if op in set('MD=X') and m.anchor < 0:
                m.anchor = k

        if m.anchor < 0:
            m.pos = 0
            return m

        end = m.pos
        for e in m.ops[m.anchor:]:
            e.tPos = end
            if e.op in set('MDNS'):
                end += e.len
        end = m.pos
        if m.anchor > 0:
            for e in m.ops[m.anchor-1::-1]:
                if e.op in set('MDNS'):
                    end -= e.len
                e.tPos = end

        m.qStart = m.ops[m.anchor].qPos
        m.tStart = m.ops[m.anchor].tPos
        for e in reversed(m.ops):
            if e.op in set('SH'): continue
            m.qEnd = e.qPos + e.len - 1
            m.tEnd = e.tPos + e.len - 1
            break
        if m.qEnd <= m.qStart or m.tStart != m.pos or m.tEnd <= m.tStart:
            raise InvalidSAMError('cannot resolve cigar: ' + cigar)
        last_op = m.ops[-1]
        m.l_qseq = last_op.qPos + last_op.len * (last_op.op != 'H')
        return m

class PSLRecord:
    def __init__(self):
        self.matches = 0
        self.misMatches = 0
        self.repMatches = 0
        self.nCount = 0
        self.qNumInsert = 0
        self.qBaseInsert = 0
        self.tNumInsert = 0
        self.tBaseInsert = 0
        self.strand = None
        self.qName = None
        self.qSize = 0
        self.qStart = 0
        self.qEnd = 0
        self.tName = None
        self.tSize = 0
        self.tStart = 0
        self.tEnd = 0
        self.blockCount = 0
        self.blockSizes = []
        self.qStarts = []
        self.tStarts = []
        self.qBlocks = []
        self.tBlocks = []

    def __str__(self):
        return '\t'.join([
            str(self.matches),
            str(self.misMatches),
            str(self.repMatches),
            str(self.nCount),
            str(self.qNumInsert),
            str(self.qBaseInsert),
            str(self.tNumInsert),
            str(self.tBaseInsert),
            self.strand or '',
            self.qName or '',
            str(self.qSize),
            str(self.qStart),
            str(self.qEnd),
            self.tName or '',
            str(self.tSize),
            str(self.tStart),
            str(self.tEnd),
            str(self.blockCount),
            ','.join(str(x) for x in self.blockSizes),
            ','.join(str(x) for x in self.qStarts),
            ','.join(str(x) for x in self.tStarts)
        ])

class SAM2PSLConverter:
    def __init__(self, input_file):
        self.input_file = input_file
        self.target_lengths = {}

    def convert_line(self, fields):
        sam = SAMRecord(fields)
        m = PosCigar.resolve_cigar_pos(sam.POS - 1, sam.CIGAR)

        if len(sam.SEQ) != m.l_qseq and m.l_qseq > 0:
            raise InvalidSAMError('SEQ length not calculated correctly')

        psl = PSLRecord()
        psl.misMatches = sam.NM
        psl.strand = "-" if (sam.FLAG & SAMRecord.BAM_FREVERSE) > 0 else "+"
        psl.qName = sam.QNAME
        psl.qSize = len(sam.SEQ)
        psl.qStart = m.qStart
        psl.qEnd = m.qEnd + 1
        psl.tName = sam.RNAME
        psl.tSize = self.target_lengths[sam.RNAME]
        psl.tStart = m.tStart
        psl.tEnd = m.tEnd + 1

        for b in sam.SEQ:
            if b == 'N' or b == 'n':
                psl.nCount += 1
        for e in m.ops:
            if e.op in set('MD=X'):
                psl.matches += e.len
            if e.op == 'I':
                psl.qNumInsert += 1
                psl.qBaseInsert += e.len
            if e.op == 'D':
                psl.tNumInsert += 1
                psl.tBaseInsert += e.len
            if e.op == 'M':
                psl.blockCount += 1
                psl.blockSizes.append(e.len)
                psl.qStarts.append(e.qPos)
                psl.tStarts.append(e.tPos)

        if psl.strand == '-':
            psl.qStart = psl.qSize - psl.qEnd
            psl.qEnd = psl.qSize - psl.qStart
            for i, size in enumerate(psl.blockSizes):
                qEndOnPlus = psl.qStarts[i] + size
                psl.qStarts[i] = psl.qSize - qEndOnPlus
            psl.blockSizes.reverse()
            psl.qStarts.reverse()
            psl.tStarts.reverse()

        return psl

    def convert_lines(self, lines):
        for line in lines:
            if line.startswith('@'):
                if not line.startswith('@SQ'):
                    continue
                cols = line.split()
                if (len(cols) < 3
                    or not cols[1].startswith('SN:')
                    or not cols[2].startswith('LN:')):
                    continue
                self.target_lengths[cols[1][3:]] = int(cols[2][3:])
            else:
                fields = line.split()
                # TODO: reconfirm if it's ok to skip unmapped alignments
                if len(fields) < 11 or fields[2] == '*':
                    continue
                yield self.convert_line(fields)

    def write(self, output_file):
        with open(self.input_file) as r, open(output_file, 'w') as w:
            w.write(
                'matches\tmisMatches\trepMatches\tnCount\tqNumInsert\t'
                'qBaseInsert\ttNumInsert\ttBaseInsert\tstrand\tqName\t'
                'qSize\tqStart\tqEnd\ttName\ttSize\tStart\ttEnd\t'
                'blockCount\tblockSizes\tqStarts\ttStarts\n'
            )
            for record in self.convert_lines(iter(r)):
                w.write(str(record))
                w.write('\n')

def convert(input_sam_file, output_psl_file):
    conv = SAM2PSLConverter(input_sam_file)
    conv.write(output_psl_file)
