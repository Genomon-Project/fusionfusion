import collections
import copy
import itertools
import pysam

from . import annotationFunction

class ShortRangeChimeraFilter:
    def __init__(self, ref_gene_tb, ens_gene_tb):
        self.ref_gene_tb = ref_gene_tb
        self.ens_gene_tb = ens_gene_tb
        self.interval = 500

    def run(self, splicing_file, bam_file, output_file):
        with open(splicing_file) as splicing_reader, \
             pysam.AlignmentFile(bam_file, 'r') as bam_reader, \
             pysam.AlignmentFile(output_file, 'w', header=bam_reader.header) as sam_writer:
            for chr, start, end in self.load_splicing_junctions(splicing_reader):
                pair_alns = collections.defaultdict(set)
                for aln in self.load_alignments(bam_reader, chr, start, end):
                    locus = (aln.reference_name, aln.reference_start)
                    entry = (aln.query_name, aln.is_read1)
                    if entry in pair_alns[locus]:
                        self.write_alignment(sam_writer, aln, chr, start, end)
                        pair_alns[locus].remove(entry)
                    else:
                        aln1, aln2 = self.split_alignment(aln, start)
                        self.write_alignment(sam_writer, aln1, chr, start, end)
                        self.write_alignment(sam_writer, aln2, chr, start, end)
                        next_locus = (aln.next_reference_name, aln.next_reference_start)
                        next_entry = (aln.query_name, not aln.is_read1)
                        pair_alns[next_locus].add(next_entry)
                # fetch and write pair reads
                for aln in self.load_pair_alignments(bam_reader, pair_alns):
                    self.write_alignment(sam_writer, aln, chr, start, end)

    def write_alignment(self, sam_writer, aln, chr, start, end):
        # Strip tags since they are unnecessary for succeeding processes
        aln.set_tags([])
        aln.query_name += '_{}:{}-{}'.format(chr, start, end+1)
        sam_writer.write(aln)

    def get_genes_at(self, chr, pos):
        return set(annotationFunction.get_gene_info(chr, pos, self.ref_gene_tb, self.ens_gene_tb))

    def load_splicing_junctions(self, splicing_reader):
        for line in splicing_reader:
            chr, start, end = line.split('\t')[:3]
            start_genes = self.get_genes_at(chr, start)
            end_genes = self.get_genes_at(chr, end)
            if not(start_genes.intersection(end_genes)):
                yield (chr, int(start)-1, int(end))

    def load_alignments(self, bam_reader, chr, start, end):
        def pairwise(it):
            it1, it2 = itertools.tee(it)
            next(it2, None)
            return zip(it1, it2)

        for aln in bam_reader.fetch(chr, start, start+1):
            if not aln.is_proper_pair or aln.is_secondary: continue

            for (_, curr_end), (next_start, _) in pairwise(aln.get_blocks()):
                if curr_end == start and next_start == end:
                    break
            else:
                continue
            yield aln

    def load_pair_alignments(self, bam_reader, pair_alns):
        def chunked(alns, interval):
            chrom = start = end = None
            entries = set()
            for (chr, pos), es in sorted(alns.items(), key=lambda x: x[0]):
                if chrom == chr and pos <= start + interval:
                    end = pos
                    entries.update(es)
                    continue
                if entries:
                    yield chrom, start, end, entries
                chrom = chr
                start = end = pos
                entries = set(es)
            if entries:
                yield chrom, start, end, entries

        for (chr, start, end, entries) in chunked(pair_alns, self.interval):
            for aln in bam_reader.fetch(chr, start, end + 1):
                locus = (aln.reference_name, aln.reference_start)
                entry = (aln.query_name, aln.is_read1)
                if not aln.is_secondary and entry in entries:
                    yield aln

    def split_alignment(self, aln, splicing_start):
        elems = aln.cigartuples
        blocks = aln.get_blocks()

        elem_idx = 0
        block_idx = 0
        for op, len in elems:
            if op == 0:
                _, block_end = blocks[block_idx]
                if block_end == splicing_start: break
                block_idx += 1
            elem_idx += 1

        elems1 = elems[:elem_idx+1]
        elems2 = list(itertools.dropwhile(lambda elem: elem[0] != 0, elems[elem_idx+1:]))
        count_bases = lambda ops, elems: sum(e[1] for e in elems if e[0] in ops)

        aln1 = copy.copy(aln)
        aln2 = copy.copy(aln)
        prim_flag, suppl_flag = aln.flag & ~0x100, aln.flag | 0x100
        # Consider alignment with more matched bases as primary and the other as supplementary
        if count_bases({0}, elems1) >= count_bases({0}, elems2):
            aln1.flag, aln2.flag = prim_flag, suppl_flag
        else:
            aln1.flag, aln2.flag = suppl_flag, prim_flag
        aln2.reference_start = blocks[block_idx+1][0]
        # To reset CIGAR (and related properties), .cigartuples must be set to None first
        aln1.cigartuples = None; aln1.cigartuples = elems1 + [(4, count_bases({0, 1, 4}, elems2))]
        aln2.cigartuples = None; aln2.cigartuples = [(4, count_bases({0, 1, 4}, elems1))] + elems2
        aln1.template_length, aln2.template_length = \
            (aln.template_length, 0) if aln.template_length >= 0 else (0, aln.template_length)
        return (aln1, aln2)
