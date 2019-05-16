#! /usr/bin/env python

from .run import *
import argparse

def create_parser():

    parser = argparse.ArgumentParser(prog = "fusionfusion")

    parser.add_argument("--version", action = "version", version = "fusionfusion-0.5.0b1")

    parser.add_argument("--star", metavar = "star.Chimeric.out.sam", default = None, type = str,
                        help = "the path to the chimeric sam file by STAR")

    parser.add_argument("--ms2", metavar = "ms2.bam", default = None, type = str,
                        help = "the path to the bam file by Map splice2")

    parser.add_argument("--th2", metavar = "th2.bam", default = None, type = str,
                        help = "the path to the bam file by TopHat2")

    parser.add_argument("--out", metavar = "output_dir", default = None, type = str, required=True,
                        help = "the path to the output directory")

    parser.add_argument("--reference_genome", metavar = "reference.fa", default = None, type = str, required=True,
                        help = "reference genome used for creating validation sequences")

    parser.add_argument("--grc", default = False, action = 'store_true',
                        help = "convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")

    parser.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                        help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    # parser.add_argument("--resource_dir", metavar = "resource_dir", type = str, required=True,
    #                     help = "annotation information directory")

    parser.add_argument("--pooled_control_file", default = None, type = str,
                        help = "the path to control data created by merge_control (default: %(default)s)")

    parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")


    parse_junction_group = parser.add_argument_group("parse_junction_condition",
                                                           "parameters used for parsing breakpoint containing read pairs from bam files")

    parse_junction_group.add_argument("--abnormal_insert_size", type = int, default = 500000,
                                      help = "size of abnormal insert size. used for checking the consistency of paired read of breakpoint containing reads (default: %(default)s)")

    parse_junction_group.add_argument("--min_major_clipping_size", type = int, default = 15,
                                      help = "minimum number of clipped bases for junction read (default: %(default)s)")


    filter_condition_group = parser.add_argument_group("filter_condition",
                                                            "parameters used in various filtering steps")

    filter_condition_group.add_argument("--min_read_pair_num", type = int, default = 3,
                                       help = "minimum required number of supporting junction read pairs (default: %(default)s)")

    filter_condition_group.add_argument("--min_valid_read_pair_ratio", type = float, default = 0.8,
                                       help = "minimum ratio of proper supporting junction read pairs (used for map-splice2) (default: %(default)s)")

    filter_condition_group.add_argument("--min_cover_size", type = int, default = 30,
                                        help = "region size which have to be covered by aligned short reads (default: %(default)s)")

    filter_condition_group.add_argument("--anchor_size_thres", type = int, default = 10,
                                       help = "at least an anchor size of one chimeric read have to be equal or larger than the specified value (default: %(default)s)")
     
    filter_condition_group.add_argument("--min_chimeric_size", type = int, default = 1000,
                                        help = "threshold of minimum chimeric transcript sizes (default: %(default)s)")

    filter_condition_group.add_argument("--min_allowed_contig_match_diff", type = int, default = 3,
                                        help = "if contigs are aligned on other positions with less than the specified mismatch value, \
                                                then the corresponding fusion transcipts are filtered (default: %(default)s)")
     
    filter_condition_group.add_argument("--check_contig_size_other_breakpoint", type = int, default = 300,
                                        help = "if contigs are aligned on within the specified size from the other breakpoint, \
                                                then the corresponding fusion transcripts are filtered (default: %(default)s)")

    filter_condition_group.add_argument("--filter_same_gene", default = False, action = 'store_true', 
                                        help = "if two breakpoints are on the same gene, then the fusion candidate is filtered %(default)s)")


    return parser

