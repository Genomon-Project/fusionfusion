#! /usr/bin/env python

# quit using ConfigParser... Y.S. 8/19, 2016
# import ConfigParser
 
# global param_conf

# param_conf = ConfigParser.SafeConfigParser()

global param_conf

class Param_conf(object):

    def __init__(self, reference_genome = None, resource_dir = None, debug = None,
                       abnormal_insert_size = None, min_major_clipping_size = None,
                       min_read_pair_num = None, min_valid_read_pair_ratio = None, min_cover_size = None,
                       anchor_size_thres = None, min_chimeric_size = None, min_allowed_contig_match_diff = None, 
                       check_contig_size_other_breakpoint = None, filter_same_gene = None):

        self.reference_genome = reference_genome
        self.resource_dir = resource_dir
        self.debug = debug
        self.abnormal_insert_size = abnormal_insert_size
        self.min_major_clipping_size = min_major_clipping_size
        self.min_read_pair_num = min_read_pair_num
        self.min_valid_read_pair_ratio = min_valid_read_pair_ratio
        self.min_cover_size = min_cover_size
        self.anchor_size_thres = anchor_size_thres
        self.min_chimeric_size = min_chimeric_size
        self.min_allowed_contig_match_diff = min_allowed_contig_match_diff
        self.check_contig_size_other_breakpoint = check_contig_size_other_breakpoint
        self.filter_same_gene = filter_same_gene

param_conf = Param_conf()
