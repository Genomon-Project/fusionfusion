#! /usr/bin/env python

import argparse, yaml
import parseJunctionInfo
import utils

def cluster_filter_junction(inputFilePath, outputFilePrefix):

    parseJunctionInfo.clusterJuncInfo(inputFilePath, outputFilePrefix + ".chimeric.clustered.bedpe")

    


def main(args):

    starBamFile = args.star
    ms2BamFile = args.ms2
    paramInfoFile = args.paramInfoFile
    output_dir =args.out

    try:
        with open(args.paramInfoFile, 'r') as fIN:
            paramConf = yaml.load(fIN)
    except yaml.YAMLError, exc:
        print "Error in sample information file:", exc


    ####################
    # make direcotry
    utils.make_directory(output_dir)
    ####################

    ####################
    # parsing chimeric reads from bam files
    if starBamFile is not None:

        parseJunctionInfo.parseJuncInfo_STAR(starBamFile, output_dir + "/star.chimeric.tmp.bedpe", paramConf)

        utils.sortBedpe(output_dir + "/star.chimeric.tmp.bedpe", output_dir + "/star.chimeric.bedpe")

        cluster_filter_junction(output_dir + "/star.chimeric.bedpe", output_dir + "/star")


    if ms2BamFile is not None:

        parseJunctionInfo.parseJuncInfo_STAR(ms2BamFile, output_dir + "/ms2.chimeric.tmp.bedpe", paramConf) 

        utils.sortBedpe(output_dir + "/ms2.chimeric.tmp.bedpe", output_dir + "/ms2.chimeric.bedpe")

        cluster_filter_junction(output_dir + "/ms2.chimeric.bedpe", output_dir + "/ms2")


