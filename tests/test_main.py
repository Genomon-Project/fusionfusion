#! /usr/bin/env python

from __future__ import print_function

import unittest
import os, tempfile, shutil, filecmp
import fusionfusion 
from .check_download import *

class TestMain(unittest.TestCase):

    def setUp(self):
        # prepare reference genome
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        check_download("https://storage.googleapis.com/friend1ws_package_data/common/GRCh37.fa", \
                       cur_dir + "/resource/reference_genome/GRCh37.fa")

        check_download("https://storage.googleapis.com/friend1ws_package_data/fusionfusion/MCF-7.Chimeric.out.sam", \
                       cur_dir + "/resource/star/MCF-7.Chimeric.out.sam")
 
        self.parser = fusionfusion.arg_parser.create_parser()


    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        star_chimeric_sam = cur_dir + "/resource/star/MCF-7.Chimeric.out.sam"
        output_dir = tmp_dir 
        ref_genome = cur_dir + "/resource/reference_genome/GRCh37.fa"
 
        output_file = tmp_dir + "/fusion_fusion.result.txt"
        answer_file = cur_dir + "/data/fusion/MCF-7/fusion_fusion.result.txt"

        print(' '.join(["--star", star_chimeric_sam, "--out", output_dir, "--reference_genome", ref_genome, "--grc"]))
        args = self.parser.parse_args(["--star", star_chimeric_sam, "--out", output_dir, "--reference_genome", ref_genome, "--grc", "--debug"])
        fusionfusion.run.fusionfusion_main(args)

        # self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        self.assertTrue(160 <= len(open(tmp_dir + "/fusion_fusion.result.txt", 'r').readlines()) <= 165)

        # shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()

