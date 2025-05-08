
import unittest
import argparse
import shutil
import os

import numpy as np

from count_bw import CountSingleBw, CountPairedBw

class CountBwTest(unittest.TestCase):
    def setUp(self):
        self._test_path = "count_bw_test_dir"
        if not os.path.exists(self._test_path):
            os.makedirs(self._test_path)

        self._pl_bw_path = os.path.join("large_data", "ENCFF565BWR.pl.bw")
        self._mn_bw_path = os.path.join("large_data", "ENCFF775FNU.mn.bw")

        self._bed3_path = os.path.join("example_data", "three_genes.bed3")
        self._bed6_path = os.path.join("example_data", "three_genes.bed6")
        self._bed6gene_path = os.path.join("example_data", "three_genes.bed6gene")

    def tearDown(self):
        if os.path.exists(self._test_path):
            shutil.rmtree(self._test_path)
    
    def get_count_paired_bw_simple_args(self):
        args = argparse.Namespace()
        args.subcommand = "count_paired_bw"
        args.bw_pl = self._pl_bw_path
        args.bw_mn = self._mn_bw_path
        args.region_file_path = self._bed6_path
        args.region_file_type = "bed6"
        args.override_strand = None
        args.quantification_type = "raw_count"
        args.opath = os.path.join(self._test_path, "output.npy")

        return args

    def test_count_paired_bw(self):
        args = self.get_count_paired_bw_simple_args()

        CountPairedBw.main(args)

        output = np.load(args.opath)

        self.assertEqual(output.shape, (3,))
        self.assertEqual(output[0], 123)
        self.assertEqual(output[2], 1877)

        args.quantification_type = "full_track"
        CountPairedBw.main(args)

        output = np.load(args.opath)

        self.assertEqual(output.shape, (3, 1001))

    def get_count_single_bw_simple_args(self):
        args = argparse.Namespace()
        args.subcommand = "count_single_bw"
        args.bw_path = self._mn_bw_path
        args.region_file_path = self._bed3_path
        args.region_file_type = "bed3"
        args.quantification_type = "raw_count"
        args.opath = os.path.join(self._test_path, "output.npy")

        return args

    def test_count_single_bw(self):
        args = self.get_count_single_bw_simple_args()

        CountSingleBw.main(args)

        output = np.load(args.opath)

        self.assertEqual(output.shape, (3,))
        self.assertEqual(output[0], -123)
        self.assertEqual(output[2], -216)

        args.quantification_type = "full_track"
        CountSingleBw.main(args)

        output = np.load(args.opath)

        self.assertEqual(output.shape, (3, 1001))
