
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

        self._pl_bw_path = os.path.join("RGTools", "large_files", "ENCFF565BWR.pl.bw")
        self._mn_bw_path = os.path.join("RGTools", "large_files", "ENCFF775FNU.mn.bw")

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
        args.negative_mn = False
        args.flip_mn = True

        return args

    def test_count_paired_bw(self):
        args = self.get_count_paired_bw_simple_args()

        CountPairedBw.main(args)

        output = np.load(args.opath)

        self.assertEqual(output.shape, (3, 1))
        self.assertEqual(output[0, 0], 123)
        self.assertEqual(output[2, 0], 1877)

        args.quantification_type = "full_track"
        CountPairedBw.main(args)

        output = np.load(args.opath)

        self.assertEqual(output.shape, (3, 1001))
        self.assertTrue(np.sum(output[0, :]) > 0)

        args.negative_mn = True
        CountPairedBw.main(args)

        output = np.load(args.opath)

        self.assertEqual(output.shape, (3, 1001))
        self.assertTrue(np.sum(output[0, :]) < 0)

    def test_count_single_bw_npz(self):
        args = self.get_count_single_bw_simple_args()
        args.opath = os.path.join(self._test_path, "output.npz")
        
        CountSingleBw.main(args)
        
        # Load from .npz
        data = np.load(args.opath)
        self.assertIn("arr_0", data.files)
        output = data["arr_0"]
        self.assertEqual(output.shape, (3, 1))

    def test_quantification_types(self):
        # Test RPK in single bw
        args = self.get_count_single_bw_simple_args()
        args.quantification_type = "RPK"
        CountSingleBw.main(args)
        output_rpk = np.load(args.opath)
        
        args.quantification_type = "raw_count"
        CountSingleBw.main(args)
        output_raw = np.load(args.opath)
        
        # RPK = raw / len * 1000. Regions are 1001bp.
        expected_rpk = output_raw / 1001 * 1000
        np.testing.assert_allclose(output_rpk, expected_rpk, rtol=1e-5)

    def test_heterogeneous_lengths_padding(self):
        # Create a bed file with different lengths
        hetero_bed = os.path.join(self._test_path, "hetero.bed")
        with open(hetero_bed, "w") as f:
            f.write("chr14\t75278325\t75279326\n") # 1001bp
            f.write("chr17\t45894026\t45894526\n") # 500bp
        
        args = self.get_count_single_bw_simple_args()
        args.region_file_path = hetero_bed
        args.quantification_type = "full_track"
        
        CountSingleBw.main(args)
        output = np.load(args.opath)
        
        # Max length should be 1001
        self.assertEqual(output.shape, (2, 1001))
        # Second region (500bp) should be zero-padded from index 500 onwards
        self.assertTrue(np.all(output[1, 500:] == 0))
        self.assertFalse(np.all(output[1, :500] == 0))

    def test_flip_mn(self):
        args = self.get_count_paired_bw_simple_args()
        args.quantification_type = "full_track"
        args.negative_mn = False
        
        # Test without flip
        args.flip_mn = False
        CountPairedBw.main(args)
        output_no_flip = np.load(args.opath)
        
        # Test with flip
        args.flip_mn = True
        CountPairedBw.main(args)
        output_flip = np.load(args.opath)
        
        # Region 0 is minus strand, so output_flip[1] should be flip(output_no_flip[1])
        np.testing.assert_array_equal(output_flip[0], np.flip(output_no_flip[0]))
        # Region 1 is plus strand, so no flip should occur
        np.testing.assert_array_equal(output_flip[1], output_no_flip[1])

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

        self.assertEqual(output.shape, (3, 1))
        self.assertEqual(output[0, 0], -123)
        self.assertEqual(output[2, 0], -216)

        args.quantification_type = "full_track"
        CountSingleBw.main(args)

        output = np.load(args.opath)

        self.assertEqual(output.shape, (3, 1001))
