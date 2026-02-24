
import unittest
import argparse
import shutil
import os

import pandas as pd
import numpy as np

from export import GenomicElementExport
from RGTools.BedTable import BedTable3
from RGTools.GenomicElements import GenomicElements

class GenomicElementExportTest(unittest.TestCase):
    def setUp(self):
        self.__bed3_path = os.path.join("example_data", "three_genes.bed3")
        self.__bed6gene_path = os.path.join("example_data", "three_genes.bed6gene")
        self.__fasta_path = os.path.join("RGTools", "large_files", "hg38.fa")
        self.__wdir = "GenomicElementExport_test"
        if not os.path.exists(self.__wdir):
            os.makedirs(self.__wdir)
        
        self.__sample1_npy_path = os.path.join(self.__wdir, "sample1.npy")
        self.__sample2_npy_path = os.path.join(self.__wdir, "sample2.npy")
        np.save(self.__sample1_npy_path, np.array([1,2,3]))
        np.save(self.__sample2_npy_path, np.array([4,5,6]))

        # Create test signal tracks for TREbed export
        # Regions from three_genes.bed3: 
        # chr14:75278325-75279326 (1001bp), chr17:45894026-45895027 (1001bp), chr6:170553801-170554802 (1001bp)
        # Create tracks with known peak positions
        self.__pl_track_path = os.path.join(self.__wdir, "pl_track.npz")
        self.__mn_track_path = os.path.join(self.__wdir, "mn_track.npz")
        
        # Plus strand tracks: peaks at positions 100, 200, 50 (relative to region start)
        pl_track = np.zeros((3, 1001))
        pl_track[0, 100] = 10.0  # Peak at position 100 for region 0
        pl_track[1, 200] = 10.0  # Peak at position 200 for region 1
        pl_track[2, 50] = 10.0   # Peak at position 50 for region 2
        
        # Minus strand tracks: peaks at positions 300, 150, 250 (relative to region start)
        mn_track = np.zeros((3, 1001))
        mn_track[0, 300] = 10.0  # Peak at position 300 for region 0
        mn_track[1, 150] = 10.0  # Peak at position 150 for region 1
        mn_track[2, 250] = 10.0  # Peak at position 250 for region 2
        
        np.savez(self.__pl_track_path, pl_track)
        np.savez(self.__mn_track_path, mn_track)

        self.__chrom_size_path = os.path.join(self.__wdir, "test.chrom.sizes")
        self.__chrom_size_df = pd.DataFrame({"chrom": ["chr1", "chr2", "chr3", "chr4", "chrFake", "chr6"],
                                             "size": [248956422, 242193529, 198295559, 190214555, 1000000, 170805979],
                                             },
                                            columns=["chrom", "size"],
                                            )
        self.__chrom_size_df.to_csv(self.__chrom_size_path,
                                    sep="\t",
                                    header=False,
                                    index=False,
                                    )

        super().setUp()

    def tearDown(self):
        if os.path.exists(self.__wdir):
            shutil.rmtree(self.__wdir)

        super().tearDown()

    def test_export_exogeneous_sequences(self):
        args = argparse.Namespace(
            region_file_path=self.__bed3_path,
            region_file_type="bed3",
            fasta_path=self.__fasta_path,
            oheader=os.path.join(self.__wdir, "test"),
            oformat="ExogeneousSequences",
        )
        GenomicElementExport.export_exogeneous_sequences(args)

        ofile = os.path.join(self.__wdir, "test.fa")
        with open(ofile, "r") as handle:
            lines = handle.readlines()
            self.assertEqual(len(lines), 6)
            self.assertEqual(lines[1][:10], "CCCCATCCCC")
    
    def test_export_count_table(self):

        args = argparse.Namespace(
            region_file_path=self.__bed3_path,
            region_file_type="bed3",
            opath=os.path.join(self.__wdir, "test.count.csv"),
            oformat="CountTable",
            region_id_type="default",
            sample_name=["sample1", "sample2"],
            stat_npy=[self.__sample1_npy_path, self.__sample2_npy_path],
        )
        GenomicElementExport.export_count_table(args)

        output_df = pd.read_csv(args.opath, 
                                index_col=0,
                                )

        self.assertEqual(output_df.shape, (3, 2))
        self.assertEqual(output_df.iloc[0, 0], 1)
        self.assertEqual(output_df.iloc[1, 0], 2)
        self.assertEqual(output_df.iloc[2, 0], 3)

        # test gene_symbol id type
        args.region_file_path = self.__bed6gene_path
        args.region_file_type = "bed6gene"
        args.region_id_type = "gene_symbol"
        GenomicElementExport.export_count_table(args)

        output_df = pd.read_csv(args.opath, 
                                index_col=0,
                                )

        self.assertEqual(output_df.shape, (3, 2))
        self.assertEqual(output_df.iloc[0, 0], 1)
        self.assertEqual(output_df.iloc[2, 0], 3)
        self.assertEqual(output_df.index[1], "gene3")

    def test_export_chrom_filtered_ge(self):
        args = argparse.Namespace(
            region_file_path=self.__bed3_path,
            region_file_type="bed3",
            chrom_size=self.__chrom_size_path,
            opath=os.path.join(self.__wdir, "test.chrom.filtered.ge"),
            oformat="ChromFilteredGE",
        )
        GenomicElementExport.export_chrom_filtered_ge(args)

        output_bt = BedTable3()
        output_bt.load_from_file(args.opath)
        self.assertEqual(len(output_bt), 1)

    def test_export_trebed(self):
        args = argparse.Namespace(
            region_file_path=self.__bed3_path,
            region_file_type="bed3",
            pl_sig_track=self.__pl_track_path,
            mn_sig_track=self.__mn_track_path,
            opath=os.path.join(self.__wdir, "test.trebed"),
            oformat="TREbed",
        )
        GenomicElementExport.export_trebed(args)

        # Load the output TREbed file
        trebed_bt = GenomicElements.BedTableTREBed(enable_sort=False)
        trebed_bt.load_from_file(args.opath)

        self.assertEqual(len(trebed_bt), 3)

        # Get the regions to verify TSS positions
        # Region 0: chr14:75278325-75279326, peak at relative pos 100 -> fwdTSS = 75278325 + 100 = 75278425
        # Region 0: peak at relative pos 300 -> revTSS = 75278325 + 300 = 75278625
        regions = list(trebed_bt.iter_regions())
        
        # Verify first region (chr14:75278325-75279326)
        self.assertEqual(regions[0]["chrom"], "chr14")
        self.assertEqual(regions[0]["start"], 75278325)
        self.assertEqual(regions[0]["end"], 75279326)
        self.assertEqual(regions[0]["name"], "chr14:75278325-75279326")
        # fwdTSS: start (75278325) + peak position (100) = 75278425
        self.assertEqual(regions[0]["fwdTSS"], 75278325 + 100)
        # revTSS: start (75278325) + peak position (300) = 75278625
        self.assertEqual(regions[0]["revTSS"], 75278325 + 300)

        # Verify second region (chr17:45894026-45895027)
        # Region 1: peak at relative pos 200 -> fwdTSS = 45894026 + 200 = 45894226
        # Region 1: peak at relative pos 150 -> revTSS = 45894026 + 150 = 45894176
        self.assertEqual(regions[1]["chrom"], "chr17")
        self.assertEqual(regions[1]["start"], 45894026)
        self.assertEqual(regions[1]["end"], 45895027)
        self.assertEqual(regions[1]["fwdTSS"], 45894026 + 200)
        self.assertEqual(regions[1]["revTSS"], 45894026 + 150)

        # Verify third region (chr6:170553801-170554802)
        # Region 2: peak at relative pos 50 -> fwdTSS = 170553801 + 50 = 170553851
        # Region 2: peak at relative pos 250 -> revTSS = 170553801 + 250 = 170554051
        self.assertEqual(regions[2]["chrom"], "chr6")
        self.assertEqual(regions[2]["start"], 170553801)
        self.assertEqual(regions[2]["end"], 170554802)
        self.assertEqual(regions[2]["fwdTSS"], 170553801 + 50)
        self.assertEqual(regions[2]["revTSS"], 170553801 + 250)

    def test_export_trebed_negative_track(self):
        # Create negative tracks
        neg_mn_track_path = os.path.join(self.__wdir, "neg_mn_track.npz")
        
        # All regions 1001bp
        mn_track = np.zeros((3, 1001))
        mn_track[0, 300] = -15.0  # Peak magnitude at 300
        mn_track[1, 150] = -20.0  # Peak magnitude at 150
        mn_track[2, 250] = -10.0  # Peak magnitude at 250
        
        np.savez(neg_mn_track_path, mn_track)
        
        args = argparse.Namespace(
            region_file_path=self.__bed3_path,
            region_file_type="bed3",
            pl_sig_track=self.__pl_track_path,
            mn_sig_track=neg_mn_track_path,
            opath=os.path.join(self.__wdir, "test_neg.trebed"),
            oformat="TREbed",
        )
        GenomicElementExport.export_trebed(args)

        # Verify TSS positions
        trebed_bt = GenomicElements.BedTableTREBed(enable_sort=False)
        trebed_bt.load_from_file(args.opath)
        regions = list(trebed_bt.iter_regions())
        
        # Region 0: chr14:75278325-75279326
        # mn peak at 300 -> revTSS = 75278325 + 300 = 75278625
        self.assertEqual(regions[0]["revTSS"], 75278325 + 300)

    def test_export_merged_ge(self):
        left_region_path = os.path.join(self.__wdir, "merge_left.bed3")
        right_region_path = os.path.join(self.__wdir, "merge_right.bed3")

        left_bt = BedTable3(enable_sort=False)
        left_bt.load_from_dataframe(pd.DataFrame({
            "chrom": ["chr2", "chr2"],
            "start": [100, 200],
            "end": [110, 210],
        }))
        left_bt.write(left_region_path)

        right_bt = BedTable3(enable_sort=False)
        right_bt.load_from_dataframe(pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [50, 300],
            "end": [60, 310],
        }))
        right_bt.write(right_region_path)

        left_stat_path = os.path.join(self.__wdir, "left.stat.npy")
        right_stat_path = os.path.join(self.__wdir, "right.stat.npy")
        left_track_path = os.path.join(self.__wdir, "left.track.npy")
        right_track_path = os.path.join(self.__wdir, "right.track.npy")
        np.save(left_stat_path, np.array([100, 200]))
        np.save(right_stat_path, np.array([10, 20]))
        np.save(left_track_path, np.array([[1] * 10, [2] * 10]))
        np.save(right_track_path, np.array([[3] * 10, [4] * 10]))

        oheader = os.path.join(self.__wdir, "merged")
        args = argparse.Namespace(
            left_region_file_path=left_region_path,
            right_region_file_path=right_region_path,
            region_file_type="bed3",
            anno_name=["stat", "track"],
            left_anno_path=[left_stat_path, left_track_path],
            right_anno_path=[right_stat_path, right_track_path],
            oheader=oheader,
            oformat="MergedGE",
        )
        GenomicElementExport.export_merged_ge(args)

        merged_region_path = oheader + ".bed3"
        merged_bt = BedTable3()
        merged_bt.load_from_file(merged_region_path)
        merged_df = merged_bt.to_dataframe()
        self.assertEqual(merged_df.iloc[0]["chrom"], "chr1")
        self.assertEqual(merged_df.iloc[1]["chrom"], "chr1")
        self.assertEqual(merged_df.iloc[2]["chrom"], "chr2")
        self.assertEqual(merged_df.iloc[3]["chrom"], "chr2")

        stat_arr = np.load(oheader + ".stat.npy")
        track_arr = np.load(oheader + ".track.npy")
        np.testing.assert_array_equal(stat_arr.reshape(-1,), np.array([10, 20, 100, 200]))
        np.testing.assert_array_equal(track_arr[:, 0], np.array([3, 4, 1, 2]))
