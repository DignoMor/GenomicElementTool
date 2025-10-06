
import unittest
import argparse
import shutil
import os

import pandas as pd
import numpy as np

from export import GenomicElementExport
from RGTools.BedTable import BedTable3

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
            genome_path=self.__fasta_path,
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
