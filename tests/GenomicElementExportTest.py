
import unittest
import argparse
import shutil
import os

from export import GenomicElementExport

class GenomicElementExportTest(unittest.TestCase):
    def setUp(self):
        self.__bed3_path = os.path.join("example_data", "three_genes.bed3")
        self.__fasta_path = os.path.join("RGTools", "large_files", "hg38.fa")
        self.__wdir = "GenomicElementExport_test"
        if not os.path.exists(self.__wdir):
            os.makedirs(self.__wdir)
        
        super().setUp()

    def tearDown(self):
        if os.path.exists(self.__wdir):
            shutil.rmtree(self.__wdir)

        super().tearDown()

    def test_export_exogeneous_sequences(self):
        args = argparse.Namespace(
            region_path=self.__bed3_path,
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

