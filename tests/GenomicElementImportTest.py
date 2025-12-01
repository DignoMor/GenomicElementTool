
import unittest
import argparse
import shutil
import os

import numpy as np

from GenomicElementImport import GenomicElementImport
from RGTools.GenomicElements import GenomicElements

class GenomicElementImportTest(unittest.TestCase):
    def setUp(self):
        self.__bed3_path = os.path.join("example_data", "three_genes.bed3")
        self.__bed6gene_path = os.path.join("example_data", "three_genes.bed6gene")
        self.__wdir = "GenomicElementImport_test"
        if not os.path.exists(self.__wdir):
            os.makedirs(self.__wdir)
        
        super().setUp()

    def tearDown(self):
        if os.path.exists(self.__wdir):
            shutil.rmtree(self.__wdir)

        super().tearDown()

    def test_import_list_npy(self):
        # Create a test list file with 3 lines (matching three_genes.bed3)
        list_file = os.path.join(self.__wdir, "test_list.txt")
        with open(list_file, 'w') as f:
            f.write("region1\n")
            f.write("region2\n")
            f.write("region3\n")

        args = argparse.Namespace(
            region_file_path=self.__bed3_path,
            region_file_type="bed3",
            inpath=list_file,
            opath=os.path.join(self.__wdir, "test_region_list.npy"),
            informat="list",
            dtype="str"
        )
        GenomicElementImport.import_list(args)

    def test_import_int_list(self):
        # Create a test list file with 3 lines (matching three_genes.bed3)
        list_file = os.path.join(self.__wdir, "test_list.txt")
        with open(list_file, 'w') as f:
            f.write("1\n")
            f.write("2\n")
            f.write("3\n")

        args = argparse.Namespace(
            region_file_path=self.__bed3_path,
            region_file_type="bed3",
            inpath=list_file,
            opath=os.path.join(self.__wdir, "test_region_list.npy"),
            informat="list",
            dtype="np.int32"
        )
        GenomicElementImport.import_list(args)

        # Load and verify the saved annotation
        loaded_arr = np.load(args.opath)
        self.assertEqual(loaded_arr.dtype, np.int32)
        self.assertEqual(len(loaded_arr), 3)
        self.assertEqual(loaded_arr[0], 1)
        self.assertEqual(loaded_arr[1], 2)
        self.assertEqual(loaded_arr[2], 3)

    def test_import_list_npz(self):
        # Create a test list file with 3 lines (matching three_genes.bed3)
        list_file = os.path.join(self.__wdir, "test_list.txt")
        with open(list_file, 'w') as f:
            f.write("region1\n")
            f.write("region2\n")
            f.write("region3\n")

        args = argparse.Namespace(
            region_file_path=self.__bed3_path,
            region_file_type="bed3",
            inpath=list_file,
            opath=os.path.join(self.__wdir, "test_region_list.npz"),
            informat="list",
            dtype="str"
        )
        GenomicElementImport.import_list(args)

        # Verify the output file exists
        self.assertTrue(os.path.exists(args.opath))
        
        # Load and verify the saved annotation
        loaded_data = np.load(args.opath, allow_pickle=True)
        # npz files contain arrays with keys, get the first (and only) array
        loaded_arr = loaded_data[loaded_data.files[0]]
        self.assertEqual(len(loaded_arr), 3)
        self.assertEqual(loaded_arr[0], "region1")
        self.assertEqual(loaded_arr[1], "region2")
        self.assertEqual(loaded_arr[2], "region3")

    def test_import_list_mismatch_count(self):
        # Create a test list file with wrong number of lines
        list_file = os.path.join(self.__wdir, "test_list.txt")
        with open(list_file, 'w') as f:
            f.write("region1\n")
            f.write("region2\n")
            # Only 2 lines, but bed3 has 3 regions

        args = argparse.Namespace(
            region_file_path=self.__bed3_path,
            region_file_type="bed3",
            inpath=list_file,
            opath=os.path.join(self.__wdir, "test_region_list.npy"),
            informat="list",
            dtype="str"
        )
        
        with self.assertRaises(ValueError) as context:
            GenomicElementImport.import_list(args)
        
        self.assertIn("does not match", str(context.exception))

    def test_import_list_invalid_output_type(self):
        # Create a test list file
        list_file = os.path.join(self.__wdir, "test_list.txt")
        with open(list_file, 'w') as f:
            f.write("region1\n")
            f.write("region2\n")
            f.write("region3\n")

        args = argparse.Namespace(
            region_file_path=self.__bed3_path,
            region_file_type="bed3",
            inpath=list_file,
            opath=os.path.join(self.__wdir, "test_region_list.txt"),  # Invalid extension
            dtype="str",
            informat="list",
        )
        
        with self.assertRaises(ValueError) as context:
            GenomicElementImport.import_list(args)
        
        self.assertIn("Invalid output file type", str(context.exception))


