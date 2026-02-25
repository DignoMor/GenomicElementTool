
import unittest
import argparse
import requests
import shutil
import os

import numpy as np

from motif_search import MotifSearch
from filter_motif_score import FilterMotifScore
from RGTools.GenomicElements import GenomicElements

class FilterMotifScoreTest(unittest.TestCase):
    def setUp(self):
        self._test_path = "filter_motif_score_test_dir"

        if not os.path.exists(self._test_path):
            os.makedirs(self._test_path)

        self._hg38_fasta_path = os.path.join("RGTools", "large_files", "hg38.fa")
        self._bed6_path = os.path.join("example_data", "three_genes.bed6")
        self._meme_motif_path = os.path.join(self._test_path, "test_motifs.meme")
        url = "https://meme-suite.org/meme/doc/examples/sample-dna-motif.meme"
        response = requests.get(url)

        # Save the file locally
        with open(self._meme_motif_path, "wb") as f:
            f.write(response.content)

        # Generate motif scores first for filtering tests
        args = argparse.Namespace()
        args.subcommand = "motif_search"
        args.fasta_path = self._hg38_fasta_path
        args.region_file_path = self._bed6_path
        args.region_file_type = "bed6"
        args.motif_file = self._meme_motif_path
        args.output_header = os.path.join(self._test_path, "three_genes.motif_search")
        args.estimate_background_freq = True
        args.strand = "-"
        MotifSearch.main(args)

    def tearDown(self):
        if os.path.exists(self._test_path):
            shutil.rmtree(self._test_path)
        super().tearDown()

    def get_filter_motif_score_simple_args(self):
        args = argparse.Namespace()

        args.subcommand = "filter_motif_score"
        args.region_file_path = self._bed6_path
        args.region_file_type = "bed6"
        args.motif_search_npy = os.path.join(self._test_path, "three_genes.motif_search.crp.npy")
        args.output_header = os.path.join(self._test_path, "three_genes.crp.filtered")
        args.filter_base = 649
        args.min_score = 0.0
        args.max_score = np.inf

        return args
        
    def test_filter_motif_score(self):
        args = self.get_filter_motif_score_simple_args()
        FilterMotifScore.main(args)

        filtered_ge = GenomicElements(region_file_path=args.output_header + ".bed",
                                      region_file_type=args.region_file_type,
                                      fasta_path=None, 
                                      )
        
        filtered_ge.load_region_anno_from_npy("motif", 
                                              args.output_header + ".motif.npy",
                                              anno_type="track",
                                              )
        
        self.assertEqual(len(filtered_ge.get_track_list("motif")), 1)

if __name__ == "__main__":
    unittest.main()

