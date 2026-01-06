
import unittest
import argparse
import requests
import shutil
import os

import numpy as np

from motif_search import MotifSearch

from RGTools.GenomicElements import GenomicElements

class MotifSearchTest(unittest.TestCase):
    def setUp(self):
        self._test_path = "motif_search_test_dir"

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

    def tearDown(self):
        if os.path.exists(self._test_path):
            shutil.rmtree(self._test_path)
        super().tearDown()
    
    def get_motif_search_simple_args(self):
        args = argparse.Namespace()

        args.subcommand = "motif_search"
        args.fasta_path = self._hg38_fasta_path
        args.region_file_path = self._bed6_path
        args.region_file_type = "bed6"
        args.motif_file = self._meme_motif_path
        args.output_header = os.path.join(self._test_path, "three_genes.motif_search")
        args.estimate_background_freq = True
        args.strand = "+"

        return args

    def test_main(self):
        args = self.get_motif_search_simple_args()

        # Call the main function of the MotifSearch class
        MotifSearch.main(args)

        output_ge = GenomicElements(region_file_path=args.region_file_path,
                                    region_file_type=args.region_file_type,
                                    fasta_path=args.fasta_path, 
                                    )
        
        output_ge.load_region_anno_from_npy("CRP", 
                                            os.path.join(args.output_header + ".crp.npy"), 
                                            )

        output_ge.load_region_anno_from_npy("LexA", 
                                            os.path.join(args.output_header + ".lexA.npy"), 
                                            )

        crp_track = output_ge.get_anno_arr("CRP")
        lexA_track = output_ge.get_anno_arr("LexA")

        self.assertEqual(crp_track.shape[0], 3)
        self.assertEqual(lexA_track.shape[1], 1001)

        # Test that scores are reasonable (not all zeros or inf)
        self.assertFalse(np.all(crp_track == 0))
        self.assertFalse(np.all(np.isinf(crp_track)))
    
    def test_main_reverse_strand(self):
        """Test motif search with reverse strand (strand='-')"""
        args = self.get_motif_search_simple_args()
        args.strand = "-"

        MotifSearch.main(args)

        output_ge = GenomicElements(region_file_path=args.region_file_path,
                                    region_file_type=args.region_file_type,
                                    fasta_path=args.fasta_path, 
                                    )
        
        output_ge.load_region_anno_from_npy("CRP", 
                                            os.path.join(args.output_header + ".crp.npy"), 
                                            )

        crp_track = output_ge.get_anno_arr("CRP")

        self.assertEqual(crp_track.shape[0], 3)
        self.assertEqual(crp_track.shape[1], 1001)
        
        # Test that scores are reasonable
        self.assertFalse(np.all(crp_track == 0))
        self.assertFalse(np.all(np.isinf(crp_track)))
    
    def test_main_both_strands(self):
        """Test motif search with both strands (strand='both')"""
        args = self.get_motif_search_simple_args()
        args.strand = "both"

        MotifSearch.main(args)

        output_ge = GenomicElements(region_file_path=args.region_file_path,
                                    region_file_type=args.region_file_type,
                                    fasta_path=args.fasta_path, 
                                    )
        
        output_ge.load_region_anno_from_npy("CRP", 
                                            os.path.join(args.output_header + ".crp.npy"), 
                                            )

        crp_track = output_ge.get_anno_arr("CRP")

        self.assertEqual(crp_track.shape[0], 3)
        self.assertEqual(crp_track.shape[1], 1001)
        
        # Test that scores are reasonable
        self.assertFalse(np.all(crp_track == 0))
        self.assertFalse(np.all(np.isinf(crp_track)))
        
        # Test that 'both' gives scores >= forward strand scores
        args_fwd = self.get_motif_search_simple_args()
        args_fwd.strand = "+"
        MotifSearch.main(args_fwd)
        output_ge_fwd = GenomicElements(region_file_path=args_fwd.region_file_path,
                                        region_file_type=args_fwd.region_file_type,
                                        fasta_path=args_fwd.fasta_path, 
                                        )
        output_ge_fwd.load_region_anno_from_npy("CRP", 
                                                 os.path.join(args_fwd.output_header + ".crp.npy"), 
                                                 )
        crp_track_fwd = output_ge_fwd.get_anno_arr("CRP")
        
        # 'both' should give max of forward and reverse, so should be >= forward
        self.assertTrue(np.all(crp_track >= crp_track_fwd))
    
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
        args = self.get_motif_search_simple_args()
        args.strand = "-"
        MotifSearch.main(args)

        args = self.get_filter_motif_score_simple_args()
        MotifSearch.filter_motif_score_main(args)

        filtered_ge = GenomicElements(region_file_path=args.output_header + ".bed",
                                      region_file_type=args.region_file_type,
                                      fasta_path=None, 
                                      )
        
        filtered_ge.load_region_anno_from_npy("motif", 
                                              args.output_header + ".motif.npy",
                                              )
        
        self.assertEqual(filtered_ge.get_anno_arr("motif").shape[0], 1)
