
from RGTools.GenomicElements import GenomicElements

import numpy as np
import pandas as pd

class GenomicElementImport:
    @staticmethod
    def set_parser(parser):
        subparsers = parser.add_subparsers(dest="informat", 
                                           help="Input format. [{}]".format(", ".join(GenomicElementImport.get_informat_options())),
                                           required=True,
                                           )

        parser_list = subparsers.add_parser("list", 
                                            help="Import from a list file and output a GenomicElements npy/npz file.",
                                            )

        GenomicElementImport.set_parser_list(parser_list)

    @staticmethod
    def set_parser_list(parser):
        GenomicElements.set_parser_genomic_element_region(parser)

        parser.add_argument("--inpath", "-I", 
                            help="Input path of the list file.",
                            required=True,
                            )
        
        parser.add_argument("--opath", 
                            help="Output path of the region file.",
                            required=True,
                            )
        return parser

    @staticmethod
    def get_informat_options():
        return ["list"]

    @staticmethod
    def import_list(args):
        '''
        Import regions from a list file.
        '''
        input_ge = GenomicElements(region_file_path=args.region_file_path,
                                   region_file_type=args.region_file_type,
                                   fasta_path=None,
                                   )
        with open(args.inpath, 'r') as f:
            region_list = [line.strip() for line in f.readlines() if line.strip()]

        if not input_ge.get_num_regions() == len(region_list):
            raise ValueError(f"Number of regions in the input file {len(region_list)} does not match the number of regions in the region file {input_ge.get_num_regions()}")

        # Convert list of strings to numpy array
        region_arr = np.array(region_list, dtype=object)
        
        # Load the region list as an annotation
        input_ge.load_region_anno_from_arr("region_list", region_arr)

        if args.opath.endswith(".npz"): 
            input_ge.save_anno_npz("region_list", args.opath)
        elif args.opath.endswith(".npy"):
            input_ge.save_anno_npy("region_list", args.opath)
        else:
            raise ValueError(f"Invalid output file type: {args.opath}")

    @staticmethod
    def main(args):
        if args.informat == "list":
            GenomicElementImport.import_list(args)
        else:
            raise ValueError(f"Invalid input format: {args.informat}")
