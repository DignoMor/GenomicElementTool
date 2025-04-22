#!/usr/bin/env python

import argparse

from count_bw import CountBw
from pad_region import PadRegion
from bed2tss_bed import Bed2TssBed
from one_hot import OneHot

class GenomicElementTool:
    @staticmethod
    def set_parser(parser):
        subparsers = parser.add_subparsers(dest="subcommand")

        parser_count_bw = subparsers.add_parser("count_bw",
                                                help="Count signal in bigwig files.",
                                                )

        CountBw.set_parser(parser_count_bw)

        parser_pad_region = subparsers.add_parser("pad_region",
                                                  help="Pad regions. This program differs from "
                                                       "padding_bed.py in that it conserve the " 
                                                       "order of elements in Genomic Elements files.",
                                                  )
        
        PadRegion.set_parser(parser_pad_region)

        parser_bed2tssbed = subparsers.add_parser("bed2tssbed",
                                                  help="Convert bed file to TSS bed file.",
                                                  )
        
        Bed2TssBed.set_parser(parser_bed2tssbed)

        parser_onehot = subparsers.add_parser("onehot",
                                              help="One-hot encode the sequence. Only support elements of the same size.",
                                              )
        
        OneHot.set_parser(parser_onehot)

    @staticmethod
    def main(args):
        if args.subcommand == "count_bw":
            CountBw.main(args)
        elif args.subcommand == "pad_region":
            PadRegion.main(args)
        elif args.subcommand == "bed2tssbed":
            Bed2TssBed.main(args)
        elif args.subcommand == "onehot":
            OneHot.main(args)
        else:
            raise ValueError("Unknown subcommand: {}".format(args.subcommand))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Genomic element tool.")
    GenomicElementTool.set_parser(parser)
    args = parser.parse_args()
    GenomicElementTool.main(args)
