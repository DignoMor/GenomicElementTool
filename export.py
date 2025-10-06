
from RGTools.GenomicElements import GenomicElements


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

class GenomicElementExport:
    @staticmethod
    def set_parser(parser):
        subparsers = parser.add_subparsers(dest="oformat")

        parser_exogeneous_sequences = subparsers.add_parser("ExogeneousSequences", 
                                                          help="Export exogeneous sequences.",
                                                          )
        GenomicElementExport.set_parser_exogeneous_sequences(parser_exogeneous_sequences)

        parser_count_table = subparsers.add_parser("CountTable", 
                                                  help="Export count table.",
                                                  )
        GenomicElementExport.set_parser_count_table(parser_count_table)

        parser_heatmap = subparsers.add_parser("Heatmap", 
                                                  help="Export heatmap.",
                                                  )
        GenomicElementExport.set_parser_heatmap(parser_heatmap)

        parser_chrom_filtered_ge = subparsers.add_parser("ChromFilteredGE",
                                                         help="Export GenomicElements with only chromosome given.", 
                                                         )
        GenomicElementExport.set_parser_chrom_filtered_ge(parser_chrom_filtered_ge)

    @staticmethod
    def set_parser_exogeneous_sequences(parser):
        GenomicElements.set_parser_genome(parser)
        GenomicElements.set_parser_genomic_element_region(parser)
        parser.add_argument("--oheader", 
                            help="Header of the output file.",
                            required=True,
                            )
        return parser

    @staticmethod
    def set_parser_count_table(parser):
        GenomicElements.set_parser_genomic_element_region(parser)
        parser.add_argument("--opath", 
                            help="Output path of the count table.",
                            required=True,
                            )

        parser.add_argument("--sample_name", 
                            help="Sample name.",
                            action="append",
                            required=True,
                            )
        
        parser.add_argument("--stat_npy", 
                            help="Path to the stat npy file.",
                            action="append",
                            required=True,
                            )
        
        parser.add_argument("--region_id_type", 
                            help="Type of the region id (default, gene_symbol).",
                            required=True,
                            default="default",
                            type=str,
                            choices=["default", "gene_symbol"],
                            )
        
        return parser
    
    @staticmethod
    def set_parser_heatmap(parser):
        GenomicElements.set_parser_genomic_element_region(parser)
        parser.add_argument("--track_npy", 
                            help="Path to the track npy file.",
                            required=True,
                            )

        parser.add_argument("--opath", 
                            help="Output path of the heatmap.",
                            required=True,
                            )
        return parser

    @staticmethod
    def set_parser_chrom_filtered_ge(parser):
        GenomicElements.set_parser_genomic_element_region(parser)
        parser.add_argument("--chrom_size", 
                            help="Chromosome size file.",
                            required=True,
                            )
        parser.add_argument("--opath", 
                            help="Output path of the filtered GenomicElements.",
                            required=True,
                            )
        return parser

    @staticmethod
    def get_oformat_options():
        return ["ExogeneousSequences", "CountTable", "Heatmap", "ChromFilteredGE"]

    @staticmethod
    def export_exogeneous_sequences(args):
        ge = GenomicElements(args.region_file_path, 
                             args.region_file_type, 
                             args.genome_path, 
                             )
        ge.export_exogeneous_sequences(args.oheader + ".fa")

    @staticmethod
    def region2region_id(region, region_id_type):
        '''
        Return the region id.

        Keyword arguments:
        - region: Region object.
        - region_id_type: Type of the region id (default, gene_symbol).

        Returns:
        - region_id: Region id.
        '''
        if region_id_type == "default":
            return region["chrom"] + ":" + str(region["start"]) + "-" + str(region["end"])
        elif region_id_type == "gene_symbol":
            if "gene_symbol" not in region.get_fields():
                raise ValueError(f"gene_symbol not supported for this region file type.")
            return region["gene_symbol"]
        else:
            raise ValueError(f"Invalid region id type: {region_id_type}")

    @staticmethod
    def export_count_table(args):
        ge = GenomicElements(args.region_file_path, 
                             args.region_file_type, 
                             None, 
                             )

        for sample_name, stat_npy in zip(args.sample_name, args.stat_npy):
            ge.load_region_anno_from_npy(sample_name, stat_npy)
        
        region_names = []
        for region in ge.get_region_bed_table().iter_regions():
            region_names.append(GenomicElementExport.region2region_id(region, args.region_id_type))

        output_df = pd.DataFrame(columns=args.sample_name, 
                                 index=region_names,
                                 )

        for sample_name in args.sample_name:
            output_df[sample_name] = ge.get_anno_arr(sample_name)

        output_df.to_csv(args.opath, 
                         )

    @staticmethod
    def export_heatmap(args):
        ge = GenomicElements(args.region_file_path, 
                             args.region_file_type, 
                             None, 
                             )
        ge.load_region_anno_from_npy("track", args.track_npy)
        track_arr = ge.get_anno_arr("track")

        width = track_arr.shape[1]

        fig, ax = plt.subplots(2, 1, 
                               figsize=(4, 8),
                               height_ratios=[1, 0.2],
                               )

        # Figure out vmin and vmax
        vmin=0
        top_1p_signal_per_track = np.percentile(track_arr, 99, axis=1)
        vmax = np.percentile(top_1p_signal_per_track, 80)

        # Figure out sorting
        sort_idx = np.argsort(track_arr.max(axis=1))

        imshow_pos = ax[0].imshow(track_arr[sort_idx, :], 
                                  cmap="Reds",
                                  aspect="auto",
                                  vmin=vmin,
                                  vmax=vmax,
                                  )

        ax[0].set_yticks([])
        ax[0].set_xticks([])
        ax[1].plot(np.arange(width), 
                   track_arr.mean(axis=0),
                   color="black",
                   )
        ax[1].set_ylabel("Mean signal")
        ax[1].set_xlabel("Position")
        ax[1].set_xticks([0, width/2, width])
        ax[1].set_xticklabels([-width//2, 0, width//2])

        fig.tight_layout()
        fig.savefig(args.opath)

    @staticmethod
    def export_chrom_filtered_ge(args):
        ge = GenomicElements(args.region_file_path, 
                             args.region_file_type, 
                             None, 
                             )
        chrom_size_df = pd.read_csv(args.chrom_size, 
                                    sep="\t", 
                                    header=None, 
                                    names=["chrom", "size"], 
                                    )
        chrom_size_dict = dict(zip(chrom_size_df["chrom"], chrom_size_df["size"]))

        filter_logical = [c in chrom_size_dict.keys() for c in ge.get_region_bed_table().get_chrom_names()]
        output_bt = ge.get_region_bed_table().apply_logical_filter(filter_logical)
        output_bt.write(args.opath)

    @staticmethod
    def main(args):
        if args.oformat == "ExogeneousSequences":
            GenomicElementExport.export_exogeneous_sequences(args)
        elif args.oformat == "CountTable":
            GenomicElementExport.export_count_table(args)
        elif args.oformat == "Heatmap":
            GenomicElementExport.export_heatmap(args)
        elif args.oformat == "ChromFilteredGE":
            GenomicElementExport.export_chrom_filtered_ge(args)
        else:
            raise ValueError(f"Invalid output format: {args.oformat}")
