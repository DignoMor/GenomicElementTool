
from RGTools.GenomicElements import GenomicElements
from RGTools.utils import str2bool

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

        parser_trebed = subparsers.add_parser("TREbed",
                                             help="Export TREbed file with fwdTSS and revTSS annotations.",
                                             )
        GenomicElementExport.set_parser_trebed(parser_trebed)

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
                            action="append",
                            required=True,
                            )
        
        parser.add_argument("--title", 
                            help="Title of the heatmap.",
                            type=str,
                            action="append",
                            required=True,
                            )

        parser.add_argument("--negative", 
                            help="Whether the track is negative.",
                            type=str2bool,
                            action="append",
                            required=True,
                            )

        parser.add_argument("--per_track_max_percentile", 
                            help="Percentile used to determine the maximum value per track.",
                            type=int,
                            default=99,
                            )

        parser.add_argument("--vmax_percentile", 
                            help="Percentile used to determine the vmax.",
                            type=int,
                            default=50,
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
    def set_parser_trebed(parser):
        GenomicElements.set_parser_genomic_element_region(parser)
        parser.add_argument("--pl_sig_track", 
                            help="Path to plus strand GROcap/PROcap signal track npy/npz file.",
                            required=True,
                            )
        parser.add_argument("--mn_sig_track", 
                            help="Path to minus strand GROcap/PROcap signal track npy/npz file.",
                            required=True,
                            )
        parser.add_argument("--opath", 
                            help="Output path for the TREbed file.",
                            required=True,
                            )
        return parser

    @staticmethod
    def get_oformat_options():
        return ["ExogeneousSequences", "CountTable", "Heatmap", "ChromFilteredGE", "TREbed"]

    @staticmethod
    def export_exogeneous_sequences(args):
        ge = GenomicElements(args.region_file_path, 
                             args.region_file_type, 
                             args.fasta_path, 
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
        if len(args.sample_name) != len(args.stat_npy):
            raise ValueError(f"Number of sample names ({len(args.sample_name)}) must match number of stat_npy files ({len(args.stat_npy)})")
        
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
                         index=True,
                         )

    @staticmethod
    def get_heatmap_vmin_vmax(track_arr, per_track_max_percentile, 
                              vmax_percentile):
        '''
        Get the vmin and vmax for the heatmap.
        This function first determines the maximum value per track, 
        then determines the vmax based on percentile of the determined 
        maximum values.

        Keyword arguments:
        - track_arr: Track array. Must be positive.
        - per_track_max_percentile: Percentile used to determine the maximum value per track.
        - vmax_percentile: Percentile used to determine the vmax.

        Returns:
        - vmin: Vmin.
        - vmax: Vmax.
        '''
        vmin=0
        top_1p_signal_per_track = np.percentile(track_arr, per_track_max_percentile, axis=1)
        vmax = np.percentile(top_1p_signal_per_track, vmax_percentile)

        if vmax==0: 
            vmax = 1

        return vmin, vmax

    @staticmethod
    def plot_heatmap_image(ax, track_arr, sort_idx, plot_cmap, 
                           title, per_track_max_percentile, 
                           vmax_percentile):
        '''
        Plot a heatmap with imshow.

        Keyword arguments:
        - ax: Axes object.
        - track_arr: Track array. Must be positive.
        - sort_idx: Sort index.
        - plot_cmap: Plot cmap.
        - title: Title.
        - per_track_max_percentile: Percentile used to determine the maximum value per track.
        - vmax_percentile: Percentile used to determine the vmax.

        Returns:
        - imshow_pos: Imshow position.
        '''
        vmin, vmax = GenomicElementExport.get_heatmap_vmin_vmax(track_arr, 
                                                                per_track_max_percentile, 
                                                                vmax_percentile, 
                                                                )

        imshow_pos = ax.imshow(track_arr[sort_idx, :], 
                               cmap=plot_cmap,
                               aspect="auto",
                               vmin=vmin,
                               vmax=vmax,
                               )
        
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_title(title)

        return imshow_pos

    @staticmethod
    def plot_heatmap_mean(ax, track_arr, negative_sig):
        '''
        Plot the mean signal.
        '''
        width = track_arr.shape[1]
        ax.plot(np.arange(width), 
                - track_arr.mean(axis=0) if negative_sig else track_arr.mean(axis=0),
                color="black",
                )
        ax.set_ylabel("Mean signal")
        ax.set_xlabel("Position")
        ax.set_xticks([0, width/2, width])
        ax.set_xticklabels([-width//2, 0, width//2])

    @staticmethod
    def export_heatmap(args):
        if len(args.title) != len(args.track_npy):
            raise ValueError(f"Number of titles ({len(args.title)}) must match number of track_npy files ({len(args.track_npy)})")
        if len(args.title) != len(args.negative):
            raise ValueError(f"Number of titles ({len(args.title)}) must match number of negative flags ({len(args.negative)})")
        
        ge = GenomicElements(args.region_file_path, 
                             args.region_file_type, 
                             None, 
                             )
        for track_title, track_npy in zip(args.title, args.track_npy):
            ge.load_region_anno_from_npy(track_title, track_npy)

        track_arr_list = [np.abs(ge.get_anno_arr(track_title)) for track_title in args.title]

        fig, ax = plt.subplots(2, len(args.title), 
                               figsize=(4 * len(args.title), 8),
                               height_ratios=[1, 0.2],
                               squeeze=False,
                               )

        # Figure out sorting
        sort_idx = np.argsort(np.concatenate(track_arr_list, axis=1).max(axis=1))

        for ind, track_title in enumerate(args.title):
            track_arr = np.abs(ge.get_anno_arr(track_title))

            if args.negative[ind]:
                plot_cmap = "Blues"
            else:
                plot_cmap = "Reds"

            imshow_pos = GenomicElementExport.plot_heatmap_image(ax[0, ind], 
                                                                 track_arr, 
                                                                 sort_idx, 
                                                                 plot_cmap, 
                                                                 track_title,
                                                                 args.per_track_max_percentile, 
                                                                 args.vmax_percentile,
                                                                 )

            GenomicElementExport.plot_heatmap_mean(ax[1, ind], 
                                                   track_arr, 
                                                   args.negative[ind],
                                                   )

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
    def export_trebed(args):
        '''
        Export TREbed file with forward and reverse TSS annotations.
        
        Examines GROcap/PROcap signal tracks to find TSS positions:
        - fwdTSS: Absolute genomic position of maximum signal in plus strand track
        - revTSS: Absolute genomic position of maximum signal in minus strand track
        '''
        ge = GenomicElements(args.region_file_path, 
                             args.region_file_type, 
                             None, 
                             )
        
        ge.load_region_anno_from_npy("pl_track", args.pl_sig_track)
        ge.load_region_anno_from_npy("mn_track", args.mn_sig_track)
        
        pl_track_arr = ge.get_anno_arr("pl_track")
        mn_track_arr = ge.get_anno_arr("mn_track")
        
        if pl_track_arr.shape != mn_track_arr.shape:
            raise ValueError(f"Plus and minus strand tracks must have the same shape. "
                           f"Plus: {pl_track_arr.shape}, Minus: {mn_track_arr.shape}")
        
        # Create TREbed output
        trebed_bt = GenomicElements.BedTableTREBed(enable_sort=False)
        
        output_dict_list = []
        for i, region in enumerate(ge.get_region_bed_table().iter_regions()):
            pl_track = pl_track_arr[i]
            mn_track = mn_track_arr[i]
            
            # Find TSS positions (position with maximum signal, relative to region start)
            # fwdTSS: forward TSS from plus strand track
            # revTSS: reverse TSS from minus strand track
            fwdTSS_rel_pos = np.argmax(pl_track) if len(pl_track) > 0 else 0
            revTSS_rel_pos = np.argmax(mn_track) if len(mn_track) > 0 else 0
            
            # Convert to absolute genomic coordinates
            fwdTSS = region["start"] + fwdTSS_rel_pos
            revTSS = region["start"] + revTSS_rel_pos
            
            # Create region name (default: chrom:start-end)
            region_name = f"{region['chrom']}:{region['start']}-{region['end']}"
            
            output_dict = {
                "chrom": region["chrom"],
                "start": region["start"],
                "end": region["end"],
                "name": region_name,
                "fwdTSS": fwdTSS,
                "revTSS": revTSS,
            }
            output_dict_list.append(output_dict)
        
        trebed_bt.load_from_dataframe(pd.DataFrame(output_dict_list))
        trebed_bt.write(args.opath)

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
        elif args.oformat == "TREbed":
            GenomicElementExport.export_trebed(args)
        else:
            raise ValueError(f"Invalid output format: {args.oformat}")
