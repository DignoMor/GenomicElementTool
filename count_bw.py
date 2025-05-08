

import numpy as np

from RGTools.GenomicElements import GenomicElements
from RGTools.BwTrack import SingleBwTrack, PairedBwTrack
from RGTools.utils import str2bool

class CountSingleBw:
    @staticmethod
    def set_parser(parser):
        GenomicElements.set_parser_genomic_element_region(parser)
        parser.add_argument("--bw_path",
                            help="Bigwig file path.",
                            required=True,
                            type=str,
                            )

        parser.add_argument("--quantification_type",
                            help="Type of quantification.",
                            type=str,
                            default="raw_count",
                            choices=SingleBwTrack.get_supported_quantification_type(),
                            )

        parser.add_argument("--opath",
                            help="Output path for counting.",
                            required=True,
                            type=str,
                            )

    @staticmethod
    def main(args):
        genomic_elements = GenomicElements(region_path=args.region_file_path,
                                           region_file_type=args.region_file_type,
                                           fasta_path=None, 
                                           )
        region_bt = genomic_elements.get_region_bed_table()

        bw_track = SingleBwTrack(bw_path=args.bw_path)

        output_list = []
        for region in region_bt.iter_regions():

            output_list.append(bw_track.count_single_region(region["chrom"],
                                                            region["start"],
                                                            region["end"],
                                                            output_type=args.quantification_type,
                                                            min_len_after_padding=1, 
                                                            ),
                               )

        output_arr = np.array(output_list)

        np.save(args.opath, output_arr)

class CountPairedBw:
    @staticmethod
    def set_parser(parser):
        GenomicElements.set_parser_genomic_element_region(parser)
        parser.add_argument("--bw_pl",
                            help="Plus strand bigwig file.",
                            required=True,
                            type=str,
                            )

        parser.add_argument("--bw_mn",
                            help="Minus strand bigwig file.",
                            type=str,
                            default=None,
                            )

        parser.add_argument("--override_strand",
                            help="Override the strand information in the input file (None if use the input strand info).",
                            type=str,
                            default=None,
                            )

        parser.add_argument("--quantification_type",
                            help="Type of quantification.",
                            type=str,
                            default="raw_count",
                            choices=PairedBwTrack.get_supported_quantification_type(),
                            )

        parser.add_argument("--opath",
                            help="Output path for counting.",
                            required=True,
                            type=str,
                            )

    @staticmethod
    def main(args):
        genomic_elements = GenomicElements(region_path=args.region_file_path,
                                           region_file_type=args.region_file_type,
                                           fasta_path=None, 
                                           )
        region_bt = genomic_elements.get_region_bed_table()

        bw_track = PairedBwTrack(bw_pl_path=args.bw_pl,
                                 bw_mn_path=args.bw_mn,
                                 )

        output_list = []
        for region in region_bt.iter_regions():
            if args.override_strand:
                strand = args.override_strand
            else:
                strand = "." if args.region_file_type == "bed3" else region["strand"]

            output_list.append(bw_track.count_single_region(region["chrom"],
                                                            region["start"],
                                                            region["end"],
                                                            strand, 
                                                            output_type=args.quantification_type,
                                                            min_len_after_padding=1, 
                                                            ),
                               )

        output_arr = np.array(output_list)

        np.save(args.opath, output_arr)
    
