
from RGTools.GenomicElements import GenomicElements

import pandas as pd

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
    def get_oformat_options():
        return ["ExogeneousSequences", "CountTable"]

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
    def main(args):
        if args.oformat == "ExogeneousSequences":
            GenomicElementExport.export_exogeneous_sequences(args)
        elif args.oformat == "CountTable":
            GenomicElementExport.export_count_table(args)
        else:
            raise ValueError(f"Invalid output format: {args.oformat}")
