
from RGTools.GenomicElements import GenomicElements

class GenomicElementExport:
    @staticmethod
    def set_parser(parser):
        GenomicElements.set_parser_genome(parser)
        GenomicElements.set_parser_genomic_element_region(parser)
        
        parser.add_argument("--oheader", 
                            help="Header of the output file.",
                            required=True,
                            )

        parser.add_argument("--oformat", 
                            help="Format of the output file.",
                            required=True,
                            )

    @staticmethod
    def get_oformat_options():
        return ["ExogeneousSequences"]

    @staticmethod
    def export_exogeneous_sequences(args):
        ge = GenomicElements(args.region_path, 
                             args.region_file_type, 
                             args.genome_path, 
                             )
        ge.export_exogeneous_sequences(args.oheader + ".fa")

    @staticmethod
    def main(args):
        if args.oformat == "ExogeneousSequences":
            GenomicElementExport.export_exogeneous_sequences(args)
        else:
            raise ValueError(f"Invalid output format: {args.oformat}")
