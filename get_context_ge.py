
from RGTools.GenomicElements import GenomicElements


class GetContextGe:
    @staticmethod
    def set_parser(parser):
        method_subparsers = parser.add_subparsers(dest="method")

        nearest_parser = method_subparsers.add_parser(
            "nearest",
            help="Select nearest context region by closest-end distance.",
        )

        GenomicElements.set_parser_genomic_element_region(nearest_parser)
        nearest_parser.add_argument(
            "--context_file_path",
            help="Path to the context region file.",
            required=True,
            type=str,
        )
        nearest_parser.add_argument(
            "--context_file_type",
            help="Type of the context region file. "
                 "Valid types: {}".format(
                     list(GenomicElements.get_region_file_suffix2class_dict().keys())
                 ),
            required=True,
            type=str,
            choices=GenomicElements.get_region_file_suffix2class_dict().keys(),
        )
        nearest_parser.add_argument(
            "--opath",
            help="Path to the output region file.",
            required=True,
            type=str,
        )

    @staticmethod
    def _closest_end_distance(region, context_region):
        # Minimal absolute distance between any pair of interval ends.
        return min(
            abs(region["start"] - context_region["start"]),
            abs(region["start"] - context_region["end"]),
            abs(region["end"] - context_region["start"]),
            abs(region["end"] - context_region["end"]),
        )

    @staticmethod
    def _run_nearest(args):
        input_ge = GenomicElements(
            region_file_path=args.region_file_path,
            region_file_type=args.region_file_type,
            fasta_path=None,
        )
        context_ge = GenomicElements(
            region_file_path=args.context_file_path,
            region_file_type=args.context_file_type,
            fasta_path=None,
        )

        input_regions = list(input_ge.get_region_bed_table().iter_regions())
        context_bt = context_ge.get_region_bed_table()
        context_regions = list(context_bt.iter_regions())

        context_by_chrom = {}
        for context_region in context_regions:
            context_by_chrom.setdefault(context_region["chrom"], []).append(context_region)

        output_regions = []
        for region in input_regions:
            same_chrom_contexts = context_by_chrom.get(region["chrom"], [])
            if len(same_chrom_contexts) == 0:
                raise ValueError(
                    f"No context regions found on chromosome '{region['chrom']}' "
                    f"for input region {region['chrom']}:{region['start']}-{region['end']}."
                )

            nearest_context = min(
                same_chrom_contexts,
                key=lambda context_region: GetContextGe._closest_end_distance(region, context_region),
            )
            output_regions.append(nearest_context)

        output_bt = context_bt._clone_empty()
        output_bt.load_from_bed_regions(output_regions)
        output_bt.write(args.opath)

    @staticmethod
    def main(args):
        if args.method == "nearest":
            GetContextGe._run_nearest(args)
        else:
            raise ValueError(f"Unknown get_context_ge method: {args.method}")
