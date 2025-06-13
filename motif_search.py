
import numpy as np

from RGTools.GenomicElements import GenomicElements
from RGTools.MemeMotif import MemeMotif

from RGTools.utils import str2bool
from RGTools.utils import reverse_complement as RC

class MotifSearch:
    @staticmethod
    def set_parser(parser):
        GenomicElements.set_parser_genome(parser)
        GenomicElements.set_parser_genomic_element_region(parser)

        parser.add_argument("--motif_file",
                            help="Motif file in MEME format.",
                            required=True,
                            )
    
        parser.add_argument("--output_header", 
                            help="Header for the output file. Output " \
                                 "will be saved as <output_header>.<motif_name>.npy",
                            type=str,
                            default="motif_search",
                            )
        
        parser.add_argument("--estimate_background_freq",
                            help="Estimate background frequency from the sequence.",
                            type=str2bool, 
                            default=True,
                            )
        
        parser.add_argument("--reverse_complement",
                            help="Reverse complement the sequence while matching for motifs.",
                            type=str2bool,
                            default=False,
                            )

    @staticmethod
    def main(args):
        genomic_elements = GenomicElements(region_path=args.region_file_path,
                                           region_file_type=args.region_file_type,
                                           fasta_path=args.fasta_path, 
                                           )
        motif_dataset = MemeMotif(args.motif_file)

        seq_list = genomic_elements.get_all_region_seqs()

        for motif in motif_dataset.get_motif_list():
            motif_pwm = motif_dataset.get_motif_pwm(motif)
            motif_alphabet = motif_dataset.get_alphabet()

            # Estimate background frequency
            if not args.estimate_background_freq:
                bg_freq = motif_dataset.get_bg_freq()
                bg_freq = np.array(bg_freq, dtype=np.float64)

                # still estimate N frequency from the sequence
                full_str = "".join([s.upper() for s in seq_list])
                if "N" in full_str:
                    motif_alphabet += "N"
                    n_count = np.char.count(full_str, "N")
                    bg_freq = np.concatenate((bg_freq, 
                                              [n_count / (len(full_str) - n_count)], 
                                              ))
                    bg_freq /= np.sum(bg_freq)
            else:
                if not args.reverse_complement:
                    full_str = "".join([s.upper() for s in seq_list])
                else:
                    full_str = "".join([RC(s.upper()) for s in seq_list])
                if "N" in full_str:
                    motif_alphabet += "N"
                motif_pwm = np.concatenate((motif_pwm, np.zeros((motif_pwm.shape[0], 1), dtype=np.float64)), axis=1)
                bg_freq = np.array([int(np.char.count(full_str, a)) for a in motif_alphabet], 
                                   dtype=np.float64)
                bg_freq /= np.sum(bg_freq)

            # Add pseudo count to the PWM
            total_counts = motif_dataset.get_motif_num_source_sites(motif) 
            total_count_matrix = motif_pwm * total_counts
            motif_pwm = (total_count_matrix + 1) / (total_counts + motif_pwm.shape[1])
            
            output_score_anno_list = []
            for seq in seq_list:
                motif_score_track = MemeMotif.search_one_motif(seq, 
                                                               motif_alphabet, 
                                                               motif_pwm, 
                                                               reverse_complement=args.reverse_complement, 
                                                               )
                output_score_anno_list.append(motif_score_track)
            output_score_anno_arr = np.array(output_score_anno_list)

            genomic_elements.load_region_anno_from_arr(motif, 
                                                       output_score_anno_arr, 
                                                       )
            genomic_elements.save_anno_npy(motif, 
                                           args.output_header + "." + motif + ".npy", 
                                           )

    @staticmethod
    def set_filter_motif_score_args(parser):
        GenomicElements.set_parser_genomic_element_region(parser)

        parser.add_argument("--motif_search_npy",
                            help="Numpy file containing motif search scores.",
                            required=True,
                            )

        parser.add_argument("--output_header",
                            help="Output header for filtered GenomicElements bed file.",
                            required=True,
                            )
    
        parser.add_argument("--filter_base", 
                            help="Base index for filtering motif search scores.",
                            type=int,
                            required=True,
                            )
        
        parser.add_argument("--min_score",
                            help="Minimum score for filtering motif search scores.",
                            type=float,
                            default=-np.inf,
                            )

        parser.add_argument("--max_score", 
                            help="Maximum score for filtering motif search scores.",
                            type=float,
                            default=np.inf,
                            )
        
    @staticmethod
    def filter_motif_score_main(args):
        genomic_elements = GenomicElements(region_path=args.region_file_path,
                                           region_file_type=args.region_file_type,
                                           fasta_path=None, 
                                           )
        genomic_elements.load_region_anno_from_npy("motif", args.motif_search_npy)

        motif_track = genomic_elements.get_anno_arr("motif")
        filter_base_score = motif_track[:, args.filter_base]
        filter = (filter_base_score > args.min_score) & (filter_base_score < args.max_score)

        output_ge = genomic_elements.apply_logical_filter(filter, 
                                                          args.output_header + ".bed",
                                                          )
                         
        output_ge.save_anno_npy("motif", 
                                args.output_header + ".motif.npy",
                                )
