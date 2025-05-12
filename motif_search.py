
import numpy as np

from RGTools.GenomicElements import GenomicElements
from RGTools.MemeMotif import MemeMotif

from RGTools.utils import str2bool

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

    @staticmethod
    def calculate_pwm_score(seq, pwm, alphabet="ACGT", bg_freq=None):
        '''
        Calculate the score of a sequence based on a given position weight matrix (PWM).
        
        Keyword arguments:
        - seq: sequence to score, must be the same length as the pwm.
        - pwm: position weight matrix, a 2D numpy array where each row 
              corresponds to a position in the sequence and each column 
              corresponds to a character in the alphabet. 
              shape = (num_positions, num_alphabet_chars)
        - bg_freq: background frequency of the alphabet, a 1D numpy array

        Returns:
        - score: the score of the sequence based on the PWM.
        '''
        if not len(seq) == pwm.shape[0]:
            raise ValueError("Length of sequence must be the same as the length of the PWM.")
        
        if not bg_freq: 
            bg_freq = np.ones(len(alphabet)) / len(alphabet)

        alphabet2idx = {char: idx for idx, char in enumerate(alphabet)}
        score = 0

        for i, char in enumerate(seq):
            score += np.log10(pwm[i, alphabet2idx[char]]) - \
                np.log10(bg_freq[alphabet2idx[char]])

        return score

    @staticmethod
    def search_one_motif(seq, motif_alphabet, motif_pwm, bg_freq=None):
        '''
        Search for a single motif in a sequence.
        Returns array of weighted score based on pwm.

        Keyword arguments:
        - seq: sequence to search
        - motif_pwm: PWM of motif to search for

        Returns:
        - output_arr: array of scores for each position in the sequence.
                      If motif length is odd, then (motif_len - 1) / 2 
                      zeros will be padded to the result on both ends. 
                      If the length is even, then motif_len / 2 - 1 zeros will be
                      padded to the left and motif_len / 2 zeros will be padded 
                      to the right.
        '''
        seq = seq.upper()

        output_arr = np.zeros(len(seq), 
                              dtype=np.float64, 
                              )

        motif_len = motif_pwm.shape[0]
        seq_len = len(seq)

        for i in range(seq_len - motif_len + 1):
            target_seq = seq[i:i + motif_len]
            score = MotifSearch.calculate_pwm_score(target_seq, motif_pwm, 
                                                    motif_alphabet, bg_freq)

            output_arr[i + int(motif_len / 2)] = score
        
        output_arr[:int(motif_len / 2)] = output_arr.min()
        output_arr[-int(motif_len / 2)+1:] = output_arr.min()

        return output_arr

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
                bg_freq = motif_dataset.get_motif_bg_freq(motif)
                bg_freq = np.array(bg_freq, dtype=np.float64)
            else:
                full_str = "".join([s.upper() for s in seq_list])
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
                motif_score_track = MotifSearch.search_one_motif(seq, motif_alphabet, motif_pwm)
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
                            required=True,
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
        filter = filter_base_score > args.min_score

        output_ge = genomic_elements.apply_logical_filter(filter, 
                                                          args.output_header + ".bed",
                                                          )
                         
        output_ge.save_anno_npy("motif", 
                                args.output_header + ".motif.npy",
                                )
