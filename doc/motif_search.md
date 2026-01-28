---
title: motif_search Subcommand
description: Search for motifs in genomic sequences using MEME format Position Weight Matrices
---

# motif_search Subcommand

Search for motifs in genomic sequences using Position Weight Matrices (PWMs) from MEME format files. This subcommand scans each genomic region and computes motif match scores at every position, outputting score tracks for each motif.

## Usage

```bash
GenomicElementTool.py motif_search [OPTIONS]
```

## Required Arguments

- `--fasta_path` (str)
  - Path to the genome FASTA file
  - Required: Yes

- `--region_file_path` (str)
  - Path to the region file (BED format)
  - Required: Yes

- `--region_file_type` (str)
  - Type of the region file
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, `bedTRE`, etc.
  - See `GenomicElements` documentation for full list

- `--motif_file` (str)
  - Path to motif file in MEME format
  - Required: Yes
  - Must contain one or more motifs in standard MEME format

## Optional Arguments

- `--output_header` (str)
  - Header prefix for output files
  - Default: `"motif_search"`
  - Output files are saved as `<output_header>.<motif_name>.npy`
  - Each motif in the MEME file generates a separate output file

- `--estimate_background_freq` (bool)
  - Estimate background nucleotide frequencies from the input sequences
  - Default: `True`
  - If `True`: Calculates background frequencies from concatenated sequences
  - If `False`: Uses background frequencies from the MEME file
  - Note: N frequencies are always estimated from sequences when present

- `--strand` (str)
  - Search reverse complement sequences for motif matches. Choices: "+", "-", "both".
  - Default: `"+"`
  - If `"+"`: Search the fasta sequence with the given pwm.
  - If `"-"`: search the reverse complemented fasta sequence with the given pwm.
    (The `output[i]` is the matching score of `RC(seq[i, i+l])`)
  - If `"both"`: `output[i]` is the higher score between matching `seq[i, i+l]` 
    and `RC(seq[i, i+l])`.

## Output

For each motif in the MEME file, the program generates:

1. **NumPy array file**: `<output_header>.<motif_name>.npy`
   - Shape: `(num_regions, region_length)`
   - Contains motif match scores for each position in each region
   - Scores are log-odds (log10 ratios)
   - Higher scores indicate better motif matches
   - Positions where the motif doesn't fit (near sequence ends) are set to minimum score
   - Can be loaded later using `GenomicElement.load_region_anno_from_npy()` as a 
     Genomic Element annotation.

## How It Works

1. **Load genomic regions**: Reads regions from the specified BED file
2. **Extract sequences**: Retrieves sequences from the FASTA file for each region
3. **Load motifs**: Parses all motifs from the MEME format file
4. **For each motif**:
   - Loads the Position Weight Matrix (PWM)
   - Estimates or uses background frequencies
   - Adds pseudo-counts to the PWM (Laplace smoothing)
   - Searches each sequence for motif matches using sliding window
   - Computes log-odds scores at each position
   - Saves scores to a NumPy array file

## Score Calculation

- Scores are computed as log-odds ratios: `log10(P(motif|sequence) / P(motif|background))`
- Pseudo-counts are added to prevent zero probabilities: `(counts + 1) / (total + alphabet_size)`
- For motifs of length L, the last L-1 positions in each sequence receive minimum scores (motif doesn't fit). 
  This is to maintain the signal_track `dimension == sequence length` invariable.

## Examples

### Basic usage

```bash
GenomicElementTool.py motif_search \
    --fasta_path /path/to/genome.fa \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --motif_file motifs.meme \
    --output_header my_motif_search
```

This will create files like:
- `my_motif_search.Motif1.npy`
- `my_motif_search.Motif2.npy`
- etc. (one per motif in the MEME file)

### Search reverse complement sequences

```bash
GenomicElementTool.py motif_search \
    --fasta_path /path/to/genome.fa \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --motif_file motifs.meme \
    --strand + \
    --output_header rc_motif_search
```

### Use MEME file background frequencies

```bash
GenomicElementTool.py motif_search \
    --fasta_path /path/to/genome.fa \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --motif_file motifs.meme \
    --estimate_background_freq False \
    --output_header motif_search
```

## Loading Results

To load the motif search results:

```python
from RGTools.GenomicElements import GenomicElements

ge = GenomicElements(
    region_file_path="regions.bed6",
    region_file_type="bed6",
    fasta_path=None  # Not needed for loading annotations
)

# Load motif scores
ge.load_region_anno_from_npy("Motif1", "my_motif_search.Motif1.npy")

# Get scores as list of numpy arrays
scores = ge.get_anno_list("Motif1")
# Length: num_regions, each element: array of shape (region_length,)

# Find best match position in first region
best_pos = np.argmax(scores[0])
best_score = scores[0][best_pos]
```

## Dependencies

- numpy
- RGTools.GenomicElements
- RGTools.MemeMotif
- BioPython (for sequence handling)

## Notes

- All sequences are converted to uppercase before processing
- N nucleotides in sequences are handled by extending the alphabet
- The program processes all motifs in the MEME file sequentially
- Output filenames use motif names exactly as they appear in the MEME file
- For large datasets, consider processing motifs separately or using parallelization

