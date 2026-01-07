---
title: export Subcommand
description: Export Genomic Element data to various output formats
---

# export Subcommand

The `export` subcommand provides functionality to export Genomic Element data to various output formats including FASTA sequences, count tables, heatmaps, and filtered region files.

## Usage

```bash
GenomicElementTool.py export <oformat> [OPTIONS]
```

The export subcommand supports five output formats:
- `ExogeneousSequences`: Export genomic sequences to FASTA format
- `CountTable`: Export statistical annotations as a count table (CSV)
- `Heatmap`: Generate heatmap visualizations from track annotations
- `ChromFilteredGE`: Filter regions by chromosome and export as BED file
- `TREbed`: Annotate regions with forward and reverse TSS from GROcap/PROcap signals

## ExogeneousSequences

Export genomic regions as sequences in FASTA format.

### Usage

```bash
GenomicElementTool.py export ExogeneousSequences [OPTIONS]
```

### Required Arguments

- `--fasta_path` (str)
  - Path to the genome FASTA file
  - Required: Yes

- `--region_file_path` (str)
  - Path to the region file (BED format)
  - Required: Yes

- `--region_file_type` (str)
  - Type of the region file
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, `TREbed`, etc.
  - See `GenomicElements` documentation for full list

- `--oheader` (str)
  - Header prefix for the output FASTA file
  - Required: Yes
  - Output file will be saved as `<oheader>.fa`

### Output

- **FASTA file**: `<oheader>.fa`
  - Contains one sequence per genomic region
  - Sequence headers use format: `>chrom:start-end`
  - Sequences are extracted from the reference genome based on region coordinates

### Example

```bash
GenomicElementTool.py export ExogeneousSequences \
    --fasta_path /path/to/genome.fa \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --oheader my_sequences
```

This creates `my_sequences.fa` with sequences for each region in the BED file.

## CountTable

Export statistical annotations (stat arrays) as a count table in CSV format. This is useful for downstream analysis with tools like DESeq2 or edgeR.

### Usage

```bash
GenomicElementTool.py export CountTable [OPTIONS]
```

### Required Arguments

- `--region_file_path` (str)
  - Path to the region file (BED format)
  - Required: Yes

- `--region_file_type` (str)
  - Type of the region file
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, etc.

- `--opath` (str)
  - Output path for the count table CSV file
  - Required: Yes

- `--sample_name` (str)
  - Sample name for the annotation
  - Required: Yes
  - Can be specified multiple times (use `--sample_name` once per sample)
  - Must match the number of `--stat_npy` arguments

- `--stat_npy` (str)
  - Path to a stat annotation NumPy file (`.npy`)
  - Required: Yes
  - Can be specified multiple times (use `--stat_npy` once per sample)
  - Must match the number of `--sample_name` arguments
  - Each file should contain a 1D array with shape `(num_regions,)`

### Optional Arguments

- `--region_id_type` (str)
  - Type of identifier to use for regions in the output table
  - Default: `"default"`
  - Choices: `"default"`, `"gene_symbol"`
  - `"default"`: Uses format `chrom:start-end`
  - `"gene_symbol"`: Uses gene symbol from region file (requires `bed6gene` or similar format with gene_symbol field)

### Output

- **CSV file**: Contains a count table with:
  - Rows: One per genomic region (identified by `region_id_type`)
  - Columns: One per sample (named by `--sample_name`)
  - Values: Statistical values from the stat annotation arrays
  - Index: Region identifiers (first column when read with `index_col=0`)

### Example

```bash
GenomicElementTool.py export CountTable \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --opath count_table.csv \
    --sample_name sample1 --sample_name sample2 \
    --stat_npy sample1_stats.npy --stat_npy sample2_stats.npy
```

With gene symbols:

```bash
GenomicElementTool.py export CountTable \
    --region_file_path genes.bed6gene \
    --region_file_type bed6gene \
    --opath count_table.csv \
    --region_id_type gene_symbol \
    --sample_name sample1 --sample_name sample2 \
    --stat_npy sample1_stats.npy --stat_npy sample2_stats.npy
```

## Heatmap

Generate heatmap visualizations from track annotations. Creates publication-quality heatmaps with mean signal plots.

### Usage

```bash
GenomicElementTool.py export Heatmap [OPTIONS]
```

### Required Arguments

- `--region_file_path` (str)
  - Path to the region file (BED format)
  - Required: Yes

- `--region_file_type` (str)
  - Type of the region file
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, etc.

- `--opath` (str)
  - Output path for the heatmap image file
  - Required: Yes
  - Recommended formats: `.png`, `.pdf`, `.svg`

- `--track_npy` (str)
  - Path to a track annotation NumPy file (`.npy`)
  - Required: Yes
  - Can be specified multiple times (use `--track_npy` once per track)
  - Must match the number of `--title` and `--negative` arguments
  - Each file should contain a 2D array with shape `(num_regions, region_length)`

- `--title` (str)
  - Title for the heatmap track
  - Required: Yes
  - Can be specified multiple times (use `--title` once per track)
  - Must match the number of `--track_npy` and `--negative` arguments

- `--negative` (bool)
  - Whether the track represents negative signal (e.g., minus strand)
  - Required: Yes
  - Can be specified multiple times (use `--negative` once per track)
  - Must match the number of `--track_npy` and `--title` arguments
  - If `True`: Uses "Blues" colormap and inverts mean signal plot
  - If `False`: Uses "Reds" colormap

### Optional Arguments

- `--per_track_max_percentile` (int)
  - Percentile used to determine the maximum value per track for normalization
  - Default: `99`
  - Range: 0-100
  - Higher values allow more extreme outliers to influence the color scale

- `--vmax_percentile` (int)
  - Percentile used to determine the overall vmax across all tracks
  - Default: `50`
  - Range: 0-100
  - Controls the color scale normalization across multiple tracks

### Output

- **Image file**: Contains a multi-panel heatmap with:
  - **Top row**: Heatmap visualizations (one per track)
    - Regions are sorted by maximum signal across all tracks
    - Color scale is normalized using percentile-based method
    - Colormap: "Reds" for positive tracks, "Blues" for negative tracks
  - **Bottom row**: Mean signal plots (one per track)
    - Shows average signal across all regions at each position
    - X-axis labels show relative position (centered at 0)
    - Negative tracks are inverted in the plot

### How It Works

1. **Load tracks**: Loads all track annotation arrays
2. **Normalize**: Determines color scale using percentile-based method:
   - For each track, finds the `per_track_max_percentile` percentile value
   - Takes the `vmax_percentile` of these maximum values as the global vmax
3. **Sort regions**: Sorts by maximum signal across all concatenated tracks
4. **Generate plots**: Creates heatmaps and mean signal plots for each track

### Example

Single track:

```bash
GenomicElementTool.py export Heatmap \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --opath heatmap.png \
    --track_npy plus_strand_track.npy \
    --title "Plus Strand" \
    --negative False
```

Paired tracks (plus and minus strand):

```bash
GenomicElementTool.py export Heatmap \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --opath heatmap.png \
    --track_npy plus_track.npy --track_npy minus_track.npy \
    --title "Plus Strand" --title "Minus Strand" \
    --negative False --negative True
```

With custom normalization:

```bash
GenomicElementTool.py export Heatmap \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --opath heatmap.png \
    --track_npy track.npy \
    --title "Signal" \
    --negative False \
    --per_track_max_percentile 95 \
    --vmax_percentile 75
```

## ChromFilteredGE

Filter genomic regions to keep only those from specified chromosomes and export as a BED file.

### Usage

```bash
GenomicElementTool.py export ChromFilteredGE [OPTIONS]
```

### Required Arguments

- `--region_file_path` (str)
  - Path to the region file (BED format)
  - Required: Yes

- `--region_file_type` (str)
  - Type of the region file
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, etc.

- `--chrom_size` (str)
  - Path to chromosome size file (tab-separated, no header)
  - Required: Yes
  - Format: Two columns: `chromosome_name<TAB>size`
  - Only regions from chromosomes listed in this file will be kept

- `--opath` (str)
  - Output path for the filtered BED file
  - Required: Yes

### Output

- **BED file**: Contains only regions from chromosomes present in the chromosome size file
  - Preserves the original BED format and all columns
  - Regions from chromosomes not in the size file are filtered out

### Example

```bash
GenomicElementTool.py export ChromFilteredGE \
    --region_file_path all_regions.bed6 \
    --region_file_type bed6 \
    --chrom_size hg38.chrom.sizes \
    --opath filtered_regions.bed6
```

Chromosome size file format (`hg38.chrom.sizes`):
```
chr1    248956422
chr2    242193529
chr3    198295559
...
```

## TREbed

Examine the GROcap/PROcap signal track to annotate regions
with both fwdTSS and revTSS. The output file is a TREbed
file as defined by `PINTS`.

### Usage

```bash
GenomicElementTool.py export TREbed [OPTIONS]
```

### Required Arguments

- `--region_file_path` (str)
  - Path to the region file (BED format)
  - Required: Yes

- `--region_file_type` (str)
  - Type of the region file
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, etc.

- `--pl_sig_track` (str)
  - path to plus strand GROcap/PROcap __signal track__ 
    npy/npz file.
  - Required: Yes

- `--mn_sig_track` (str)
  - path to minus strand GROcap/PROcap __signal track__ 
    npy/npz file.
  - Required: Yes

- `--opath` (str)
  - Output path for the TREbed file
  - Required: Yes

### Output

- **TREbed file**: Contains regions annotated with forward and reverse TSS positions
  - Format follows the TREbed specification as defined by `PINTS`
  - Each region includes `fwdTSS` and `revTSS` annotations based on GROcap/PROcap signal peaks
  - Preserves original region coordinates with added TSS annotations

### Example

```bash
GenomicElementTool.py export TREbed \
    --region_file_path all_regions.bed6 \
    --region_file_type bed6 \
    --pl_sig_track all_regions.GROcap.pl.npz \
    --mn_sig_track all_regions.GROcap.mn.npz \
    --opath all_regions.TREbed
```

The pl and mn __signal tracks__ can be obtained using 
`GenomicElementTool.py count_paired_bw`.

## Common Patterns

### Loading Annotations Before Export

Annotations (stat or track arrays) must be created using other GenomicElementTool subcommands before exporting:

```bash
# Step 1: Count signals
GenomicElementTool.py count_single_bw \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --bw_file signal.bw \
    --opath signal_track.npy

# Step 2: Export as heatmap
GenomicElementTool.py export Heatmap \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --opath heatmap.png \
    --track_npy signal_track.npy \
    --title "ChIP-seq Signal" \
    --negative False
```

### Working with Multiple Samples when exporting CountTable

For count tables with multiple samples:

```bash
GenomicElementTool.py export CountTable \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --opath all_samples.csv \
    --sample_name control1 --sample_name control2 \
    --sample_name treatment1 --sample_name treatment2 \
    --stat_npy control1.npy --stat_npy control2.npy \
    --stat_npy treatment1.npy --stat_npy treatment2.npy
```

## Dependencies

- numpy
- pandas
- matplotlib
- RGTools.GenomicElements
- RGTools.utils

## Notes

- **Region order**: All export formats preserve the order of regions as they appear in the input file
- **Coordinate convention**: All coordinates follow BED format (0-based, half-open)
- **Annotation arrays**: 
  - Stat arrays must have shape `(num_regions,)`
  - Track arrays must have shape `(num_regions, region_length)`
- **List arguments**: When using `--append` arguments (like `--sample_name`, `--track_npy`), ensure all lists have matching lengths
- **File formats**: 
  - Count tables are saved as CSV (comma-separated)
  - Heatmaps can be saved in any format supported by matplotlib (PNG, PDF, SVG, etc.)
  - Filtered regions maintain the original BED format

