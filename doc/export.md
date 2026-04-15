---
title: export Subcommand
description: Export Genomic Element data to various output formats
---

# export Subcommand

The `export` subcommand provides functionality to export Genomic Element data to various output formats including FASTA sequences, count tables, heatmaps, filtered region files, and merged genomic elements.

## Usage

```bash
GenomicElementTool.py export <oformat> [OPTIONS]
```

The export subcommand supports the following output formats:
- `ExogeneousSequences`: Export genomic sequences to Exogeneous Sequences (FASTA format).
- `WTES`: Export genomic sequences to Wild Type ES ready for oligo library construction pipeline (FASTA format).
- `CountTable`: Export statistical annotations as a count table (CSV)
- `Heatmap`: Generate heatmap visualizations from track annotations
- `ChromFilteredGE`: Filter regions by chromosome and export as BED file
- `MaskedGE`: Filter regions by a mask annotation and export a filtered Genomic Element dataset
- `TREbed`: Annotate regions with forward and reverse TSS from GROcap/PROcap signals
- `MergedGE`: Merge multiple Genomic Element files into one
- `bed6poly`: Bed6 file with an extra column for polymorphism information (e.g., SNPs).

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

- `--opath` (str)
  - Fasta output file path.
  - Required: Yes

### Output

- **FASTA file**: `<opath>`
  - Contains one sequence per genomic region
  - Sequence headers use format: `>chrom:start-end`
  - Sequences are extracted from the reference genome based on region coordinates

### Example

```bash
GenomicElementTool.py export ExogeneousSequences \
    --fasta_path /path/to/genome.fa \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --opath my_sequences.fa
```

This writes `my_sequences.fa` with sequences for each region in the BED file.

## WTES

Export genomic regions as sequences in FASTA format for Wild Type ES ready for oligo library construction pipeline.

### Usage

```bash
GenomicElementTool.py export WTES [OPTIONS]
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

- `--num_replicates` (int)
  - Number of replicates for the Wild Type ES
  - Required: Yes

- `--opath` (str)
  - Output path for the FASTA file
  - Required: Yes

### Output

- **FASTA file**: `opath`
  - Contains `num_replicates` sequences per genomic region
  - Sequence headers use format: `>chrom:start-end_<replicate_number>` for each replicate
  - Sequences are extracted from the reference genome based on region coordinates
  - For WT elements, replicates simply means multiple copies of the same sequence.

### Example

```bash
GenomicElementTool.py export WTES \
    --fasta_path /path/to/genome.fa \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --num_replicates 5 \
    --opath my_wtes.fa
```

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

## MaskedGE

Filter genomic regions using a boolean mask annotation and export the selected regions as a BED file. Optionally, export one or more annotations sliced by the same mask.

### Usage

```bash
GenomicElementTool.py export MaskedGE [OPTIONS]
```

### Required Arguments

- `--region_file_path` (str)
  - Path to the region file (BED format)
  - Required: Yes

- `--region_file_type` (str)
  - Type of the region file
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, etc.

- `--mask_npy` (str)
  - Path to a mask annotation NumPy file (`.npy` or `.npz`)
  - Required: Yes
  - Must contain one boolean value per input region
  - `True` values are kept in the output

- `--opath` (str)
  - Output path for the filtered BED file
  - Required: Yes

### Optional Arguments

- `--anno_name` (str)
  - Annotation name to load and export after masking
  - Can be specified multiple times
  - Must align with `--anno_npy` and `--anno_type`

- `--anno_npy` (str)
  - Path to annotation file (`.npy` or single-array `.npz`)
  - Can be specified multiple times
  - Must align with `--anno_name` and `--anno_type`

- `--anno_type` (str)
  - Annotation type for each input annotation
  - Can be specified multiple times
  - Choices: `track`, `stat`, `mask`, `array`
  - Must align with `--anno_name` and `--anno_npy`

- `--anno_oheader` (str)
  - Output prefix for masked annotation files
  - Required when any annotation arguments are provided
  - Output files are written as `<anno_oheader>.<anno_name>.npy`

### Output

- **Output Region file**: Contains only regions where the mask annotation is `True`
  - Preserves the original BED format and all columns
  - Keeps the original region order among retained regions
- **Optional masked annotation files**:
  - One `.npy` file per requested annotation
  - Annotation arrays are sliced by the same mask used for regions

### Example

```bash
GenomicElementTool.py export MaskedGE \
    --region_file_path all_regions.bed6 \
    --region_file_type bed6 \
    --mask_npy high_confidence_mask.npy \
    --opath high_confidence_regions.bed6
```

With annotation export:

```bash
GenomicElementTool.py export MaskedGE \
    --region_file_path all_regions.bed6 \
    --region_file_type bed6 \
    --mask_npy high_confidence_mask.npy \
    --opath high_confidence_regions.bed6 \
    --anno_name sample_stat --anno_npy sample_stat.npy --anno_type stat \
    --anno_name sample_track --anno_npy sample_track.npy --anno_type track \
    --anno_oheader high_confidence_regions
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

## MergedGE

Merge two Genomic Element DataSet into a single 
Genomic Element DataSet. 
This combines the regions from both input files 
into one output Region file.

### Usage

```bash
GenomicElementTool.py export MergedGE [OPTIONS]
```

### Required Arguments

- `--left_region_file_path` (str)
  - Path to the first region file.
  - Required: Yes

- `--right_region_file_path` (str)
  - Path to the second region file.
  - Required: Yes

- `--region_file_type` (str)
  - Type of the region files.
  - Required: Yes

- `--anno_name` (str)
  - Name of the annotation tracks to be merged.
  - Can be multiple.

- `--left_anno_path` (str)
  - Path to the annotation track of the first GE.
  - Can be multiple, must match the number of `anno_name`.

- `--right_anno_path` (str)
  - Path to the annotation track of the second GE.
  - Can be multiple, must match the number of `anno_name`.

- `--oheader` (str)
  - Output header for the region file and annotations.
  - Required: Yes

### Output

- **Region file**: A single region file containing all regions 
  from the input files and sorted by coordinates. File type 
  is determined by `region_file_type`.

### Example

```bash
GenomicElementTool.py export MergedGE \
    --left_region_file_path file1.bed6 \
    --right_region_file_path file2.bed6 \
    --anno_name example_anno \
    --left_anno_path file1.anno.npy \
    --right_anno_path file2.anno.npy \
    --region_file_type bed6 \
    --oheader merged
```

## bed6poly

Export regions as `bed6poly`, a BED6-plus file with an additional polymorphism column 
seperated by `/`.

The program utilizes `RGTools.SNP_utils` for ensembl API to get the polymorphism information from rsids. If the rsid is not found, the program will raise an exception
or drop the region. 
However, if the position does not match, the program will not stop but throw a warning.


### Usage

```bash
GenomicElementTool.py export bed6poly [OPTIONS]
```

### Required Arguments

- `--region_file_path` (str)
  - Path to the input region file (BED format)
  - Required: Yes

- `--region_file_type` (str)
  - Type of the region file
  - Required: Yes
  - Valid types: `bed6` only.

- `--opath` (str)
  - Output path for the `bed6poly` file
  - Required: Yes

### Optional Arguments

- `--genome_version` (str)
  - Genome assembly used for SNP lookup through Ensembl REST API
  - Default: `hg38`
  - Valid values: `hg38`, `GRCh38`, `hg19`, `GRCh37`

- `--rsid_not_found_handling` (str)
  - Handling of rsids not found in Ensembl REST API
  - Default: `raise`
  - Valid values: `raise`, `drop`

### Output

- **bed6poly file**: Contains one row per input region with:
  - BED6 core columns: `chrom`, `start`, `end`, `name`, `score`, `strand`
  - Additional polymorphism column with region-level polymorphism summary, 
    seperated by `/`.
  - Alleles follow Ensembl `allele_string` order: reference allele first, then alternate allele(s). The order of alternative alleles is not a priority ranking.
  - Original region ordering preserved

Example output:
```
chr1	161133654	161133655	rs10797093	2.677e-08	+	T/G
chr1	161136563	161136564	rs11265557	2.242e-08	+	T/C/G
chr1	161140556	161140557	rs12041364	2.206e-08	+	G/A/C
chr1	161141572	161141573	rs11591206	2.326e-08	+	C/A/G/T
```

### Example

```bash
GenomicElementTool.py export bed6poly \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --genome_version hg38 \
    --opath regions.bed6poly
```

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

