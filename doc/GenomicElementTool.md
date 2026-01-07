---
title: GenomicElementTool
description: Command-line tool for working with genomic elements
---

# GenomicElementTool

GenomicElementTool is a command-line tool for working with genomic elements. It provides a suite of subcommands for various genomic data processing tasks including signal counting, region manipulation, sequence encoding, motif searching, and data import/export.

## Overview

GenomicElementTool is designed to work with the Genomic Element data structure, which provides a standardized way to store and manipulate genomic regions along with their associated annotations and sequences.

## Usage

```bash
GenomicElementTool.py <subcommand> [OPTIONS]
```

To see available subcommands and options:

```bash
GenomicElementTool.py --help
```

To see help for a specific subcommand:

```bash
GenomicElementTool.py <subcommand> --help
```

## Genomic Element Data Structure

Genomic Element is the key data structure this tool is based on. A Genomic Element contains 3 key components:

### Reference Genome File

The reference genome file of a Genomic Element is in FASTA format. This file contains the reference genome sequences used to extract sequences for the genomic regions.

### Region File

A region file defines the coordinates of the "elements". Region files are one of the supported variants of BED format. The coordinates follow the same half-open convention as BED files (0-based, end-exclusive).

Supported region file types include:
- `bed3`: Standard 3-column BED format (chrom, start, end)
- `bed6`: Standard 6-column BED format (chrom, start, end, name, score, strand)
- `bed6gene`: BED6 format with gene-specific annotations
- `bed3gene`: BED3 format with gene-specific annotations
- `narrowPeak`: ENCODE narrowPeak format
- `TREbed`: Transcription regulatory element BED format
- `bedGraph`: BEDGraph format

### Annotations

Annotations are NumPy arrays stored in `.npy` files. The first dimension of the array always equals the number of regions. Depending on the second dimension size, an annotation can be of different types:

- **track**: Store signal tracks for each element. The second dimension is larger than 1. These are used to store continuous signals across regions (e.g., ChIP-seq signals, motif scores, etc.).
- **stat**: Store statistics for each element. The second dimension is of size 1. These are used to store single values per region (e.g., expression levels, peak scores, etc.).

## Subcommands

GenomicElementTool provides the following subcommands:

### count_single_bw

Count signal in a single BigWig file.

```bash
GenomicElementTool.py count_single_bw [OPTIONS]
```

Extracts signal values from a BigWig track for specified genomic regions.

### count_paired_bw

Count signal in paired BigWig files.

```bash
GenomicElementTool.py count_paired_bw [OPTIONS]
```

Extracts signal values from a pair of BigWig tracks (e.g., plus and minus strand tracks) for specified genomic regions.

### pad_region

Pad regions while conserving the order of elements in Genomic Elements files.

```bash
GenomicElementTool.py pad_region [OPTIONS]
```

This program differs from standard padding tools in that it preserves the order of elements in Genomic Elements files.

### bed2tssbed

Convert BED file to TSS (Transcriptional Start Site) BED file.

```bash
GenomicElementTool.py bed2tssbed [OPTIONS]
```

Converts a BED file containing gene regions to a BED file containing TSS coordinates.

### onehot

One-hot encode sequences. Only supports elements of the same size.

```bash
GenomicElementTool.py onehot [OPTIONS]
```

Converts DNA sequences to one-hot encoded representations. All regions must have the same length.

### motif_search

Search for motifs in sequences.

```bash
GenomicElementTool.py motif_search [OPTIONS]
```

Searches for motifs in genomic sequences using Position Weight Matrices (PWMs) from MEME format files. See [motif_search.md](motif_search.md) for detailed documentation.

### track2tss_bed

Produce TSS BED file from track.

```bash
GenomicElementTool.py track2tss_bed [OPTIONS]
```

Generates TSS coordinates from signal tracks (e.g., identifying TSS locations based on signal peaks).

### filter_motif_score

Filter motif search scores.

```bash
GenomicElementTool.py filter_motif_score [OPTIONS]
```

Filters and processes motif search score outputs, typically used to identify significant motif matches.

### export

Export to other formats.

```bash
GenomicElementTool.py export [OPTIONS]
```

Exports Genomic Element data to various output formats.

### import

Import from other formats.

```bash
GenomicElementTool.py import [OPTIONS]
```

Imports data from various formats into the Genomic Element structure.

## Getting Help

For detailed help on any subcommand:

```bash
GenomicElementTool.py <subcommand> --help
```

For example:

```bash
GenomicElementTool.py motif_search --help
GenomicElementTool.py count_single_bw --help
```

## Examples

### Basic workflow

1. Import data into Genomic Element format:
   ```bash
   GenomicElementTool.py import --input_file data.bed6 --region_file_type bed6
   ```

2. Count signals from BigWig files:
   ```bash
   GenomicElementTool.py count_single_bw --bw_file signal.bw --region_file regions.bed6
   ```

3. Search for motifs:
   ```bash
   GenomicElementTool.py motif_search --fasta_path genome.fa --region_file_path regions.bed6 --region_file_type bed6 --motif_file motifs.meme
   ```

4. Export results:
   ```bash
   GenomicElementTool.py export --output_file results.txt
   ```

## Dependencies

GenomicElementTool relies on the RGTools library and requires:
- Python 3.x
- NumPy
- BioPython (for sequence handling)
- pyBigWig (for BigWig file support)
- RGTools package (GenomicElements, BedTable, BwTrack, MemeMotif, etc.)

## Notes

- All genomic coordinates follow the BED convention: 0-based, end-exclusive (half-open intervals)
- Region files must be sorted or compatible with the expected format
- Annotation arrays are stored in NumPy's `.npy` format for efficient I/O
- The tool preserves element order when processing, which is important for maintaining correspondence between regions and annotations


