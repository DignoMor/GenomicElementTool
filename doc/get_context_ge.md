---
title: get_context_ge Subcommand
description: Generate context genomic elements 
---

# get_context_ge Subcommand

Select a context region for each element from 
a pool of context (e.g. promoters).

## Usage

```bash
GenomicElementTool.py get_context_ge <method> [OPTIONS]
```

The following methods are supported:

- `nearest`: select context region by minimal distance
- `windowed_argmax`: select context region with maximum 
  provided stat in a window.

## `nearest`

select context region by minimal distance. Distance 
is calculated by the distance between the closest 
2 ends of the elements.

### Required Arguments

- `--region_file_path` (str)
  - Path to the region file 
  - Required: Yes

- `--region_file_type` (str)
  - Type of the region file
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, `TREbed`, etc.
  - See `GenomicElements` documentation for full list

- `--context_file_path` (str)
  - Path to the context element file 
  - Required: Yes

- `--context_file_type` (str)
  - Type of the context element file
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, `TREbed`, etc.
  - See `GenomicElements` documentation for full list

- `--opath` (str)
  - Path to the output region file.
  - Required: Yes

### Example

```bash
GenomicElementTool.py get_context_ge nearest \
  --region_file_path enhancer.bed3 \
  --region_file_type bed3 \
  --context_file_path promoter.bed3 \
  --context_file_type bed3 \
  --opath per_enhancer_context.bed3

```

## `windowed_argmax`

Select context region with the maximum provided stat in a window.

### Required Arguments

- `--region_file_path` (str)
  - Path to the region file. The file should be the 
    element of interest file padded to the desired window
    size.
  - Required: Yes

- `--region_file_type` (str)
  - Type of the region file
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, `TREbed`, etc.
  - See `GenomicElements` documentation for full list

- `--context_file_path` (str)
  - Path to the context element file 
  - Required: Yes

- `--context_file_type` (str)
  - Type of the context element file
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, `TREbed`, etc.
  - See `GenomicElements` documentation for full list

- `--context_stat_path` (str)
  - Path to the stat file for the context GE
  - Required: Yes

- `--opath` (str)
  - Path to the output region file.
  - Required: Yes

### Example

```bash
GenomicElementTool.py bed2tssbed \
  --region_file_path enhancer.bed3 \
  --region_file_type bed3 \
  --opath stdout \
  --output_site center | \
GenomicElementTool.py pad_region \
  --region_file_path stdin \
  --region_file_type bed3 \
  --upstream_pad 500000 \
  --downstream_pad 499999 \
  --opath stdout \
  --ignore_strand True | \
GenomicElementTool.py get_context_ge windowed_argmax \
  --region_file_path stdin \
  --region_file_type bed3 \
  --context_file_path promoter.bed3 \
  --context_file_type bed3 \
  --context_stat_path promoter.GROcap_count.npy \
  --opath stdout | \
GenomicElementTool.py onehot \
  --fasta_path hg38.fa \
  --region_file_path stdin \
  --region_file_type bed3 \
  --opath enhancer.context.npy
```
