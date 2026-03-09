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

### Example

```bash
GeomicElementTool.py get_context_ge nearest \
  --region_file_path enhancer.bed3 \
  --region_file_type bed3 \
  --context_file_path promoter.bed3 \
  --context_file_type bed3 \
  --opath per_enhancer_context.bed3

```

