---
title: mask_op Subcommand
description: Operation for masks of the same Genomic Element
---

# `mask_op` Subcommand

The `mask_op` subcommand performs logical operations on mask
annotations from the same Genomic Element dataset.
It is useful for combining filters generated from different
criteria (for example, signal threshold, motif presence, or QC filters).

## Usage

```bash
GenomicElementTool.py mask_op <operation> [OPTIONS]
```

The `mask_op` subcommand supports the following operations:

- `intersect`: element-wise AND (keeps entries that are `True` in all masks)
- `union`: element-wise OR (keeps entries that are `True` in any mask)
- `opposite`: element-wise NOT (inverts `True`/`False` values of one mask)

## `intersect`

Compute the logical intersection of all provided masks.

### Usage

```bash
GenomicElementTool.py mask_op intersect [OPTIONS]
```

### Required Arguments

- `--region_file_path` (str)
  - Path to the input region file corresponding to the mask arrays.
  - Required: Yes

- `--region_file_type` (str)
  - Type of the input region file.
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, `bed3gene`, `narrowPeak`, `TREbed`, `bedGraph`

- `--mask_npy` (str)
  - Path to an input mask array file (`.npy` or single-array `.npz`)
  - Required: Yes
  - Can be specified multiple times (minimum 2)

- `--opath` (str)
  - Output path of the intersected mask (`.npy`)
  - Required: Yes

### Behavior

- element-wise AND

### Example

```bash
GenomicElementTool.py mask_op intersect \
  --region_file_path regions.bed3 \
  --region_file_type bed3 \
  --mask_npy accessible.npy \
  --mask_npy motif_present.npy \
  --opath strict_filter.npy
```

## `union`

Compute the logical union of all provided masks.

### Usage

```bash
GenomicElementTool.py mask_op union [OPTIONS]
```

### Required Arguments

- `--region_file_path` (str)
  - Path to the input region file corresponding to the mask arrays.
  - Required: Yes

- `--region_file_type` (str)
  - Type of the input region file.
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, `bed3gene`, `narrowPeak`, `TREbed`, `bedGraph`

- `--mask_npy` (str)
  - Path to an input mask array file (`.npy` or single-array `.npz`)
  - Required: Yes
  - Can be specified multiple times (minimum 2)

- `--opath` (str)
  - Output path of the unioned mask (`.npy`)
  - Required: Yes

### Behavior

- element-wise OR

### Example

```bash
GenomicElementTool.py mask_op union \
  --region_file_path regions.bed3 \
  --region_file_type bed3 \
  --mask_npy active_in_celltypeA.npy \
  --mask_npy active_in_celltypeB.npy \
  --opath combined_activity.npy
```

## `opposite`

Compute the logical opposite of one input mask.

### Usage

```bash
GenomicElementTool.py mask_op opposite [OPTIONS]
```

### Required Arguments

- `--region_file_path` (str)
  - Path to the input region file corresponding to the mask array.
  - Required: Yes

- `--region_file_type` (str)
  - Type of the input region file.
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, `bed3gene`, `narrowPeak`, `TREbed`, `bedGraph`

- `--mask_npy` (str)
  - Path to an input mask array file (`.npy` or single-array `.npz`)
  - Required: Yes
  - Must be specified exactly once

- `--opath` (str)
  - Output path of the inverted mask (`.npy`)
  - Required: Yes

### Behavior

- element-wise NOT

### Example

```bash
GenomicElementTool.py mask_op opposite \
  --region_file_path regions.bed3 \
  --region_file_type bed3 \
  --mask_npy blacklisted_regions.npy \
  --opath keep_regions.npy
```

## Output

- **Mask array file** (`.npy`)
  - `mask` annotation type for Genomic Elements
  - Contains one boolean value per region (commonly stored as `(num_regions,)` or `(num_regions, 1)`)

## Typical Workflows

Combine two independent filters:

```bash
GenomicElementTool.py mask_op intersect \
  --region_file_path regions.bed6 \
  --region_file_type bed6 \
  --mask_npy pass_qc.npy \
  --mask_npy has_signal.npy \
  --opath final_mask.npy
```

Then export filtered regions:

```bash
GenomicElementTool.py export MaskedGE \
  --region_file_path regions.bed6 \
  --region_file_type bed6 \
  --mask_npy final_mask.npy \
  --opath filtered_regions.bed6
```