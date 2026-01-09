---
title: count_single_bw and count_paired_bw Subcommands
description: Quantify signal from bigwig files across genomic regions
---

# Signal Quantification Subcommands

These subcommands allow for the quantification of signal from bigwig files across specified genomic regions. There are two versions: `count_single_bw` for a single bigwig file and `count_paired_bw` for paired (stranded) bigwig files.

## count_single_bw Subcommand

Quantify signal from a single bigwig file.

### Usage

```bash
GenomicElementTool.py count_single_bw [OPTIONS]
```

### Required Arguments

- `--bw_path` (str)
  - Path to the bigwig file
  - Required: Yes

- `--region_file_path` (str)
  - Path to the region file (BED format)
  - Required: Yes

- `--region_file_type` (str)
  - Type of the region file
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, `bedTRE`, etc.
  - See `GenomicElements` documentation for full list

- `--opath` (str)
  - Output path for the counting results (NumPy array file)
  - Required: Yes

### Optional Arguments

- `--quantification_type` (str)
  - Type of quantification to perform
  - Default: `"raw_count"`
  - Choices:
    - `raw_count`: Sum of signal in the region
    - `RPK`: Reads Per Kilobase (sum of signal / region length * 1000)
    - `full_track`: Returns the raw signal track for each region

## count_paired_bw Subcommand

Quantify signal from paired plus and minus strand bigwig files.

### Usage

```bash
GenomicElementTool.py count_paired_bw [OPTIONS]
```

### Required Arguments

- `--bw_pl` (str)
  - Plus strand bigwig file path
  - Required: Yes

- `--bw_mn` (str)
  - Minus strand bigwig file path
  - Required: Yes

- `--negative_mn` (bool)
  - Whether to output the minus strand signal as negative
  - Required: Yes

- `--flip_mn` (bool)
  - Whether to flip the minus strand signal
  - Required: Yes

- `--region_file_path` (str)
  - Path to the region file (BED format)
  - Required: Yes

- `--region_file_type` (str)
  - Type of the region file
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, `bedTRE`, etc.

- `--opath` (str)
  - Output path for the counting results (NumPy array file)
  - Required: Yes

### Optional Arguments

- `--override_strand` (str)
  - Override the strand information in the input file. Choices: `+`, `-`, `.`
  - Default: `None` (uses input strand info)

- `--quantification_type` (str)
  - Type of quantification to perform
  - Default: `"raw_count"`
  - Choices: `raw_count`, `RPK`, `full_track`

## Output

Both subcommands generate a NumPy array file at the specified `--opath`. The format depends on the file extension:
- `.npy`: Standard NumPy array file.
- `.npz`: Compressed NumPy archive (saved as the default array name).

- **Shape**: `(num_regions, 1)` for `raw_count` and `RPK`, 
    or `(num_regions, max(region_length))` for `full_track`.
- **Content**: Quantified signal for each region in the input file.

The results can be loaded later using `GenomicElements.load_region_anno_from_npy()`.

## Examples

### Basic Single BigWig Quantification

```bash
GenomicElementTool.py count_single_bw \
    --bw_path signal.bw \
    --region_file_path regions.bed \
    --region_file_type bed3 \
    --opath counts.npy
```

### Paired BigWig Full Track Extraction

```bash
GenomicElementTool.py count_paired_bw \
    --bw_pl plus.bw \
    --bw_mn minus.bw \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --quantification_type full_track \
    --negative_mn True \
    --flip_mn False \
    --opath signals.npy
```

### Paired BigWig extraction that write the mn signal from 5' to 3'

```bash
GenomicElementTool.py count_paired_bw \
    --bw_pl plus.bw \
    --bw_mn minus.bw \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --quantification_type full_track \
    --negative_mn True \
    --flip_mn True \
    --opath signals.npy
```

## Dependencies

- numpy
- pyBigWig
- RGTools.GenomicElements
- RGTools.BwTrack

## Notes

### `--negative_mn` and `--flip_mn` in `count_paired_bw`

When quantifying paired pl and mn bigWigs, the user have 
the option to flip the output (write signal track in 5' - 3' 
order) and to negate the output (so that mn strand signal 
are stored as negative values). This is irrelevant 
to the negativity of signal stored in bigWig files since 
absolute values are taken when reading them. For quantification 
type that is flip-invariant (e.g. `raw_count` and `rpk`), 
The `--flip_mn` flag will not alter the output and `--negative_mn` 
will negate it.

