---
title: pad_region Subcommand
description: Expand or shrink genomic regions by a specified amount
---

# pad_region Subcommand

The `pad_region` subcommand allows for expanding or shrinking genomic regions in a BED file. Padding is performed relative to the strand of each region unless specified otherwise.

## Usage

```bash
GenomicElementTool.py pad_region [OPTIONS]
```

## Required Arguments

- `--region_file_path` (str)
  - Path to the input region file (e.g., BED format).
  - Required: Yes

- `--region_file_type` (str)
  - Type of the input region file.
  - Required: Yes
  - Valid types: `bed3`, `bed6`, `bed6gene`, `bedTRE`, etc.
  - See `GenomicElements` documentation for the full list of supported types.

- `--upstream_pad` (int)
  - Amount to extend to the upstream of the region.
  - Positive values expand the region, while negative values shrink it.
  - Required: Yes

- `--downstream_pad` (int)
  - Amount to extend to the downstream of the region.
  - Positive values expand the region, while negative values shrink it.
  - Required: Yes

- `--opath` (str)
  - Output path for the padded BED file.
  - Required: Yes

## Optional Arguments

- `--ignore_strand` (bool)
  - Whether to ignore the strand information in the input file.
  - If `True`, or if the input file does not contain strand information (e.g., `bed3`, `TREbed`, `bedGraph`), all regions are treated as being on the positive (`+`) strand. In this case, upstream is always towards the smaller coordinate.
  - **Note**: For file types without a strand column (like `bed3`), you **must** set this to `True`; otherwise, the tool will raise an error.
  - Default: `False`

- `--method_resolving_invalid_region` (str)
  - Method to resolve regions that become invalid after padding (e.g., start > end).
  - Choices:
    - `fallback`: Use the original region if the padded region is invalid.
    - `raise`: Raise an `InvalidBedRegionException`.
    - `drop`: Remove the region from the output.
  - Default: `fallback`

## Examples

### Expand regions by 100bp on both sides

```bash
GenomicElementTool.py pad_region \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --upstream_pad 100 \
    --downstream_pad 100 \
    --opath regions.padded.bed6
```

### Shrink regions by 50bp from the 5' end (upstream)

```bash
GenomicElementTool.py pad_region \
    --region_file_path regions.bed6 \
    --region_file_type bed6 \
    --upstream_pad -50 \
    --downstream_pad 0 \
    --opath regions.shrunk.bed6
```

### Expand regions ignoring strand information

```bash
GenomicElementTool.py pad_region \
    --region_file_path regions.bed3 \
    --region_file_type bed3 \
    --upstream_pad 100 \
    --downstream_pad 100 \
    --ignore_strand True \
    --opath regions.padded.bed3
```

## Dependencies

- pandas
- RGTools.GenomicElements
- RGTools.exceptions
- RGTools.utils

## Notes

### Handling Files Without Strand Information

For genomic region files that do not contain strand information (such as `bed3`, `TREbed`, `bedGraph`, or `bed3gene`), the `pad_region` subcommand requires the `--ignore_strand True` flag. 

If `--ignore_strand` is set to `False` (the default) for these file types, the tool will raise an `InvalidStrandnessException` because it cannot find a `strand` field to determine the upstream and downstream directions. When `--ignore_strand True` is used, all regions are treated as being on the positive (`+`) strand, meaning:
- **Upstream** padding extends towards smaller coordinates.
- **Downstream** padding extends towards larger coordinates.

