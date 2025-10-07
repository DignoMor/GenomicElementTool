
OPATH="heatmap_example"
mkdir -p ${OPATH}

cat RGTools/large_files/ENCFF156JSS.bedTRE | \
    awk 'BEGIN{FS="\t"; OFS="\t"}{print $1, $2, $3}' > \
    ${OPATH}/pints_regions.bed

./GenomicElementTool.py bed2tssbed \
    --region_file_path ${OPATH}/pints_regions.bed \
    --region_file_type bed3 \
    --output_site center \
    --opath ${OPATH}/pints_regions.center.bed3

./GenomicElementTool.py pad_region \
    --region_file_path ${OPATH}/pints_regions.center.bed3 \
    --region_file_type bed3 \
    --upstream_pad 500 \
    --downstream_pad 499 \
    --ignore_strand True \
    --method_resolving_invalid_region drop \
    --opath ${OPATH}/pints_regions.1kb.bed3

./GenomicElementTool.py count_single_bw \
    --region_file_path ${OPATH}/pints_regions.1kb.bed3 \
    --region_file_type bed3 \
    --bw_path RGTools/large_files/ENCFF565BWR.pl.bw \
    --quantification_type full_track \
    --opath ${OPATH}/pints_regions.1kb.pl.track.npy

./GenomicElementTool.py count_single_bw \
    --region_file_path ${OPATH}/pints_regions.1kb.bed3 \
    --region_file_type bed3 \
    --bw_path RGTools/large_files/ENCFF775FNU.mn.bw \
    --quantification_type full_track \
    --opath ${OPATH}/pints_regions.1kb.mn.track.npy

./GenomicElementTool.py export Heatmap \
    --region_file_path ${OPATH}/pints_regions.1kb.bed3 \
    --region_file_type bed3 \
    --track_npy ${OPATH}/pints_regions.1kb.pl.track.npy \
    --title "GROcap.pl" \
    --per_track_max_percentile 99 \
    --vmax_percentile 20 \
    --opath ${OPATH}/pints_regions.1kb.pl.heatmap.png

./GenomicElementTool.py export Heatmap \
    --region_file_path ${OPATH}/pints_regions.1kb.bed3 \
    --region_file_type bed3 \
    --title "GROcap.mn" \
    --track_npy ${OPATH}/pints_regions.1kb.mn.track.npy \
    --negative True \
    --opath ${OPATH}/pints_regions.1kb.mn.heatmap.png
