computeMatrix reference-point \
--referencePoint center \
-b 1500 -a 1500 \
--binSize 10 \
--sortRegions descend \
--sortUsing max \
--skipZeros \
--regionsFileName CEBPa_WT_IS15_Common.bed CEBPa_IS15_unique_filtered_dectedtedBefore_wa.bed CEBPa_IS15_unique_filtered_neverDetected_FilteredCEBPa_F.bed \
--scoreFileName WT_24_ChIP_merged.bw WT_48_ChIP_merged.bw IS15_24_ChIP_merged.bw IS15_48_ChIP_merged.bw \
-p 30 \
--smartLabels \
-o CEBPa_cM_b10_1500bp.gz \
--outFileSortedRegions CEBPa_cM_b10_1500bp_sorted.bed



plotHeatmap \
--matrixFile CEBPa_cM_b10_1500bp.gz \
--plotFileFormat png \
--colorList 'white,#f05756' 'white,#f05756' 'white,#652d90' 'white,#652d90' \
--missingDataColor 1 \
-out CEBPa_RefPoint_cM_b10_1500bp.png \
--dpi 300
