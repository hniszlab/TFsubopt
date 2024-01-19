#!/usr/bin/env bash


computeMatrix reference-point \
--referencePoint center \
-b 1500 -a 1500 \
--binSize 10 \
--regionsFileName ./bed/OUT/ZIP13K2_NGN2_24h_WT_merged.bed ./bed/OUT/ZIP13K2_NGN2_48h_WT_merged.bed ./bed/OUT/ZIP13K2_NGN2_24h_AroPERFECT_merged.bed ./bed/OUT/ZIP13K2_NGN2_48h_AroPERFECT_merged.bed ./bed/OUT/ZIP13K2_NGN2_24h_AroLITE_merged.bed ./bed/OUT/ZIP13K2_NGN2_48h_AroLITE_merged.bed \
--scoreFileName ./bw/ZIP13K2_NGN2_24h_WT_merged.bw ./bw/ZIP13K2_NGN2_48h_WT_merged.bw ./bw/ZIP13K2_NGN2_24h_WT_merged.bw ./bw/ZIP13K2_NGN2_48h_AroPERFECT_merged.bw ./bw/ZIP13K2_NGN2_24h_AroLITE_merged.bw ./bw/ZIP13K2_NGN2_48h_AroLITE_merged.bw \
-p 8 \
--smartLabels \
-o NGN2_RefPoint_cM_b10_1500bp_sorted.gz \
--outFileSortedRegions NGN2_RefPoint_cM_b10_1500bp_sorted.bed

plotHeatmap \
--matrixFile NGN2_RefPoint_cM_b10_1500bp_sorted.gz \
--plotFileFormat pdf \
--colorList 'white,#f05756' 'white,#f05756' 'white,#652d90' 'white,#652d90' 'white,#1f72b8' 'white,#1f72b8'\
--missingDataColor 1 \
-out NGN2_RefPoint_cM_b10_1500bp_sorted.pdf \
--dpi 300 \
--sortUsingSamples 3 4
