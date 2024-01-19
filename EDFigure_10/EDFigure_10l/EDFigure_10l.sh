computeMatrix scale-regions -S ZIP13K2_NGN2_AroLITE_12h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_AroLITE_12h_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_AroLITE_24h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_AroLITE_24h_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_AroPERFECT_12_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_AroPERFECT_12h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_AroPERFECT_24h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_AroPERFECT_24h_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_WT_12h_r1_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_WT_12h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_WT_24h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_WT_24h_TCReads_cov_CPM_reverse.bw ZIP13K2_TCReads_cov_CPM_forward.bw ZIP13K2_TCReads_cov_CPM_reverse.bw -R hg38_genes_refseq_plus.bed -b 2000 -a 2000 --skipZeros -o matrix_ttslamseq_scaled2k_plus.gz --outFileNameMatrix matrix_ttslamseq_scaled2k_plus.tab --outFileSortedRegions matrix_ttslamseq_scaled2k_plus.bed -p 80 -m 2000

plotProfile -m matrix_ttslamseq_scaled2k_plus.gz --plotFileFormat svg -o out_2k_plus.svg
plotProfile -m matrix_ttslamseq_scaled2k_minus.gz --plotFileFormat svg -o out_2k_minus.svg



computeMatrix scale-regions -S ZIP13K2_NGN2_AroLITE_12h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_AroLITE_12h_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_AroLITE_24h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_AroLITE_24h_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_AroPERFECT_12_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_AroPERFECT_12h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_AroPERFECT_24h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_AroPERFECT_24h_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_WT_12h_r1_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_WT_12h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_WT_24h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_WT_24h_TCReads_cov_CPM_reverse.bw ZIP13K2_TCReads_cov_CPM_forward.bw ZIP13K2_TCReads_cov_CPM_reverse.bw -R hg38_genes_minus.bed -b 2000 -a 2000 --skipZeros -o matrix_ttslamseq_scaled_minus.gz --outFileNameMatrix matrix_ttslamseq_scaled_minus.tab --outFileSortedRegions matrix_ttslamseq_scaled_minus2.bed -p 80 -m 5000

plotProfile -m matrix_ttslamseq_scaled_minus.gz --plotFileFormat svg -o out_minus.svg


computeMatrix scale-regions -S ZIP13K2_NGN2_AroLITE_12h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_AroLITE_12h_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_AroLITE_24h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_AroLITE_24h_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_AroPERFECT_12_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_AroPERFECT_12h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_AroPERFECT_24h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_AroPERFECT_24h_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_WT_12h_r1_TCReads_cov_CPM_reverse.bw ZIP13K2_NGN2_WT_12h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_WT_24h_TCReads_cov_CPM_forward.bw ZIP13K2_NGN2_WT_24h_TCReads_cov_CPM_reverse.bw ZIP13K2_TCReads_cov_CPM_forward.bw ZIP13K2_TCReads_cov_CPM_reverse.bw -R hg38_genes_plus.bed -b 2000 -a 2000 --skipZeros -o matrix_ttslamseq_scaled_plus.gz --outFileNameMatrix matrix_ttslamseq_scaled_plus.tab --outFileSortedRegions matrix_ttslamseq_scaled_plus.bed -p 80 -m 5000

plotProfile -m matrix_ttslamseq_scaled_plus.gz --plotFileFormat svg -o out_plus.svg
