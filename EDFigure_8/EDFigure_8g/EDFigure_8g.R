library(Seurat)
library(cowplot)
library(hdf5r)
library(kableExtra)
library(biomaRt)
library(limma)
library(topGO)
library(org.Hs.eg.db)
library(sva)
library(scran)
library(tidyverse)
library(DoubletFinder)
library(gprofiler2)
library(rliger)
library(future)
library(harmony)
library(pheatmap)
library(clustree)
library(venn)
library(enrichR)
library(rafalib)
library(scPred)
library(loomR)
library(SeuratWrappers)
library(velocyto.R)
library("Nebulosa")
library(ggalluvial)
library(patchwork)
plan("multicore", workers = 10)
options(future.globals.maxSize = 28 * 1024 ^ 3)
set.seed(12345)

alldataF <- readRDS("Main_seuratObject.rds")


# Compute differentiall expression
DEGS_clusters <- FindAllMarkers(alldataF, logfc.threshold = 0.2, test.use = "wilcox", 
                                min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, 
                                assay = "RNA")


top10allclustersdf <- DEGS_clusters %>% group_by(cluster) %>% top_n(-10, p_val)
write.csv(top10allclustersdf, "top10allclusters.csv")
top10allclusters <- top10allclustersdf %>% select(gene) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]


markers_cluster <- c("CD37", "IFITM3", "ARPP21", "FRMD4B", "IGFBP2", "LINC01013", "RFLNB", "IRX1", "H2AC14", "H3C2", "H3C4", "H2BC7", "H2AC11", "H2AC16", "H1-2", "H2AC20", "H2AC17", "H1-5", "H1-4", "H4C5", "H4C4", "H2BC9", "H3C8", "H2BC18", "H1-3", "H2BC12", "H3C7", "H3C3", "H3C10", "H2BC11", "CDK1", "NUF2", "H2AC12", "CHST6", "MT1G", "TRIM14", "PLS3", "KCNJ12", "CLEC2D", "NLGN1", "GPR160", "BLNK", "CRIP1", "GATM", "MAP1A", "BACH2", "CXCR4", "IARS2", "CORO2A", "ENDOD1", "DDX54", "AFF3", "PRKCZ", "BRI3BP", "RABGAP1L", "DTX1", "CFD", "AZU1", "PDE4D", "RETN", "ANXA1", "TFEC", "CSF3R", "WLS", "AIF1", "ELANE", "NUCB2", "ASPH", "SAP30", "CSTA", "BLVRB", "MS4A3", "MED30", "TRGC2", "HGF", "FUT7", "NFIL3", "ANKRD22", "GPI", "ZNF467", "RAB27A", "CD79A", "CCDC26", "SIPA1L2", "PPA1", "NRIP1", "CDK6", "VPREB1", "PSME2", "CCT3", "ZNF608", "PXDN", "UGT3A2", "LDHB", "MZB1", "OFD1", "PPFIBP1", "LRRC28", "HLA-DPB1", "MLEC", "MARCKSL1", "TOP1MT", "SRM", "BMI1", "KCNQ5", "LEF1", "CD14", "S100A9", "NCF1", "S100A8", "KCTD12", "JAML", "VNN2", "ITGB2", "VASP", "CYBB", "FCER1G", "LST1", "PTPRC", "LILRA5", "FCN1", "LILRB3", "TYROBP", "BCL2A1", "PIM1", "S100A12", "C5AR1", "LSP1", "CDC42EP3", "CD55", "MPEG1"
)

top5markets <- c("CD37", "IFITM3", "ARPP21", "FRMD4B", "IGFBP2", "H2AC14", "H3C2", "H3C4", "H2BC7", "H2AC11", "CHST6", "MT1G", "TRIM14", "PLS3", "KCNJ12", "CLEC2D
CFD", "AZU1", "PDE4D", "RETN", "ANXA1", "CD79A", "CCDC26", "SIPA1L2", "PPA1", "NRIP1", "CD14", "S100A9", "NCF1", "S100A8", "KCTD12")

top10markers <- c("CD37", "IFITM3", "ARPP21", "FRMD4B", "IGFBP2", "LINC01013", "RFLNB", "IRX1", "H2AC14", "H3C2", "H3C4", "H2BC7", "H2AC11", "H2AC16", "H1-2", "H2AC20", "H2AC17", "H1-5", "CHST6", "MT1G", "TRIM14", "PLS3", "KCNJ12", "CLEC2D", "NLGN1", "GPR160", "BLNK", "CRIP1", "CFD", "AZU1", "PDE4D", "RETN", "ANXA1", "TFEC", "CSF3R", "WLS", "AIF1", "ELANE", "CD79A", "CCDC26", "SIPA1L2", "PPA1", "NRIP1", "CDK6", "VPREB1", "PSME2", "CCT3", "ZNF608", "CD14", "S100A9", "NCF1", "S100A8", "KCTD12", "JAML", "VNN2", "ITGB2", "VASP", "CYBB")


p <- DotPlot(alldataF, features = top10allclusters,
             cols = clustercols,
             dot.scale = 10,
             #split.by = "RNA_snn_res.0.5",
             group.by = "RNA_snn_res.0.5") +
  RotatedAxis()
p


df<- p$data
exp_mat<-df %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()

percent_mat<-df %>% 
  select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()


library(viridis)
library(Polychrome)
Polychrome::swatch(inferno(20, direction = -1))


quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))


percent_mat <- t(percent_mat)
exp_mat <- t(exp_mat)


## any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20, direction = -1)[c(1,10, 20)])


cell_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}


order_genes <- c("H2BC7", "H2AC14", "H1-5", "H3C2", "H2AC11", "H3C4", "H3C3", "H2AC16", "MAD2L1", "H2AC20", "MT1X", "BLNK", "AFF3", "LCN6", "PRKCZ", "OFD1", "VPREB3", "IRAG2", "STK39", "BLK", "SRM", "NME1", "IL7R", "IRX1", "UGT3A2", "RAG1", "CD79A", "DANCR", "ELK3", "SPATS2L", "TSPOAP1", "TNFAIP8", "PLEK", "CCDC26", "RGS18", "ATP8B4", "ANKRD28", "P4HB", "HERC5", "SIPA1L2", "ANXA1", "APLP2", "SRGN", "CFD", "MNDA", "SERPINB1", "LYZ", "AIF1", "S100A11", "HCST", "S100A8", "S100A9", "CYBB", "ITGB2", "NCF1", "LST1", "PIM1", "PTPRC", "FCER1G")


## To make the size of the dot in the heatmap body comparable to the legend, I used fixed
## size unit.(2, "mm) rather than min(unit.c(w, h).
cluster_anno<-  c("Early-inter", "Inicial B-cell", "Early", "Inter/late", "Inter-late", "Late Macrophage")

column_ha<- HeatmapAnnotation(
  cluster_anno = cluster_anno,
  col = list(cluster_anno = setNames(brewer.pal(6, "Paired"), unique(cluster_anno))
  ),
  na_col = "grey"
)


layer_fun1 = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= sqrt(pindex(percent_mat, i, j)/100)  * unit(2, "mm"),
              gp = gpar(fill = col_fun(pindex(exp_mat, i, j)), col = NA))}

lgd_list1 = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt",
          graphics = list(
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0  * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                             gp = gpar(fill = "black")))
  ))


hp<- Heatmap(exp_mat,
             heatmap_legend_param=list(title="expression"),
             column_title = "clustered dotplot", 
             col=col_fun,
             rect_gp = gpar(type = "none"),
             layer_fun = layer_fun,
             row_names_gp = gpar(fontsize = 5),
             split = cluster_anno,
             border = "black",
             width = unit(20, "cm") , height = unit(6, "cm"),
             column_order = order_genes,
             row_order = c("3", "4", "1", "6", "5", "7")
             #row_annotation = column_ha
)

draw( hp, annotation_legend_list = lgd_list)
