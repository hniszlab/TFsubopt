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

CellsTypeCluster <- group_by(meta.data, type, RNA_snn_res.0.5) %>% summarise(count = n())
CellsTypeCluster <- CellsTypeCluster %>%
  group_by(type) %>%
  mutate(per =  count/sum(count)) %>% 
  ungroup
CellsTypeCluster

CellsTypeCluster$type <- factor(CellsTypeCluster$type, levels = c("WT", "AroPerfect_IS15", "AroPerfect_IS10"))
p <- ggplot(CellsTypeCluster, aes(x = type, stratum = RNA_snn_res.0.5, alluvium = RNA_snn_res.0.5, y = per, fill = RNA_snn_res.0.5, label = RNA_snn_res.0.5)) +
  geom_flow(alpha = .6) +
  geom_stratum(color = 'black') +
  theme_classic() +
  scale_fill_manual(values=clustercols)
p
