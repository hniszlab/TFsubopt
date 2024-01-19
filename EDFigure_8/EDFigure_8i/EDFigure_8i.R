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


#Diff expression 5, 6, 7 bw samples
# select all cells in cluster 1
alldata7 <- SetIdent(alldataF, value = "RNA_snn_res.0.5")
alldata7 <- subset(alldata7,  idents = 7)
alldata7 <- SetIdent(alldata7, value = "type")
alldata7 <- subset(alldata7,  idents = c("AroPerfect_IS15", "WT"))
cell_selection <- SetIdent(alldata7, value = "type")
# Compute differentiall expression
DGE_cell_selection <- FindAllMarkers(cell_selection, logfc.threshold = 0.2, test.use = "DESeq2", 
                                     min.pct = 0.1, min.diff.pct = 0.2, only.pos = F, max.cells.per.ident = 200, 
                                     assay = "RNA")


DGE_cell_selection = DGE_cell_selection %>% filter(cluster == "WT")
top25_cell_of7ontype <- DGE_cell_selection %>% group_by(cluster) %>% top_n(-25, p_val_adj)
write.csv(top25_cell_of7ontype, "top25_cell_of7ontypeWT_IS15_desq2.csv")
write.csv(DGE_cell_selection, "DGE_cell_selection.csv")



DGE_cell_selection["log10pval"] <- -1 * log10(DGE_cell_selection$p_val)
DGE_cell_selection["log10Adjpval"] <- -1 * log10(DGE_cell_selection$p_val_adj)


# Create new categorical column ------------------------------------------------ 
DGE_cell_selection <- DGE_cell_selection %>%
  mutate(gene_type = case_when(avg_log2FC >= 0.5 & p_val_adj <= 0.05 ~ "up",
                               avg_log2FC <= -0.5 & p_val_adj <= 0.05 ~ "down",
                               TRUE ~ "ns"))  


DGE_cell_selection %>%
  distinct(gene_type) %>%
  pull()  



# Add colour, size and alpha (transparency) to volcano plot --------------------
cols <- c("up" = "#26b3ff", "down" = "#26b3ff", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)



myfeatures <- c("FCGR2A", "FCGR2B", "FCGR1A","FCGR1B", "HSPA6", "CEACAM1", "CEACAM8")



#genes.to.labelT = DGE_cell_selection  %>% top_n(-25, p_val_adj)
#genes.to.labelT = genes.to.labelT %>% arrange(desc(avg_log2FC))
#genes.to.label = as.character(unique(genes.to.labelT$gene))




genes.to.label2 <- DGE_cell_selection %>% filter(gene %in% myfeatures)


graph <- ggplot(DGE_cell_selection, aes(
  x = avg_log2FC, 
  y = log10Adjpval,
  fill = gene_type,   
  size = gene_type,
  alpha = gene_type,
  #label = gene
)) +
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") +
  labs(x = 'Avg Log2 FC', y = '-Log10 p-value (adj)') +
  ggtitle('Wild type vs AroPerfect IS15') +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    aspect.ratio = 2/1.5) +
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) # Modify point transparency
# geom_text_repel(data = . %>%
#                   mutate(gene = ifelse(gene %in% genes.to.label, gene, "")),
#                   aes(label = gene),
#                   min.segment.length = 0,
#                   seed = 42,
#                   box.padding = 0.5,
#                 max.overlaps = Inf
#                   )


graph <- LabelPoints(plot = graph, points = myfeatures, repel = T, min.segment.length = 0, box.padding = 0.5, max.overlaps = 2000000000)