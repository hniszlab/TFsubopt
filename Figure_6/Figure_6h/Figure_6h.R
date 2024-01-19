library('tidyverse')
#library('DESeq2')
library('EnhancedVolcano')
#library(ComplexHeatmap)
library(RColorBrewer)
library(dendextend)
library(viridis)
library(gridExtra)
library(grid)
library(Seurat)
library(cowplot)

KOvsWTGFP = read_csv("ZIP13K2xxWTvsNGN2xxWT.complete.csv") %>% distinct(Id, .keep_all= TRUE) %>%  column_to_rownames("Id")
AroPLUSGFPvsWTGFP = read_csv("NGN2xxAroLITEvsNGN2xxWT.complete.csv") %>% distinct(Id, .keep_all= TRUE) %>% column_to_rownames("Id")
AroPERFECTGFPvsWTGFP = read_csv("NGN2xxAroPERFECTvsNGN2xxWT.complete.csv") %>% distinct(Id, .keep_all= TRUE) %>% column_to_rownames("Id")
NGN2_targets = read_csv("ZIP13K2_WTvsNGN2_WT.csv")
NGN2_targets_Chip = read_csv("ChipNGN2_targets.csv")


Annotation_HG38 = read_csv("Annotation_HG38.csv")

NGN2_targets <- NGN2_targets %>% left_join(Annotation_HG38, by = c("Id" = "GeneID"))
NGN2_targets_Chip <- NGN2_targets_Chip %>% left_join(Annotation_HG38, by = c("ID" = "GeneID"))


KOvsWTGFP$log2FoldChange <- KOvsWTGFP$log2FoldChange * -1


#hgene = (NGN2_targets$Id)
hgene = c("TMEM97", "SERTM1")
keyvals.colour <- ifelse(
  row.names(KOvsWTGFP) %in% hgene, 'blue','gray')
keyvals.colour[is.na(keyvals.colour)] <- 'gray'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'NGN2 Target'
names(keyvals.colour)[keyvals.colour == 'gray'] <- 'Other'

p1 <- EnhancedVolcano(KOvsWTGFP,
                      lab = KOvsWTGFP$GeneSymbol,
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = 'KO versus WTGFP',
                      subtitle = NULL,
                      caption = NULL,
                      titleLabSize = 10,
                      pCutoff = 0.01,
                      FCcutoff = 1.5,
                      pointSize = 1.0,
                      labSize = 3,
                      arrowheads = FALSE,
                      selectLab = hgene,
                      col = c('gray50', 'gray50', 'gray50', 'gray50'),
                      colCustom = keyvals.colour,
                      colAlpha = 1,
                      axisLabSize = 10,
                      drawConnectors = T,
                      widthConnectors = 0.75,
                      legendPosition = 'top',
                      legendLabSize = 8,
                      legendIconSize = 2.0,
                      border = 'full',
                      boxedLabels = F)
p2 <- EnhancedVolcano(AroPERFECTGFPvsWTGFP,
                      lab = AroPERFECTGFPvsWTGFP$GeneSymbol,
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = 'AroPERFECT versus WTGFP',
                      subtitle = NULL,
                      caption = NULL,
                      titleLabSize = 10,
                      pCutoff = 0.01,
                      FCcutoff = 1.5,
                      pointSize = 1.0,
                      labSize = 3,
                      arrowheads = FALSE,
                      selectLab = hgene,
                      #col = c('gray50', 'gray50', 'gray50', 'gray50'),
                      colCustom = keyvals.colour,
                      colAlpha = 1,
                      axisLabSize = 10,
                      drawConnectors = T,
                      widthConnectors = 0.75,
                      legendPosition = 'top',
                      legendLabSize = 8,
                      legendIconSize = 2.0,
                      border = 'full',
                      boxedLabels = F)
p3 <- EnhancedVolcano(AroPLUSGFPvsWTGFP,
                      lab = AroPLUSGFPvsWTGFP$GeneSymbol,
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = 'AroPLUS versus WTGFP',
                      subtitle = NULL,
                      caption = NULL,
                      titleLabSize = 10,
                      pCutoff = 0.01,
                      FCcutoff = 1.5,
                      pointSize = 1.0,
                      labSize = 3,
                      arrowheads = FALSE,
                      selectLab = hgene,
                      #col = c('gray50', 'gray50', 'gray50', 'gray50'),
                      colCustom = keyvals.colour,
                      colAlpha = 1,
                      axisLabSize = 10,
                      drawConnectors = T,
                      widthConnectors = 0.75,
                      legendPosition = 'top',
                      legendLabSize = 8,
                      legendIconSize = 2.0,
                      border = 'full',
                      boxedLabels = F)

pdf(file = "volcano_NGN2_CHIPseqTargets.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 4)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()