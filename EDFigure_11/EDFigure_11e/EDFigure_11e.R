
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

MYODtargets  = read_csv("Targets.csv")

C2C12xxWTvsMYOD1xxWT = read_csv("C2C12xxWTvsMYOD1xxWT.complete.csv") %>% distinct(Id, .keep_all= TRUE) %>%  column_to_rownames("Id")
MYOD1xxAroPerfectvsMYOD1xxWT = read_csv("MYOD1xxAroPerfectvsMYOD1xxWT.complete.csv") %>% distinct(Id, .keep_all= TRUE) %>%  column_to_rownames("Id")
MYOD1xxAroLitevsMYOD1xxWT = read_csv("MYOD1xxAroLitevsMYOD1xxWT.complete.csv") %>% distinct(Id, .keep_all= TRUE) %>%  column_to_rownames("Id")
MYOD1xxAroLiteCvsMYOD1xxWT = read_csv("MYOD1xxAroLiteCvsMYOD1xxWT.complete.csv") %>% distinct(Id, .keep_all= TRUE) %>%  column_to_rownames("Id")
MYOD1xxAroPerfectCvsMYOD1xxWT = read_csv("MYOD1xxAroPerfectCvsMYOD1xxWT.complete.csv") %>% distinct(Id, .keep_all= TRUE) %>%  column_to_rownames("Id")



Annotation_mm10 = read_csv("annotation_mm10.csv")

MYODtargets <- MYODtargets %>% left_join(Annotation_mm10, by = c("Id" = "gene_id"))


C2C12xxWTvsMYOD1xxWT$log2FoldChange <- C2C12xxWTvsMYOD1xxWT$log2FoldChange * -1


Listgenes <- c("Atp1b2", "Mcam", "Parm1", "Stap2", "Ccdc136")

#  "Megf10", "Rap1gap", "Epha2", "Ubash3b", , "Mcam", "Jag2", , "Dgkg", "Rdh5", "Slc35g1", "Lman1l", "Htr2b", "Slc31a2", "Edem2", "Tmx4", #"Mcam", "Zfp185", "Col4a5", "Sytl1", "Slc25a19", "Pam", "Ripor3", "Ltbp3")

#hgene = (MYODtargets$Id)
hgene = c(Listgenes)
keyvals.colour <- ifelse(
  row.names(C2C12xxWTvsMYOD1xxWT) %in% hgene, 'blue','gray')
keyvals.colour[is.na(keyvals.colour)] <- 'gray'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'MYOD1 Target'
names(keyvals.colour)[keyvals.colour == 'gray'] <- 'Other'


p1 <- EnhancedVolcano(C2C12xxWTvsMYOD1xxWT,
                      lab = C2C12xxWTvsMYOD1xxWT$Gene_name,
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = 'MYOD1 WT vs C2C12 wild type',
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
                      boxedLabels = F)+ 
  ggplot2::coord_cartesian(xlim=c(-15, 15)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-15, 15, 5))
p2 <- EnhancedVolcano(MYOD1xxAroLitevsMYOD1xxWT,
                      lab = MYOD1xxAroLitevsMYOD1xxWT$Gene_name,
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = 'MYOD1 AroLite vs MYOD1 WT',
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
                      boxedLabels = F)+ 
  ggplot2::coord_cartesian(xlim=c(-15, 15)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-15, 15, 5))
p3 <- EnhancedVolcano(MYOD1xxAroPerfectvsMYOD1xxWT,
                      lab = MYOD1xxAroPerfectvsMYOD1xxWT$Gene_name,
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = 'MYOD1 AroPerfect vs MYOD1 WT',
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
                      boxedLabels = F)+ 
  ggplot2::coord_cartesian(xlim=c(-15, 15)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-15, 15, 5))
p4 <- EnhancedVolcano(MYOD1xxAroLiteCvsMYOD1xxWT,
                      lab = MYOD1xxAroLiteCvsMYOD1xxWT$Gene_name,
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = 'MYOD1 AroLite-C vs MYOD1 WT',
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
                      boxedLabels = F)+ 
  ggplot2::coord_cartesian(xlim=c(-15, 15)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-15, 15, 5))
p5 <- EnhancedVolcano(MYOD1xxAroPerfectCvsMYOD1xxWT,
                      lab = MYOD1xxAroPerfectCvsMYOD1xxWT$Gene_name,
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = 'MYOD1 AroPerfect-C vs MYOD1 WT',
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
                      boxedLabels = F)+ 
  ggplot2::coord_cartesian(xlim=c(-15, 15)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-15, 15, 5))


pdf(file = "volcano_MYOd1_RNAseqTargets_WGenes_w2.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 4)
grid.arrange(p1, p2, p3, p4, p5,  ncol=5) & ggplot2::coord_cartesian(xlim=c(-15, 15))
dev.off()

```