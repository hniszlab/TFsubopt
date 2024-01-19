library(biomaRt)
library(SARTools)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(pcaExplorer)
library(Seurat)
library(cowplot)
library(topGO)
library(org.Mm.eg.db)

vst <- vst(out.DESeq2$dds)
counts.trans <- assay(varianceStabilizingTransformation(out.DESeq2$dds))

sampleDists <- dist(t(vsd))

n=min(500, nrow(counts.trans))
rv = apply(counts.trans, 1, var, na.rm=TRUE)
pca = prcomp(t(counts.trans[order(rv, decreasing = TRUE), ][1:n,]))
prp <- pca$sdev^2 * 100 / sum(pca$sdev^2)
prp <- round(prp[1:3],2)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

p <- pheatmap(sampleDistMatrix,
              clustering_distance_rows=sampleDists,
              clustering_distance_cols=sampleDists,
              col=colors, 
              cellwidth=10, cellheight=10)

pcaData = as.data.frame(plot)
p <- ggplot(pcaData, aes(PC1, PC2, color=colorsg)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",prp[1],"% variance")) +
  ylab(paste0("PC2: ",prp[2],"% variance")) + 
  coord_fixed() +
  theme_bw() + 
  theme(aspect.ratio = 1, legend.position = "none", panel.grid.major = element_blank() , panel.grid.minor = element_blank())
pg<- LabelPoints(plot = p, points = labels, repel = T, min.segment.length = 0, box.padding = 0.5, max.overlaps = Inf)

p2 <- ggplot(pcaData, aes(PC1, PC3, color=colorsg)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",prp[1],"% variance")) +
  ylab(paste0("PC3: ",prp[3],"% variance")) + 
  coord_fixed() +
  theme_bw() + 
  theme(aspect.ratio = 1, legend.position = "none", panel.grid.major = element_blank() , panel.grid.minor = element_blank())
pg2<- LabelPoints(plot = p2, points = labels, repel = T, min.segment.length = 0, box.padding = 0.5, max.overlaps = Inf)

plot_grid(ncol = 2, pg, pg2)