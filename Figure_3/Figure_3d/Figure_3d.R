
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
library(org.Hs.eg.db)

vst <- vst(out.DESeq2$dds)
counts.trans <- assay(varianceStabilizingTransformation(out.DESeq2$dds))


anno_df_biomart <- get_annotation(dds = out.DESeq2$dds,
                                  biomart_dataset = "hsapiens_gene_ensembl",
                                  idtype = "ensembl_gene_id_version")


write_csv(anno_df_biomart, "annotation_hg38.csv")

}
# subset the background to include only the expressed genes
bg_ids <- rownames(out.DESeq2$dds)[rowSums(counts(out.DESeq2$dds)) > 50]


pca2go_Hs <- pca2go(vst,
                    annotation = anno_df_biomart,
                    annopkg = "org.Hs.eg.db",
                    ensToGeneSymbol = TRUE,
                    background_genes = bg_ids)


# and finally, with all the objects prepared...
pcaExplorer(dds = out.DESeq2$dds, dst = vst, annotation = anno_df_biomart, pca2go = pca2go_Mm)



plot <- SARTools::PCAPlot(counts.trans=counts.trans, group=target[,varInt], col=colors, outfile = F, ggplot_theme = theme_bw())


n=min(500, nrow(counts.trans))
rv = apply(counts.trans, 1, var, na.rm=TRUE)
pca = prcomp(t(counts.trans[order(rv, decreasing = TRUE), ][1:n,]))
prp <- pca$sdev^2 * 100 / sum(pca$sdev^2)
prp <- round(prp[1:3],2)



pcaExplorer(dds = out.DESeq2$dds,
            dst = vst,
            annotation = anno_df_biomart)
 
colors <- c("#6491c0","#6491c0","#6491c0",
            "#bd8fbd","#bd8fbd","#bd8fbd",
            "#ed8b56", "#ed8b56", "#ed8b56")#,
#  "##2e68b1","##2e68b1","##2e68b1",
#   "##bd8fbd", "##bd8fbd", "##bd8fbd")
labels <- rownames(plot)
colorsg <- setNames(colors,labels )

sampleDists <- dist(t(vsd))


save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}




sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

p <- pheatmap(sampleDistMatrix,
              clustering_distance_rows=sampleDists,
              clustering_distance_cols=sampleDists,
              col=colors, 
              cellwidth=10, cellheight=10)

save_pheatmap_pdf(p, "HOXD4_distanceMatrix.pdf")


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
ggsave("HOXD4_RNAseq_PCA.pdf")
