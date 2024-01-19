library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)
library(tidyverse)
library("readxl")


my_data_metapredict <- read_csv("MasterTable_qIDR_17042022.csv")
#my_data_metapredict=my_data_metapredict[!(my_data_metapredict$gene_symbol %in% ZF_keys), ]
#my_data_metapredict_top=my_data_metapredict[my_data_metapredict$min_pval < 0.01 & !(my_data_metapredict$gene_symbol %in% ZF_keys), ]

my_data_metapredict_top=my_data_metapredict[my_data_metapredict$min_pval < 0.01, ]

#x=data.frame(my_data_metapredict_top$gene_symbol,my_data_metapredict_top$min_pval )
x=data.frame(my_data_metapredict$gene_symbol,my_data_metapredict$min_pval )

colnames(x)=c("Gene", "metric")
newdata <- x[order(x$metric),] 
newdata$metric=-log10(newdata$metric)

## Extract the foldchanges
foldchanges <- newdata$metric

## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- newdata$Gene

foldchanges <- sort(foldchanges, decreasing = TRUE)



# GSEA using gene sets associated with BP Gene Ontology terms
gseaGO <- gseGO(geneList = foldchanges, 
                OrgDb = org.Hs.eg.db, 
                ont = 'BP', 
                minGSSize = 20, 
                pvalueCutoff = 0.05,
                verbose = FALSE,keyType = 'SYMBOL') 


gseaGO_results <- gseaGO@result



dotplot(gseaGO, showCategory=30) + ggtitle("dotplot for GSEA")



# GSEA analysis for Rna Binding Proteins 
rbp <- read.gmt("gmt/rna_binding_all.gmt.txt")
pld <- read.gmt("gmt/PLDs.gmt.txt")
tfs <- read.gmt("gmt/tfs.gmt.txt")
tmp <- read.gmt("gmt/TMPs.gmt.txt")
trnsp <- read.gmt("gmt/transporters.gmt.txt")

tfsz = tfs[!(tfs$gene %in% ZF_keys),]



msigrbp <- GSEA(foldchanges, TERM2GENE=rbp, verbose=FALSE)
msigpld <- GSEA(foldchanges, TERM2GENE=pld, verbose=FALSE)
msigtfs <- GSEA(foldchanges, TERM2GENE=tfs, verbose=FALSE)
msigtmp <- GSEA(foldchanges, TERM2GENE=tmp, verbose=FALSE)
gseaGO_results_rbp <- msigrbp@result
gseaGO_results_pld <- msigpld@result
gseaGO_results_tfs <- msigtfs@result
gseaGO_results_tmp <- msigtmp@result


#gseaplot2(msig, geneSetID = rna binding proteins, title = msig$Description[1], base_size = 6)



gseaplot2(msigrbp, geneSetID = "rna binding proteins", title = "RNA binding proteins", base_size = 6, subplots = 1:2)
gseaplot2(msigpld, geneSetID = "PLDs", title = "PLDs", base_size = 6, subplots = 1:2)
gseaplot2(msigtfs, geneSetID = "tfs", title = "Transcription factors", base_size = 6, subplots = 1:2)
gseaplot2(msigtmp, geneSetID = "TMP", title = "Transporters", base_size = 6, subplots = 1:2)




my_data_metapredict_top=my_data_metapredict[my_data_metapredict$min_pval < 0.01, ]

y=data.frame(my_data_metapredict_top$gene_symbol,my_data_metapredict_top$min_pval )
#x=data.frame(my_data_metapredict$gene_symbol,my_data_metapredict$min_pval )

colnames(y)=c("Gene", "metric")
newdata001 <- y[order(y$metric),] 
newdata001$metric=-log10(newdata001$metric)

## Extract the foldchanges
foldchanges001 <- newdata001$metric

## Name each fold change with the corresponding Entrez ID
names(foldchanges001) <- newdata001$Gene

foldchanges001 <- sort(foldchanges001, decreasing = TRUE)




msigrbp001 <- GSEA(foldchanges001, TERM2GENE=rbp, verbose=FALSE)
msigpld001 <- GSEA(foldchanges001, TERM2GENE=pld, verbose=FALSE)
msigtfs001 <- GSEA(foldchanges001, TERM2GENE=tfs, verbose=FALSE)
msigtrnsp001 <- GSEA(foldchanges001, TERM2GENE=trnsp, verbose=FALSE)
gseaGO_results_rbp001 <- msigrbp001@result
gseaGO_results_pld001 <- msigpld001@result
gseaGO_results_tfs001 <- msigtfs001@result
gseaGO_results_trnsp001 <- msigtrnsp001@result



gseaplot2(msigrbp001, geneSetID = "rna binding proteins", title = "RNA binding proteins", base_size = 6, subplots = 1:2)
gseaplot2(msigpld001, geneSetID = "PLDs", title = "PLDs", base_size = 6, subplots = 1:2)
gseaplot2(msigtfs001, geneSetID = "tfs", title = "Transcription factors", base_size = 6, subplots = 1:2)
gseaplot2(msigtrnsp001, geneSetID = "Transporters", title = "Transporters", base_size = 6, subplots = 1:2)



my_data_metapredict_only <- my_data_metapredict %>% 
  filter(my_data_metapredict$overlap_metapredict  =="yes" )


z=data.frame(my_data_metapredict_only$gene_symbol,my_data_metapredict_only$min_pval )
#x=data.frame(my_data_metapredict$gene_symbol,my_data_metapredict$min_pval )

colnames(y)=c("Gene", "metric")
newdataMETA <- y[order(y$metric),] 
newdataMETA$metric=-log10(newdataMETA$metric)

## Extract the foldchanges
foldchangesMETA <- newdataMETA$metric

## Name each fold change with the corresponding Entrez ID
names(foldchangesMETA) <- newdataMETA$Gene

foldchangesMETA <- sort(foldchangesMETA, decreasing = TRUE)



msigrbpMETA <- GSEA(foldchangesMETA, TERM2GENE=rbp, verbose=FALSE)
msigpldMETA <- GSEA(foldchangesMETA, TERM2GENE=pld, verbose=FALSE)
msigtfsMETA <- GSEA(foldchangesMETA, TERM2GENE=tfs, verbose=FALSE)
msigtrnspMETA <- GSEA(foldchangesMETA, TERM2GENE=trnsp, verbose=FALSE)
gseaGO_results_rbpMETA <- msigrbpMETA@result
gseaGO_results_pldMETA <- msigpldMETA@result
gseaGO_results_tfsMETA <- msigtfsMETA@result
gseaGO_results_trnspMETA <- msigtrnspMETA@result


gseaplot2(msigrbpMETA, geneSetID = "rna binding proteins", title = "RNA binding proteins", base_size = 6, subplots = 1:2)
gseaplot2(msigpldMETA, geneSetID = "PLDs", title = "PLDs", base_size = 6, subplots = 1:2)
gseaplot2(msigtfsMETA, geneSetID = "tfs", title = "Transcription factors", base_size = 6, subplots = 1:2)
gseaplot2(msigtrnspMETA, geneSetID = "Transporters", title = "Transporters", base_size = 6, subplots = 1:2)




gsea_c <- read.gmt("./gmt/gsea_Custom.gmt.txt")
gsea_cz = gsea_c[!(gsea_c$gene %in% ZF_keys),]
msigAll <- GSEA(foldchanges, TERM2GENE=gsea_c, verbose=FALSE)
gseaGO_results_All <- msigAll@result

msigAll001 <- GSEA(foldchanges001, TERM2GENE=gsea_c, verbose=FALSE)
msigAll001z <- GSEA(foldchanges001, TERM2GENE=gsea_cz, verbose=FALSE)
gseaGO_results_All001 <- msigAll001@result


g <- gseaplot2(msigAll, geneSetID = c("PLDs", "TFs") , title = "", base_size = 6, subplots = 1:2, rel_heights = c(1.5, 0.5, 0.5), color = c("green","blue"))
g1 <- gseaplot2(msigAll, geneSetID = c("PLDs", "RBP") , title = "", base_size = 6, subplots = 1:2, rel_heights = c(1.5, 0.5, 0.5), color = c("green","blue"))
g2 <- gseaplot2(msigAll, geneSetID = c("PLDs", "PLD99") , title = "", base_size = 6, subplots = 1:2, rel_heights = c(1.5, 0.5, 0.5), color = c("green","blue"))

c("PLDs", "TFs", "RBP", "Transporters")


plot_grid(g, g1, g2, ncol = 1, nrow = 3 )
ggsave("gsea_qIDR_Final_v2.pdf")



g <- gseaplot2(msigAll001, geneSetID = c("PLD99", "TFs") , title = "", base_size = 6, subplots = 1:2, rel_heights = c(1.5, 0.5, 0.5), color = c("green","blue"))
g1 <- gseaplot2(msigAll001z, geneSetID = c("PLD99", "TFs") , title = "", base_size = 6, subplots = 1:2, rel_heights = c(1.5, 0.5, 0.5), color = c("green","blue"))
g2 <- gseaplot2(msigAll001, geneSetID = c("PLDs", "PLD99") , title = "", base_size = 6, subplots = 1:2, rel_heights = c(1.5, 0.5, 0.5), color = c("green","blue"))

c("PLDs", "TFs", "RBP", "Transporters")


plot_grid(g, g1, g2, ncol = 1, nrow = 3 )
ggsave("gsea_qIDR_Final_v20_001.pdf")
