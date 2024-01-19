{r message=FALSE, warning=FALSE, include=FALSE}
setwd("/Users/magalhae/Desktop/Julian_RNAseqHOXD4/heatmaps/")
library(tidyverse)
#library('DESeq2')
library(ComplexHeatmap)
library(RColorBrewer)
library(dendextend)
library(viridis)

#Load data processed
load("HOXD4_heatmap.RData")



#DEG_RFC = read.csv("DEGS.txt" , header = T)
#Allvst = read_csv("AllCount-vst.csv")

longAllvst <- gather(Allvst, "group", "Expression", -ID)

# create a new column that as the name of the samples without the _
# and the number of the replicate
longAllvst = longAllvst %>% separate(group, c("tgroup", "tgroup2", NA), sep = "_", remove = F) %>% unite("tgroup", tgroup:tgroup2, sep= "_")


# this is a pipe. group_by will divide the table into groups,
# in this case by ID and replicate,
# summarize will create stats for the groups and we are selecting the mean
# and spread returns the table to the wide format.
MeansAllvst  <- longAllvst %>%
  group_by(ID, tgroup) %>%
  summarize(expression_mean = mean(Expression)) %>%
  spread(., tgroup, expression_mean)
# makes the column named ID as names of the rows

write.csv(MeansAllvst, file= "MeansAllvst.csv")


MeansAllvst <- MeansAllvst %>% column_to_rownames("ID")
MeansAllvstSC <- t(scale(t(MeansAllvst), scale = T, center = T))
write.csv(MeansAllvstSC, file= "MeansAllvstSC.csv")
hist(MeansAllvstSC)


MeansDEGSvstSC <- MeansAllvstSC[which(row.names(MeansAllvstSC) %in% DEG_RFC$Id),]
hist(MeansDEGSvstSC)
write.csv(MeansDEGSvstSC, file= "MeansDEGSvstSC.csv")


col.order <- c("Parental_WT", "HOXD4_KO", "HOXD4_WT", "HOXD4_AroPERFECT", "HOXD4_AroPLUS")

MeansDEGSvstSC <- as.data.frame(MeansDEGSvstSC)
MeansDEGSvstSC = MeansDEGSvstSC %>% dplyr::select(all_of(col.order))


MeansDEGSvstSC = MeansDEGSvstSC %>% as.matrix()
MeansDEGSvstSC <- na.omit(MeansDEGSvstSC) 
write.csv(MeansDEGSvstSC, file= "MeansDEGSvstSC.csv")



#row_dendo = hclust(dist(DEGvstMeansSC), method = "complete") # row clustering
kclus <- kmeans(MeansDEGSvstSC, 6)
split <- kclus$cluster


split2 <- factor(paste0("Cluster\n", kclus$cluster), levels=c("Cluster\n1","Cluster\n2","Cluster\n7","Cluster\n6", "Cluster\n4","Cluster\n5","Cluster\n8","Cluster\n3"))

h <- Heatmap(MeansDEGSvstSC,
             name = "z-score", #title of legend
             #col = colorRampPalette(rev(brewer.pal(n = 9, name ="Spectral")))(10),
             col = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(10),
             #cluster_rows = color_branches(row_dendo, k = 9),
             column_order = order(as.numeric(gsub("column", "", colnames(MeansDEGSvstSC)))),
             show_row_names = FALSE,
             split=split,
             cluster_row_slices = FALSE,
             width = unit(6, "cm") , height = unit(12, "cm")
)
pdf("HOXD4_RNAseq_Heatmap_test.pdf")
h = draw(h)
#dev.off()
