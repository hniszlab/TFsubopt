library(tidyverse)
#library('DESeq2')
library(ComplexHeatmap)
library(RColorBrewer)
library(dendextend)
library(viridis)


DEG_RFC = read.csv("DEG_list.csv" , header = T)
Allvst = read_csv("AllCount-vst.csv")


longAllvst <- gather(Allvst, "group", "Expression", -ID)

# create a new column that as the name of the samples without the _
# and the number of the replicate
longAllvst = longAllvst %>% separate(group, c("tgroup", "tgroup2", NA), sep = "_", remove = F)%>% unite("tgroup", tgroup:tgroup2, sep= "_")

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

col.order <- c("ZIP13K2_WT", "NGN2_WT", "NGN2_AroLITE","NGN2_AroPERFECT")

MeansDEGSvstSC <- as.data.frame(MeansDEGSvstSC)
MeansDEGSvstSC = MeansDEGSvstSC %>% dplyr::select(col.order)


MeansDEGSvstSC = MeansDEGSvstSC %>% as.matrix()

write.csv(MeansDEGSvstSC, file= "MeansDEGSvstSC.csv")


NeuronMarkers = read_csv("NeuronMarkers copy.csv")

MarkersExp <- MeansAllvstSC[which(row.names(MeansAllvstSC) %in% NeuronMarkers$GeneID),] %>% as.data.frame() %>% rownames_to_column("GeneID")

NeuronMarkers <- NeuronMarkers %>% left_join(MarkersExp, by =  "GeneID") %>% column_to_rownames("GeneName")


col.order <- c("GeneID","PMC8452516",  "ZIP13K2_WT", "NGN2_WT", "NGN2_AroLITE","NGN2_AroPERFECT")

NeuronMarkers <- as.data.frame(NeuronMarkers)
NeuronMarkers = NeuronMarkers %>% dplyr::select(all_of(col.order))

#row_dendo = hclust(dist(DEGvstMeansSC), method = "complete") # row clustering
kclus3 <- kmeans(NeuronMarkers[, 3:6], 4)
split3 <- kclus3$cluster

colours <- list("PMC8452516" = c("Progenitors" = "#c7006e", "Pluripotency" = "#d70000","CD99+ cells" = "#efdb00", "PRPH+/PHOX2B+ neurons-1" =  "#c7006e", "PRPH+/POU4F1+ neurons-1" = "#f47100", "PRPH+/POU4F1+ neurons-2" = "#ba00e5", "GPM6A+ neurons" = "#41c300", "PRPH+/PHOX2B+ neurons-2" = "#d70000"))

subMarker_df <- NeuronMarkers %>% mutate(ID = row_number())

subsetlist <- c("GPM6B", "GNG5", "SSR2", "H2AZ2", "GPM6A", "PCDH9", "TMEM97", "DLST", "MPPED2", "FGF13", "SHOX2", "BHLHE22", "EMX1", "TLE1", "TUBB1", "SERTM1", "TBX1", "POU5F1", "SOX2", "NANOG", "EOMES", "FOXG1", "NR2F1", "EN1", "NEUROD1", "HES1", "PAX6")

subMarker_df <- subMarker_df[which(row.names(subMarker_df) %in% subsetlist),]

#kclus3 <- kmeans(subMarker_df[, 3:6], 3)
#split3 <- kclus3$cluster



subset = subMarker_df$ID
labels = row.names(subMarker_df)

rowAnn <- HeatmapAnnotation(df = subMarker_df[, 2],
                            which = 'row',
                            col = colours,
                            annotation_width = unit(c(1, 1), 'cm'),
                            gap = unit(1, 'mm'))

ha = rowAnnotation(foo = anno_mark(at = subset, labels = labels, labels_gp = gpar(fontsize = 10)))

h2 <- Heatmap(subMarker_df[, 3:6],
              name = "z-score", #title of legend
              #col = colorRampPalette(rev(brewer.pal(n = 9, name ="Spectral")))(10),
              col = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(10),
              #cluster_rows = color_branches(row_dendo, k = 9),
              column_order = order(as.numeric(gsub("column", "", colnames(subMarker_df[, 3:6])))),
              cluster_rows = F, cluster_columns = F,
              show_row_names = T,
              row_names_gp = gpar(fontsize = 2),
              #split=split3,
              cluster_row_slices = FALSE,
              width = ncol(subMarker_df)*unit(5, "mm"), 
              height = nrow(subMarker_df)*unit(5, "mm"),
              right_annotation = rowAnn
)# + ha

h2 = draw(h2)

h2 <- Heatmap(subMarker_df[, 3:6],
              name = "z-score", #title of legend
              #col = colorRampPalette(rev(brewer.pal(n = 9, name ="Spectral")))(10),
              col = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(10),
              #cluster_rows = color_branches(row_dendo, k = 9),
              column_order = order(as.numeric(gsub("column", "", colnames(subMarker_df[, 3:6])))),
              cluster_rows = F, cluster_columns = F,
              show_row_names = F,
              row_names_gp = gpar(fontsize = 2),
              #split=split3,
              cluster_row_slices = FALSE,
              width = ncol(subMarker_df)*unit(5, "mm"), 
              height = nrow(subMarker_df)*unit(5, "mm"),
              right_annotation = rowAnn
) + ha
