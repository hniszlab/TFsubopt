
library(tidyverse)
#library('DESeq2')
library(ComplexHeatmap)
library(RColorBrewer)
library(dendextend)
library(viridis)
library(cluster)
library(factoextra)
library(NbClust)



#open the DEG count Table (all samples in one table)

AllBulkscaled = read.csv("MARKER_setBulk.csv" , header = T, stringsAsFactors=T) %>% column_to_rownames('GeneID')
HExp_sc = read.csv("CellMarkers_higestExpSC.csv" , header = T)
HExp_elife = read.csv("HigestExp_elife.csv" , header = T)
ClusterTypes = read.csv("ClusterTypes.csv", header = T)
Gene_cluster = read.csv("Gene_cluster.csv", header = T)



MarkerBULK <- AllBulkscaled[which(row.names(AllBulkscaled) %in% Gene_cluster$GeneID),]
MarkerBULK <- na.omit(MarkerBULK) 
write.csv(MarkerBULK, file= "MarkerBULK.csv")



AllBulkscaledm = AllBulkscaled[3:14] %>% as.matrix()




row_dendo = hclust(dist(AllBulkscaledm), method = "complete") # row clustering
kclus <- kmeans(AllBulkscaledm,5)
split <- kclus$cluster
split2 <- Gene_cluster$Cluster



cc <- data.frame(AllBulkscaled$Cluster)
colnames(cc)<- c("Cluster")
colorsP <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#f781bf")
nr_cluster=sort(unique(cc$Cluster))
names(colorsP) = nr_cluster
colz = list(Cluster = colorsP)
ha = HeatmapAnnotation(df = cc, which = "row", width = unit(1, "cm"), col = colz)




fh = function(x) fastcluster::hclust(dist(x))

heat <- Heatmap(AllBulkscaledm,
                column_title = "All Genes from Elife markerset",
                name = "z-score", #title of legend
                #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                col = viridis(5),
                #cluster_rows = color_branches(row_dendo, k = 7),
                #clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
                #column_order = order(as.numeric(gsub("column", "", colnames(HExp_el_Bulkm)))),
                show_row_names = FALSE,
                split=split,
                cluster_row_slices = T,
                width = unit(6, "cm") , height = unit(12, "cm"),
                cluster_columns = fh,
                right_annotation=ha,
                use_raster = T
)
heat2 <- Heatmap(AllBulkscaledm,
                 name = "z-score", #title of legend
                 #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                 col = viridis(5),
                 #cluster_rows = color_branches(row_dendo, k = 7),
                 #clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
                 #column_order = order(as.numeric(gsub("column", "", colnames(HExp_el_Bulkm)))),
                 show_row_names = FALSE,
                 split=split2,
                 cluster_row_slices = T,
                 width = unit(6, "cm") , height = unit(12, "cm"),
                 cluster_columns = fh,
                 right_annotation=ha,
                 use_raster = T
)
h = draw(heat)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(heat, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(heat2, newpage = FALSE)
upViewport()

upViewport()



r.dend <- row_dend(h)  #If needed, extract row dendrogram
rcl.list <- row_order(h) #Extract clusters (output is a list)


lapply(rcl.list, function(x) length(x))


clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(GeneID = rownames(AllBulkscaledm[rcl.list[[i]],]),
                    ClusterHM = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)
Gene_cluster = Gene_cluster[c("GeneID", "Cluster")]
clu_df <- clu_df %>% left_join(Gene_cluster, by = "GeneID")
write.csv(clu_df, "AllMarkers_Elife_gene_row.csv", row.names = T)



clu_dfC <- clu_df %>% count(ClusterHM, Cluster, sort = TRUE)

{r fig.height=3, fig.width=2}
graph <- ggplot(clu_dfC, aes(x = Cluster, y = n, fill = Cluster)) +
  geom_bar(stat = "identity") +
  facet_grid( ClusterHM ~ . ) +
  scale_fill_manual(values=colorsP)+
  theme_bw() 
graph







HExp_el_Bulk <- AllBulkscaled[which(row.names(AllBulkscaled) %in% HExp_elife$GeneID),]
HExp_el_Bulk <- na.omit(HExp_el_Bulk) 
write.csv(HExp_el_Bulk, file= "HExp_el_Bulk.csv")
HExp_el_Bulkm = HExp_el_Bulk[3:14] %>% as.matrix()



row_dendo = hclust(dist(HExp_el_Bulkm), method = "complete") # row clustering
kclus <- kmeans(HExp_el_Bulkm,7)
split2 <- kclus$cluster
split <- HExp_el_Bulk$Cluster



cc <- data.frame(HExp_el_Bulk$Cluster)
colnames(cc)<- c("Cluster")
colorsP <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#f781bf")
nr_cluster=sort(unique(cc$Cluster))
names(colorsP) = nr_cluster
colz = list(Cluster = colorsP)
ha = HeatmapAnnotation(df = cc, which = "row", width = unit(1, "cm"), col = colz)



fh = function(x) fastcluster::hclust(dist(x))

heat <- Heatmap(HExp_el_Bulkm,
                column_title = "Higest Exp Genes from Elife markerset",
                name = "z-score", #title of legend
                #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                col = viridis(5),
                #cluster_rows = color_branches(row_dendo, k = 7),
                #clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
                #column_order = order(as.numeric(gsub("column", "", colnames(HExp_el_Bulkm)))),
                show_row_names = FALSE,
                split=split2,
                cluster_row_slices = T,
                width = unit(6, "cm") , height = unit(12, "cm"),
                cluster_columns = fh,
                right_annotation=ha,
                use_raster = T
)
heat2 <- Heatmap(HExp_el_Bulkm,
                 name = "z-score", #title of legend
                 #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                 col = viridis(5),
                 #cluster_rows = color_branches(row_dendo, k = 7),
                 #clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
                 #column_order = order(as.numeric(gsub("column", "", colnames(HExp_el_Bulkm)))),
                 show_row_names = FALSE,
                 split=split,
                 cluster_row_slices = T,
                 width = unit(6, "cm") , height = unit(12, "cm"),
                 cluster_columns = fh,
                 right_annotation=ha,
                 use_raster = T
)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(heat, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(heat2, newpage = FALSE)
upViewport()

upViewport()




HExp_SC_Bulk <- AllBulkscaled[which(row.names(AllBulkscaled) %in% HExp_sc$GeneID),]
HExp_SC_Bulk <- na.omit(HExp_SC_Bulk) 
write.csv(HExp_SC_Bulk, file= "HExp_SC_Bulk.csv")
HExp_SC_Bulkm = HExp_SC_Bulk[3:14] %>% as.matrix()



selected <- c( "kl", "ch", "hartigan",  "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "frey", "mcclain", "dunn", "hubert", "sdindex", "dindex", "sdbw")
results <- vector("list",19)
for (i in 1:length(selected)) {
  
  results[[i]] <- try(NbClust(HExp_SC_Bulkm, min.nc=2, max.nc=15, method="kmeans", index=selected[i]))
  
}



hist(results$Best.nc[1,], breaks = max(na.omit(results$Best.nc[1,])))




row_dendo = hclust(dist(HExp_SC_Bulkm), method = "complete") # row clustering
kclus <- kmeans(HExp_SC_Bulkm,14)
#split2 <- factor(kclus$cluster, levels = c())
split2 <- kclus$cluster
split <- HExp_SC_Bulk$Cluster
cc <- data.frame(HExp_SC_Bulk$Cluster)
colnames(cc)<- c("Cluster")
colorsP <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#f781bf")
nr_cluster=sort(unique(cc$Cluster))
names(colorsP) = nr_cluster
colz = list(Cluster = colorsP)
ha = HeatmapAnnotation(df = cc, which = "row", width = unit(1, "cm"), col = colz)



fh = function(x) fastcluster::hclust(dist(x))

heat <- Heatmap(HExp_SC_Bulkm,
                column_title = "Higest Exp Genes from SC based in Elife markerset",
                name = "z-score", #title of legend
                #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                col = viridis(5),
                #cluster_rows = color_branches(row_dendo, k = 7),
                #clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
                #column_order = order(as.numeric(gsub("column", "", colnames(HExp_el_Bulkm)))),
                show_row_names = FALSE,
                split=split2,
                cluster_row_slices = T,
                width = unit(6, "cm") , height = unit(12, "cm"),
                cluster_columns = fh,
                right_annotation=ha,
                use_raster = T
)
heat2 <- Heatmap(HExp_SC_Bulkm,
                 name = "z-score", #title of legend
                 #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                 col = viridis(5),
                 #cluster_rows = color_branches(row_dendo, k = 7),
                 #clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
                 #column_order = order(as.numeric(gsub("column", "", colnames(HExp_el_Bulkm)))),
                 show_row_names = FALSE,
                 split=split,
                 cluster_row_slices = T,
                 width = unit(6, "cm") , height = unit(12, "cm"),
                 cluster_columns = fh,
                 right_annotation=ha,
                 use_raster = T
)
h <- draw(heat)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(heat, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(heat2, newpage = FALSE)
upViewport()

upViewport()




r.dend <- row_dend(h)  #If needed, extract row dendrogram
rcl.list <- row_order(h) #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))
clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(GeneID = rownames(AllBulkscaledm[rcl.list[[i]],]),
                    ClusterHM = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)
Gene_cluster = Gene_cluster[c("GeneID", "Cluster")]
clu_df <- clu_df %>% left_join(Gene_cluster, by = "GeneID")
write.csv(clu_df, "HigExp_SC_gene_row.csv", row.names = T)



clu_dfC <- clu_df %>% count(ClusterHM, Cluster, sort = TRUE)

{r fig.height=3, fig.width=2}
graph <- ggplot(clu_dfC, aes(x = Cluster, y = n, fill = Cluster)) +
  geom_bar(stat = "identity") +
  facet_grid( ClusterHM ~ . ) +
  scale_fill_manual(values=colorsP)+
  theme_bw() 
graph







-------------
  TOP 50
------------
  
  
TOP50_higestExp_elife <- read.csv("TOP50_higestExp_elife.csv", header = T)
TOP50_el_Bulk <- AllBulkscaled[which(row.names(AllBulkscaled) %in% TOP50_higestExp_elife$GeneID),]
TOP50_el_Bulk <- na.omit(TOP50_el_Bulk) 
write.csv(TOP50_el_Bulk, file= "TOP50_el_Bulk.csv
          ")
TOP50_el_Bulkm = TOP50_el_Bulk[3:14] %>% as.matrix()
#HExp_el_Bulkm = HExp_el_Bulk[c(5, 6, 7, 8, 10, 11, 13, 14) ] %>% as.matrix()



row_dendo = hclust(dist(TOP50_el_Bulkm), method = "complete") # row clustering
kclus <- kmeans(TOP50_el_Bulkm,7)
split2 <- kclus$cluster
split <- TOP50_el_Bulk$Cluster



cc <- data.frame(TOP50_el_Bulk$Cluster)
colnames(cc)<- c("Cluster")
colorsP <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#f781bf")
nr_cluster=sort(unique(cc$Cluster))
names(colorsP) = nr_cluster
colz = list(Cluster = colorsP)
ha = HeatmapAnnotation(df = cc, which = "row", width = unit(1, "cm"), col = colz)



fh = function(x) fastcluster::hclust(dist(x))

heat <- Heatmap(TOP50_el_Bulkm,
                column_title = "Top50 Exp Genes from Elife markerset",
                name = "z-score", #title of legend
                #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                col = viridis(5),
                #cluster_rows = color_branches(row_dendo, k = 7),
                #clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
                #column_order = order(as.numeric(gsub("column", "", colnames(HExp_el_Bulkm)))),
                show_row_names = FALSE,
                split=split2,
                cluster_row_slices = T,
                width = unit(6, "cm") , height = unit(12, "cm"),
                cluster_columns = fh,
                right_annotation=ha,
                use_raster = T
)
heat2 <- Heatmap(TOP50_el_Bulkm,
                 name = "z-score", #title of legend
                 #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                 col = viridis(5),
                 #cluster_rows = color_branches(row_dendo, k = 7),
                 #clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
                 #column_order = order(as.numeric(gsub("column", "", colnames(HExp_el_Bulkm)))),
                 show_row_names = FALSE,
                 split=split,
                 cluster_row_slices = T,
                 width = unit(6, "cm") , height = unit(12, "cm"),
                 cluster_columns = fh,
                 right_annotation=ha,
                 use_raster = T
)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(heat, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(heat2, newpage = FALSE)
upViewport()

upViewport()



TOP50_higestExp_SC <- read.csv("TOP50_higestExp_SC.csv", header = T)
TOP50_higestExp_SC <- AllBulkscaled[which(row.names(AllBulkscaled) %in% TOP50_higestExp_SC$GeneID),]
TOP50_higestExp_SC <- na.omit(TOP50_higestExp_SC) 
write.csv(TOP50_higestExp_SC, file= "TOP50_higestExp_SC.csv
          ")
TOP50_higestExp_SCm = TOP50_higestExp_SC[3:14] %>% as.matrix()
#HExp_el_Bulkm = HExp_el_Bulk[c(5, 6, 7, 8, 10, 11, 13, 14) ] %>% as.matrix()



row_dendo = hclust(dist(TOP50_higestExp_SCm), method = "complete") # row clustering
kclus <- kmeans(TOP50_higestExp_SCm,7)
split2 <- kclus$cluster
split <- TOP50_higestExp_SC$Cluster



cc <- data.frame(TOP50_higestExp_SC$Cluster)
colnames(cc)<- c("Cluster")
colorsP <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#f781bf")
nr_cluster=sort(unique(cc$Cluster))
names(colorsP) = nr_cluster
colz = list(Cluster = colorsP)
ha = HeatmapAnnotation(df = cc, which = "row", width = unit(1, "cm"), col = colz)


fh = function(x) fastcluster::hclust(dist(x))

heat <- Heatmap(TOP50_higestExp_SCm,
                column_title = "Top50 Exp Genes from SC based on Elife markerset",
                name = "z-score", #title of legend
                #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                col = viridis(5),
                #cluster_rows = color_branches(row_dendo, k = 7),
                #clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
                #column_order = order(as.numeric(gsub("column", "", colnames(HExp_el_Bulkm)))),
                show_row_names = FALSE,
                split=split2,
                cluster_row_slices = T,
                width = unit(6, "cm") , height = unit(12, "cm"),
                cluster_columns = fh,
                right_annotation=ha,
                use_raster = T
)
heat2 <- Heatmap(TOP50_higestExp_SCm,
                 name = "z-score", #title of legend
                 #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                 col = viridis(5),
                 #cluster_rows = color_branches(row_dendo, k = 7),
                 #clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
                 #column_order = order(as.numeric(gsub("column", "", colnames(HExp_el_Bulkm)))),
                 show_row_names = FALSE,
                 split=split,
                 cluster_row_slices = T,
                 width = unit(6, "cm") , height = unit(12, "cm"),
                 cluster_columns = fh,
                 right_annotation=ha,
                 use_raster = T
)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(heat, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(heat2, newpage = FALSE)
upViewport()

upViewport()



-------------
  TOP 10
------------
  
  
  
top10_el_DF <- read.csv("TOP_10_elife.csv", header = T)
top10_el <- AllBulkscaled[which(row.names(AllBulkscaled) %in% top10_el_DF$GeneID),]
top10_el <- na.omit(top10_el) 
write.csv(top10_el, file= "HExp_el_Bulk_TOP10.csv
          ")
top10_elm = top10_el[3:14] %>% as.matrix()
#HExp_el_Bulkm = HExp_el_Bulk[c(5, 6, 7, 8, 10, 11, 13, 14) ] %>% as.matrix()



row_dendo = hclust(dist(top10_elm), method = "complete") # row clustering
kclus <- kmeans(top10_elm,7)
split2 <- kclus$cluster
split <- top10_el$Cluster



cc <- data.frame(top10_el$Cluster)
colnames(cc)<- c("Cluster")
colorsP <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#f781bf")
nr_cluster=sort(unique(cc$Cluster))
names(colorsP) = nr_cluster
colz = list(Cluster = colorsP)
ha = HeatmapAnnotation(df = cc, which = "row", width = unit(1, "cm"), col = colz)



fh = function(x) fastcluster::hclust(dist(x))

heat <- Heatmap(top10_elm,
                column_title = "Top10 Exp Genes from Elife markerset",
                name = "z-score", #title of legend
                #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                col = viridis(5),
                #cluster_rows = color_branches(row_dendo, k = 7),
                #clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
                #column_order = order(as.numeric(gsub("column", "", colnames(HExp_el_Bulkm)))),
                show_row_names = FALSE,
                split=split2,
                cluster_row_slices = T,
                width = unit(6, "cm") , height = unit(12, "cm"),
                cluster_columns = fh,
                right_annotation=ha,
                use_raster = T
)
heat2 <- Heatmap(top10_elm,
                 name = "z-score", #title of legend
                 #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                 col = viridis(5),
                 #cluster_rows = color_branches(row_dendo, k = 7),
                 #clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
                 #column_order = order(as.numeric(gsub("column", "", colnames(HExp_el_Bulkm)))),
                 show_row_names = T,
                 split=split,
                 cluster_row_slices = T,
                 width = unit(6, "cm") , height = unit(12, "cm"),
                 cluster_columns = fh,
                 right_annotation=ha,
                 use_raster = T
)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(heat, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(heat2, newpage = FALSE)
upViewport()

upViewport()





top10_SC_DF <- read.csv("TOP_10_SC.csv", header = T)
top10_SC <- AllBulkscaled[which(row.names(AllBulkscaled) %in% top10_SC_DF$GeneID),]
top10_SC <- na.omit(top10_SC) 
write.csv(HExp_el_Bulk, file= "HExp_SC_Bulk_TOP10.csv
          ")
top10_SCm = top10_SC[3:14] %>% as.matrix()
#HExp_el_Bulkm = HExp_el_Bulk[c(5, 6, 7, 8, 10, 11, 13, 14) ] %>% as.matrix()



row_dendo = hclust(dist(top10_SCm), method = "complete") # row clustering
kclus <- kmeans(top10_SCm,7)
split2 <- kclus$cluster
split <- top10_SC$Cluster



cc <- data.frame(top10_SC$Cluster)
colnames(cc)<- c("Cluster")
colorsP <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#f781bf")
nr_cluster=sort(unique(cc$Cluster))
names(colorsP) = nr_cluster
colz = list(Cluster = colorsP)
ha = HeatmapAnnotation(df = cc, which = "row", width = unit(1, "cm"), col = colz)



fh = function(x) fastcluster::hclust(dist(x))

heat <- Heatmap(top10_SCm,
                column_title = "Top10 Exp Genes from SC based on Elife markerset",
                name = "z-score", #title of legend
                #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                col = viridis(5),
                #cluster_rows = color_branches(row_dendo, k = 7),
                #clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
                #column_order = order(as.numeric(gsub("column", "", colnames(HExp_el_Bulkm)))),
                show_row_names = T,
                split=split2,
                cluster_row_slices = T,
                width = unit(6, "cm") , height = unit(12, "cm"),
                cluster_columns = fh,
                right_annotation=ha,
                use_raster = T,
                row_names_gp = gpar(fontsize = 4, fontfamily = "sans", fontface = "bold"),
                column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                column_names_gp = grid::gpar(fontsize = 9)
)
heat2 <- Heatmap(top10_SCm,
                 name = "z-score", #title of legend
                 #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
                 col = viridis(5),
                 #cluster_rows = color_branches(row_dendo, k = 7),
                 #clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
                 #column_order = order(as.numeric(gsub("column", "", colnames(HExp_el_Bulkm)))),
                 show_row_names = T,
                 split=split,
                 cluster_row_slices = T,
                 width = unit(6, "cm") , height = unit(12, "cm"),
                 cluster_columns = fh,
                 right_annotation=ha,
                 use_raster = T,
                 row_names_gp = gpar(fontsize = 4, fontfamily = "sans", fontface = "bold"),
                 column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                 column_names_gp = grid::gpar(fontsize = 9)
)
pdf('Top10Genes.pdf', width = 12, height = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(heat, newpage = FALSE)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(heat2, newpage = FALSE)
upViewport()

dev.off()

