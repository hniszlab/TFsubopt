library(tidyverse)
library(ggthemes)
library(ggpubr)
library(bbplot)
library(circlize)


pblocks_df <- read_csv("Periodic_blocks_HG38_TFs_0-3_4-10_11-20_21-30_threshold_4_20231009.csv")


df <- pblocks_df[1:80, ]
df <- df %>% arrange(GeneFamily, desc(Score))
#df <- df %>% rowid_to_column("ID")
df$GeneName <- factor(df$GeneName,levels = rev(unique(df$GeneName)), ordered=T)


df$gene <- df$GeneName

GeneFamilyS_v = c("Homeobox"="Homeobox", "bHLH"="bHLH", "RHD"="RHD", "T-box"="T-box", "TF_Otx"="Otx","TF_bZIP"="bZIP", "HSF"="others", "AP-2"="others", "DM"="others", "ESR-like"="others", "ETS"="others", "GCM"="others", "HMG"="others", "MYB"="others", "NGFIB-like"="others", "Pou"="others", "RFX"="others", "ZBTB"="others", "zf-C2H2"="C2H2-ZF")

df$GeneFamilyS <- GeneFamilyS_v[df$GeneFamily]

df <- df %>% group_by(GeneFamilyS) %>% arrange(desc(Score))

mat1 <- df %>% column_to_rownames(var = "GeneName")

split <- dplyr::pull(mat1, GeneFamily)
split2 <- dplyr::pull(mat1, GeneFamilyS)
score <- dplyr::pull(mat1, Score)


col_fun1 = colorRamp2(c(-1, 0, 1), c("blue", "white", "#e41a1c"))
col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "#377eb8"))
col_fun3 = colorRamp2(c(-1, 0, 1), c("blue", "white", "#4daf4a"))
col_fun4 = colorRamp2(c(-1, 0, 1), c("blue", "white", "#984ea3"))
col_fun5 = colorRamp2(c(-1, 0, 1), c("blue", "white", "#ff7f00"))

GeneFamily_list = unique(split)

group_color = c("Homeobox"="#276B7D", "bHLH"="#c0524f", "RHD"="#6699cc", "T-box"="#99cc66", "TF_Otx"="#FF9933",  "Otx"="#FF9933","TF_bZIP"="#29235C","bZIP"="#29235C", "HSF"="#919396", "AP-2"="#919396", "DM"="#919396", "ESR-like"="#919396", "ETS"="#919396", "GCM"="#919396", "HMG"="#919396", "MYB"="#919396", "NGFIB-like"="#919396", "Pou"="#919396", "RFX"="#919396", "ZBTB"="#919396", "zf-C2H2"="#8064A1","C2H2-ZF"="#8064A1", "others"="#919396" )

#group_color = c("#919396", "#c0524f", "#919396", "#919396", "#919396", "#919396", "#919396", "#276B7D", "#919396", "#919396", "#919396", "#919396", "#919396", "#919396", "#29235C", "#FF9933", "#919396", "#919396", "#8064A1")
#names(group_color) <- GeneFamily_list


mat1$color2 <- group_color[mat1$GeneFamilyS]
mat1$color <- group_color[mat1$GeneFamily]




mat1 = mat1 %>% arrange_at(3, desc) %>% arrange(match(GeneFamilyS, c("Homeobox", "bHLH", "RHD", "T-box", "TF_Otx", "Otx","TF_bZIP","bZIP", "HSF", "AP-2", "DM", "ESR-like", "ETS", "GCM", "HMG", "MYB", "NGFIB-like", "Pou", "RFX", "ZBTB", "others", "zf-C2H2", "C2H2-ZF")))

mat1$GeneFamilyS <- factor(mat1$GeneFamilyS , levels = c("Homeobox", "bHLH", "RHD", "T-box", "Otx","bZIP", "others","C2H2-ZF"))


split <- dplyr::pull(mat1, GeneFamily)
split2 <- dplyr::pull(mat1, GeneFamilyS)
score <- dplyr::pull(mat1, Score)

labels <- row.names(mat1)

mat1 = mat1 %>% mutate(row_id = row_number())

colorvector <- dplyr::pull(mat1, color)
names(colorvector) <- dplyr::pull(mat1, Score)

colorvector2 <- dplyr::pull(mat1, color2)
names(colorvector2) <- dplyr::pull(mat1, Score)


colorvector3 <- dplyr::pull(mat1,color2)
names(colorvector3) <- dplyr::pull(mat1, row_id)


circos.clear()
circos.par(start.degree =0, gap.degree = 2,cell.padding = c(0, 0, 0, 0))
circos.heatmap.initialize(mat1, cluster = FALSE, split = split2)
circos.track(ylim = c(0, 30), track.height = 0.3, bg.border = NA, panel.fun = function(x, y) {
  x = mat1[CELL_META$subset, 2]
  circos.barplot(x, pos = CELL_META$cell_middle, 
                 bar_width = CELL_META$cell_width*0.9)
})
pos = circos.heatmap.get.x(c(1,2,3,8,12,19,22,25,28,29,40, 42,50,68))
circos.clear()
circos.par(start.degree =0, gap.degree = 2,cell.padding = c(0, 0, 0, 0))
circos.heatmap.initialize(mat1, cluster = FALSE, split = split2)
circos.labels(pos[, 1], x = pos[, 2], labels = mat1$gene[pos[, 3]], side = "outside", line_lwd = 0)
circos.track(ylim = c(0, 30), track.height = 0.3, bg.border = NA, panel.fun = function(x, y) {
  x = mat1[CELL_META$subset, 2]
  circos.barplot(x, pos = CELL_META$cell_middle, 
                 bar_width = CELL_META$cell_width*0.9, col= colorvector3, border = NA)
})

circos.heatmap(mat1[1], col = group_color, track.height = 0.02)
#
set_track_gap(cm_h(0))
circos.heatmap(mat1[9], col = col_fun1, track.height = 0.03, cluster = FALSE, split = split,
               cell.border = "black", cell.lwd = 0.3, cell.lty = 1)
set_track_gap(cm_h(0))
circos.heatmap(mat1[10], col = col_fun2, track.height = 0.03, cluster = FALSE, split = split,
               cell.border = "black", cell.lwd = 0.3, cell.lty = 1)
set_track_gap(cm_h(0))
circos.heatmap(mat1[12], col = col_fun3, track.height = 0.03, cluster = FALSE, split = split,
               cell.border = "black", cell.lwd = 0.3, cell.lty = 1)
set_track_gap(cm_h(0))
circos.heatmap(mat1[13], col = col_fun5, track.height = 0.03, cluster = FALSE, split = split,
               cell.border = "black", cell.lwd = 0.3, cell.lty = 1)
