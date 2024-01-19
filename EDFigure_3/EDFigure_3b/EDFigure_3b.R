library(MASS)
library(viridis)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(RColorBrewer)
#library(Seurat)
library(ggrepel)
library(readxl)

my_data_metapredict <- read_csv("MasterTable_qIDR_17042022.csv")

thre=0.01


## 1) top qIDRs
my_data_metapredict_top=my_data_metapredict[my_data_metapredict$min_pval < 0.01, ]
#my_data_metapredict_top$min_pval[which(my_data_metapredict_top$min_pval < thre)]=thre

highlight_df_3 <- my_data_metapredict_top %>% 
  filter(my_data_metapredict_top$overlap_metapredict  == "yes" )

blues <- brewer.pal(9, "Blues")
blue_range <- colorRampPalette(blues)

genes.to.label <- c("DAZ1", "EWSR1", "HNRNPA1", "EGR1")
genes.to.label2 <- my_data_metapredict_top %>% filter(gene_symbol %in% genes.to.label)

ggplot(my_data_metapredict_top, aes(x=-log10(min_pval), y= length_IDR_qIDR) ) +
  geom_point(data=my_data_metapredict_top,  aes(x=-log10(min_pval),y=length_IDR_qIDR), 
             color='black', size=0.5, alpha = 0.5) +
  geom_density2d_filled(bins=9, contour = TRUE, alpha = 0.9) + 
  theme_classic() +
  scale_fill_manual(values = blue_range(9)) + 
  theme(aspect.ratio=1,
        legend.position="none" ,
        legend.title = element_text( size=9),
        legend.text=element_text(size=7)
  ) +
  coord_cartesian(xlim = c(NA,5),ylim = c(100,300)) +
  xlab(expression(-log[10]~(p-value))) +
  ylab("Length (amino acids)") + # for the x axis label
  geom_point(data=genes.to.label2,  aes(x=-log10(min_pval),y=length_IDR_qIDR), 
             color='black', size=0.5, alpha = 1) +
  geom_label_repel(data = genes.to.label2, aes(label = gene_symbol),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   min.segment.length = 0
  )
ggsave("too_qIDR_plot_pval_.pdf")