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

##all qIDRs density
g <- ggplot(highlight_df_3, aes(x=-log10(min_pval), y= length_IDR_qIDR) ) + 
  geom_density2d_filled(bins=9) +
  coord_cartesian(xlim = c(NA,5), ylim = c(100,300)) +
  scale_fill_manual(values = blue_range(9)) +
  theme_classic() +
  ggtitle("All IDRs with significant periodicit") +
  theme(aspect.ratio=1,
        plot.title = element_text(size = 6),
        legend.position="none",
        legend.title = element_text( size=9),
        legend.text=element_text(size=7),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) + 
  xlab(expression(-log[10]~(p-value))) +
  ylab("Length (amino acids)") # for the x axis label




g2 <- ggplot(highlight_df_TFs, aes(x=-log10(min_pval), y= length_IDR_qIDR) ) + 
  geom_density2d_filled(bins=9) +
  coord_cartesian(xlim = c(NA,5), ylim = c(100,300)) +
  scale_fill_manual(values = greens_range(9)) +
  theme_classic() +
  ggtitle("Transcription Factors") +
  theme(aspect.ratio=1,
        plot.title = element_text(size = 6),
        legend.position="none",
        legend.title = element_text( size=9),
        legend.text=element_text(size=7),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) + 
  xlab(expression(-log[10]~(p-value))) +
  ylab("Length (amino acids)") # for the x axis label



g3 <- ggplot(highlight_df_PLDs, aes(x=-log10(min_pval), y= length_IDR_qIDR) ) + 
  geom_density2d_filled(bins=9) +
  coord_cartesian(xlim = c(NA,5), ylim = c(100,300)) +
  scale_fill_manual(values = purple_range(9)) +
  theme_classic() +
  ggtitle("Prion-like domains") +
  theme(aspect.ratio=1,
        plot.title = element_text(size = 6),
        legend.position="none",
        legend.title = element_text( size=9),
        legend.text=element_text(size=7),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) + 
  xlab(expression(-log[10]~(p-value))) +
  ylab("Length (amino acids)") # for the x axis label

g4 <- ggplot(highlight_df_PLDs99, aes(x=-log10(min_pval), y= length_IDR_qIDR) ) + 
  geom_density2d_filled(bins=9) +
  coord_cartesian(xlim = c(NA,5), ylim = c(100,300)) +
  scale_fill_manual(values = red_range(9)) +
  theme_classic() +
  ggtitle("Aromatic-rich prion-like domains") +
  theme(aspect.ratio=1,
        plot.title = element_text(size = 6),
        legend.position="none",
        legend.title = element_text( size=9),
        legend.text=element_text(size=7),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) + 
  xlab(expression(-log[10]~(p-value))) +
  ylab("Length (amino acids)") # for the x axis label
```
```{r}
plot_grid(g, g2, g3, g4, ncol = 2, nrow = 2 , align = "hv")
ggsave("dencityPlot_qIDR_Final.pdf")


