library(tidyverse)
library(ggplotgui)
library(plotly)
library(gridExtra)
library(ggpubr)
library(tidyverse)
library(agricolae)
library(cowplot)
library(multcompView)
library(car)
library(ggsignif)
library(multcomp)
library(ggrepel)
library(Seurat)

df <- read_csv("HG38_MetaIDR_qIDR_wOmega_wMoR.csv")
df <- read_csv("OmegaScore_df.csv")

keys <- read.csv(file = 'HG30_pep_annotation.csv')
TFnames <- read.table("TF_IDs.txt",header = F)
TFfamilies <- read_csv("AnnotationTFs.csv")
PLDs <- read.table("PLD_IDs.txt",header = F)

RBP <- read_csv("RBP_HS_keys.csv")
target <- c("zf-BED", "zf-C2H2", "zf-C2HC", "zf-CCCH", "zf-GATA", "zf-LITAF-like", "zf-MIZ", "zf-NF-X1")
ZF_keys= TFfamilies %>% filter(TFfamilies$Type  %in% target)
ZF_keys= as.vector(ZF_keys$Name)
TFfamilies <- TFfamilies[,c("EnsenbleProteinID","Type")]

names(TFfamilies)[1] <- 'PepID'
names(TFfamilies)[2] <- 'TF_Family'


df <- left_join(df, TFfamilies, by = "PepID")

#convert PLD keys to match
PLD_keys=keys$gene_symbol[match(PLDs$V1, keys$ID )]
#convert RBP keys to match
RPBs_keys=toupper(tail(as.vector(RBP$UNIQUE), -1))


df$Genetype <- ""
df$Genetype2 <- ""
df$Genetype3 <- ""
df$Genetype4 <- ""
df$Genetype[df$Gene_symbol  %in% PLD_keys] <- "PLD"
df$Genetype2[df$Gene_symbol  %in% PLD99_keys] <- "PLD99"
df$Genetype3[df$Gene_symbol  %in% TFnames$V1] <- "TF"

qidr <- read_csv("MasterTable_qIDR_17042022_wOmega_wO_wMoR100.csv")

qidr001 <- qidr[qidr$min_pval < 0.01, ]

qidr001=as.vector(qidr001$ID)

df001 <- dff %>% filter(dff$ID  %in% qidr001 | grepl("_IDR",ID) )


labels_df <- df001 %>% subset(Gene_symbol %in% labels)

labels <- c("CEBPA", "HOXD4", "HOXB1", "HOXB1", "HOXB1", "HOXB1", "HOXB1", "HOXB1", "HOXC4", "NEUROG2", "MYOD1", "MYOD1")
points <- c("9264", "3010", "162", "172", "181", "191", "12125", "12135", "8475", "3679", "938", "938")

publishedTF <- c("NEUROG1", "ATOH1", "JDP2", "AFF1", "TERF1", "HLF", "MEIS1", "GCM2", "CREBL2", "GATA5", "KLF7", "ASCL1", "AKNA", "NKX3-1", "GBX1", "PRDM3", "NEUROG2", "MEOX1", "SOX9", "ATF6", "NR1H3", "ELF5", "SOX2", "PAX7", "CBFB", "GLI1", "SOX1", "AIRE", "SOX6", "RELA", "LMX1A", "GABPA", "BACH2", "HES2", "NKX6_2", "HOXA10", "ID1", "MESP1", "FLI1", "NEUROG3", "SOX4", "PRRX1", "DLX3", "HHEX", "GFI1B", "MYB", "HOXA9", "DLX5", "NANOG")

TFRepro <- c("NEUROG1", "ATOH1", "ZNF700", "ZNF541", "ZNF616", "ZNF85", "ZNF614", "ZNF93", "EPAS1", "ZNF69", "DMRTA1", "ZNF440", "ZNF107", "ZNF684", "HESX1", "ZNF626", "THAP1", "ZNF718", "ASCL3", "JDP2", "ZNF239", "ZNF253", "E2F3", "ZBTB7A", "ZNF695", "BBX", "DDIT3", "HOXA6", "ZNF83", "AR", "ZNF367", "TSC22D3", "GON4L", "NKX1-2", "ZNF114", "AFF1", "TERF1", "HLF", "MEIS1", "ZBTB7C", "HOXB13", "ZNF593", "ZNF710", "TAL2", "TFAP2A", "ZNF667", "CREB1", "GCM2", "ZNF816", "CREBL2", "ZNF491", "GATA5", "MLXIP", "DMRT3", "ZZZ3", "ZNF585A", "KLF7", "ZBTB8B", "ZBTB40", "TFDP2", "ZNF256", "TFEB", "TFDP3", "ZNF320", "ATF7", "ZNF492", "ZNF98", "ASCL1", "CPXCR1", "KIAA0961", "AKNA", "ZNF28", "NKX3-1", "ZNF26", "ZNF512", "FOXJ1", "ZNF468", "GBX1", "ZNF141", "ZNF766", "ZNF655", "ZNF92", "ZNF135", "PRDM3", "MAX", "THAP6", "ETS1", "ATF2", "CEBPE", "NEUROG2", "MEOX1", "ZNF562", "ZNF547", "TERF2", "DRAP1", "STAT6", "ZNF570", "SOX9", "OTP", "RFX7", "C13ORF8", "ZNF160", "ZBTB6", "ZSCAN25", "ZIC3", "ATF6", "NR1H3", "TFE3", "NFIX", "ZNF260", "ELF5", "ZNF148", "SOX2", "MIER1", "ZBTB1", "FERD3L", "BOLA3", "DR1", "ZNF419", "PAX7", "ZMAT4", "ZHX3", "CBFB", "ZNF187", "POU1F1", "ZSCAN26", "EBF4", "GLI1", "SOX1", "ZBTB42", "THAP3", "AIRE", "HOXC9", "BATF", "NFYB", "SOX6", "ZNF827", "ZNF44", "RELA", "HOXD3", "FOXN2", "LMX1A", "FAM170A", "BATF2", "PRDM7", "ZNF697", "HOXB6", "IRX5", "IKZF5", "ZNF536", "OTX1", "SOX14", "NR1C1", "ZNF33B", "VSX1", "GABPA", "BACH2", "ZNF34", "GRHL2", "EGR2", "EVX2", "ZNF302", "HES2", "ZIK1", "ZNF833P", "BHLHB4", "ZNF607", "HMGA1", "ZNF333", "ZNF236", "NKX6_1", "ZNF560", "LOC91661", "DPF1", "HSFY1", "RFXANK", "ZBTB38", "ZNF99", "ZNF793", "NR1C2", "NR1H4", "SMAD4", "HSFY2", "HOXD8", "CREBZF", "NFAT5", "NKX6_2", "ZSCAN1", "ZNF17", "HOXA10", "DBX1", "ZC3H3", "PREB", "DMRTA2", "HOXB9", "NR2B3", "SPIC", "PLEK", "ISX", "TEAD4", "ZNF311", "ZNF705D", "UBP1", "ID1", "MESP1", "USF1", "CREB5", "DPRX", "ZNF703", "POU4F2", "ETV3L", "OTX2", "FLI1", "ZNF687", "CDX1", "FOSB", "BHLHB2", "TFEC", "ZNF138", "NEUROG3", "ELF4", "MYOG", "NEUROD4", "BOLA2", "ZNF324", "ZFHX2", "FOS", "SOX4", "ZNF324B", "ZBTB3", "IRF1", "ZNF200", "TEAD2", "PRRX1", "MAFK", "ALX4", "FLJ35784", "ZNF358", "ZFP92", "ZNF705A", "ZSCAN32", "ZNF434", "PCGF6", "HMBOX1", "ZNF768", "ZBTB44", "ZNF22", "BOLA2B", "MEIS3", "DLX3", "SMAD8", "THAP8", "ZSCAN5B", "MAZ", "LOC51058", "ZNF449", "PITX2", "HHEX", "ID2", "TLX2", "ZNF540", "ZNF600", "ZNF678", "ZNF692", "MEOX2", "HES1", "SPDEF", "GFI1B", "NR0B1", "PBX4", "MYB", "ETV6", "MITF", "ZNF317", "NFIC", "ZNF701", "ZNF553", "HOXA9", "POU5F2", "ZKSCAN10", "ZNF726", "DLX5", "ZNF644", "ZIM2", "ZNF774", "MYF6", "FOXD2", "FOXD1", "NANOG", "SP140L")



# reorder is close to order, but is made to change the order of the factor levels.
colors <- c("#6355a2", "#cb181d", "#478949", "#617d8b")
graph <- ggplot(df001[!is.na(df001$omega) & df001$omega < 1 & df001$IDRfrac_aromatic >= 0.1,], aes(x = Gtype, y = omega, colour = Gtype , label = Gene_symbol)) +
  #geom_point(position=position_jitterdodge(),size = 0.5, alpha = 0.2) +
  geom_jitter(shape=16, position=position_jitter(0.2),size = 0.5, alpha = 0.2) +
  geom_text_repel(data = labels_df, aes(label= Gene_symbol), check_overlap = T, box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0)  +
  geom_violin(notch = FALSE, alpha = 0.5,trim=FALSE) +
  geom_boxplot(width=0.1, outlier.shape = NA)+
  facet_grid(. ~ Type) +
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1
  ) +
  scale_colour_manual(values = colors) +
  ylim(0, 1)
graph

graph2 <- ggplot(df001[!is.na(df001$omega) & df001$omega < 1 & df001$IDRfrac_aromatic >= 0.1 & df001$Gene_symbol %in% TFRepro,], aes(x = Gtype, y = omega, colour = Gtype , label = Gene_symbol)) +
  #geom_point(position=position_jitterdodge(),size = 0.5, alpha = 0.2) +
  geom_jitter(shape=16, position=position_jitter(0.2),size = 0.5, alpha = 0.2) +
  geom_text_repel(data = labels_df, aes(label= Gene_symbol), check_overlap = T, box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0)  +
  geom_violin(notch = FALSE, alpha = 0.5,trim=FALSE) +
  geom_boxplot(width=0.1, outlier.shape = NA)+
  facet_grid(. ~ Type) +
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1
  ) +
  scale_colour_manual(values = colors) +
  ylim(0, 1)
graph2

graph3 <- ggplot(df001[!is.na(df001$omega) & df001$omega < 1 & df001$IDRfrac_aromatic >= 0.1 & df001$Gene_symbol %in% publishedTF,], aes(x = Gtype, y = omega, colour = Gtype , label = Gene_symbol)) +
  #geom_point(position=position_jitterdodge(),size = 0.5, alpha = 0.2) +
  geom_jitter(shape=16, position=position_jitter(0.2),size = 0.5, alpha = 0.2) +
  geom_text_repel(data = labels_df, aes(label= Gene_symbol), check_overlap = T, box.padding = 0.5, max.overlaps = Inf, min.segment.length = 0)  +
  geom_violin(notch = FALSE, alpha = 0.5,trim=FALSE) +
  geom_boxplot(width=0.1, outlier.shape = NA)+
  facet_grid(. ~ Type) +
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1
  ) +
  scale_colour_manual(values = colors) +
  ylim(0, 1)
graph3

#LabelPoints(plot = graph, points = points, labels = labels)

# graph3 <- ggplot(df001[!is.na(df001$omega) & df001$omega < 1 & df001$IDRfrac_aromatic >= 0.1,], aes(x = Gtype, y = Mean_of_random, colour = Gtype, label = Gene_symbol)) +
#   #geom_point(position=position_jitterdodge(),size = 0.5, alpha = 0.2) +
#   geom_jitter(shape=16, position=position_jitter(0.2),size = 0.5, alpha = 0.2) +
#     geom_violin(notch = FALSE, alpha = 0.5,trim=FALSE) +
#   geom_boxplot(width=0.1, outlier.shape = NA)+
#     facet_grid(. ~ Type) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1
#   )+
#   scale_colour_manual(values = colors) +
#    ylim(0, 1)
# 
# graph5 <- ggplot(df001[!is.na(df001$omega) & df001$omega < 1 & df001$IDRfrac_aromatic >= 0.1,], aes(x = Gtype, y = omega-Mean_of_random, colour = Gtype, label = Gene_symbol)) +
#   #geom_point(position=position_jitterdodge(),size = 0.5, alpha = 0.2) +
#   geom_jitter(shape=16, position=position_jitter(0.2),size = 0.5, alpha = 0.2) +
#     geom_violin(notch = FALSE, alpha = 0.5,trim=FALSE) +
#   geom_boxplot(width=0.1, outlier.shape = NA)+
#     facet_grid(. ~ Type) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio=1
#   ) +
#   scale_colour_manual(values = colors)
ggarrange(graph,graph3,graph2,  ncol = 3, nrow = 1, common.legend = T)
ggsave("OmegaBoxPlot_violin_22Feb2023.pdf")
