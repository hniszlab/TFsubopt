library(fgsea)
library(tidyverse)

GOBP_MYOBLAST_FUSION <- gmtPathways("GOBP_MYOBLAST_FUSION.v2022.1.Mm.gmt")
GOBP_MYOBLAST_DIFFERENTIATION <- gmtPathways("GOBP_MYOBLAST_DIFFERENTIATION.v2022.1.Mm.gmt")
GOBP_CELL_CELL_ADHESION <- gmtPathways("GOBP_CELL_CELL_ADHESION.v2022.1.Mm.gmt")
AroPC_WT <- read_csv("AroPC_WT.csv")
AroP_WT <- read_csv("AroP_WT.csv")

ranksAroPC_WT = AroPC_WT$stat
names(ranksAroPC_WT) = AroPC_WT$Gene_name

wdup = which(duplicated(names(ranksAroPC_WT)))
if (length(wdup) > 0) ranksAroPC_WT = ranksAroPC_WT[-wdup]

ranksAroP_WT = AroP_WT$stat
names(ranksAroP_WT) = AroP_WT$Gene_name

wdup2 = which(duplicated(names(ranksAroP_WT)))
if (length(wdup2) > 0) ranksAroP_WT = ranksAroP_WT[-wdup2]


fgseaRes2C <-fgseaSimple(GOBP_MYOBLAST_FUSION[1], ranksAroPC_WT, minSize=15, nperm = 1000)
fgseaRes3C <-fgseaSimple(GOBP_MYOBLAST_DIFFERENTIATION[1], ranksAroPC_WT, minSize=15, nperm = 1000)
fgseaRes4C <-fgseaSimple(GOBP_CELL_CELL_ADHESION[1], ranksAroPC_WT, minSize=15, nperm = 1000)

fgseaRes2 <-fgseaSimple(GOBP_MYOBLAST_FUSION[1], ranksAroP_WT, minSize=15, nperm = 1000)
fgseaRes3 <-fgseaSimple(GOBP_MYOBLAST_DIFFERENTIATION[1], ranksAroP_WT, minSize=15, nperm = 1000)
fgseaRes4 <-fgseaSimple(GOBP_CELL_CELL_ADHESION[1], ranksAroP_WT, minSize=15, nperm = 1000)


