
library('tidyverse')
#library('DESeq2')
library('EnhancedVolcano')
#library(ComplexHeatmap)
library(RColorBrewer)
library(dendextend)
library(viridis)
library(gridExtra)
library(grid)

AroPLUSGFPvsWTGFP = read.csv("AroPLUSGFPvsWTGFP.complete.txt.csv" , header = T)
AroPERFECTGFPvsWTGFP = read.csv("AroPERFECTGFPvsWTGFP.complete.txt.csv" , header = T)

hgene = c("HOXD1", "HOXD3", "HOXD9", "HOXD10", "HOXD11", "HOXD13", "HOXD4", "HOXD12", "HOXD8")

DEGsKOvsWTGFP = read.csv("KOvsWTGFP.txt", header = T, sep = "\t")
listDegs = unlist(DEGsKOvsWTGFP$Id, use.names = FALSE)
names(listDegs) <- unlist(DEGsKOvsWTGFP$color)

keyvals <- ifelse(
  KOvsWTGFP$log2FoldChange < -2 & KOvsWTGFP$padj < 0.01, 'blue',
  ifelse(KOvsWTGFP$log2FoldChange > 2 & KOvsWTGFP$padj < 0.01, 'blue',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == 'blue'] <- 'DEGKOvsWTGFP'
names(keyvals)[keyvals == 'gray'] <- ''
names(keyvals)[keyvals == 'blue'] <- 'DEGKOvsWTGFP'



# p1 <- EnhancedVolcano(KOvsWTGFP,
#                       lab = KOvsWTGFP$Id,
#                       x = 'log2FoldChange',
#                       y = 'padj',
#                       title = 'KO versus WTGFP',
#                       subtitle = NULL,
#                       caption = NULL,
#                       titleLabSize = 10,
#                       pCutoff = 0.01,
#                       FCcutoff = 2,
#                       pointSize = 1.0,
#                       labSize = 0,
#                       colCustom = keyvals,
#                       arrowheads = FALSE,
#                       col = c('black', 'black', 'black', 'red3'),
#                       axisLabSize = 10,
#                       drawConnectors = F,
#                       widthConnectors = 0.75,
#                       legendPosition = 'top',
#                       legendLabSize = 8,
#                       legendIconSize = 2.0,
#                       border = 'full',
#                       boxedLabels = F)
p2 <- EnhancedVolcano(AroPLUSGFPvsWTGFP,
                      lab = AroPLUSGFPvsWTGFP$Id,
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = 'KO versus WTGFP',
                      subtitle = NULL,
                      caption = NULL,
                      titleLabSize = 10,
                      pCutoff = 0.01,
                      FCcutoff = 2,
                      pointSize = 1.0,
                      labSize = 0,
                      colCustom = keyvals,
                      arrowheads = FALSE,
                      col = c('gray', 'gray', 'gray', 'red3'),
                      axisLabSize = 10,
                      drawConnectors = F,
                      widthConnectors = 0.75,
                      legendPosition = 'top',
                      legendLabSize = 8,
                      legendIconSize = 2.0,
                      border = 'full',
                      boxedLabels = F)
p3 <- EnhancedVolcano(AroPERFECTGFPvsWTGFP,
                      lab = AroPERFECTGFPvsWTGFP$Id,
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = 'KO versus WTGFP',
                      subtitle = NULL,
                      caption = NULL,
                      titleLabSize = 10,
                      pCutoff = 0.01,
                      FCcutoff = 2,
                      pointSize = 1.0,
                      labSize = 0,
                      colCustom = keyvals,
                      arrowheads = FALSE,
                      col = c('gray', 'gray', 'gray', 'red3'),
                      axisLabSize = 10,
                      drawConnectors = F,
                      widthConnectors = 0.75,
                      legendPosition = 'top',
                      legendLabSize = 8,
                      legendIconSize = 2.0,
                      border = 'full',
                      boxedLabels = F)


hgene = c("HOXD1", "HOXD3", "HOXD9", "HOXD10", "HOXD11", "HOXD13", "HOXD4", "HOXD12", "HOXD8")
keyvals <- ifelse(
  KOvsWTGFP$GeneSymbol == hgene, 'blue',
  'gray')
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == 'blue'] <- 'HOXD Genes'
names(keyvals)[keyvals == 'gray'] <- ''

# h1 <- EnhancedVolcano(KOvsWTGFP,
#                       lab = KOvsWTGFP$GeneSymbol,
#                       x = 'log2FoldChange',
#                       y = 'padj',
#                       title = 'KO versus WTGFP',
#                       subtitle = NULL,
#                       caption = NULL,
#                       titleLabSize = 10,
#                       pCutoff = 0.01,
#                       FCcutoff = 2,
#                       pointSize = 1.0,
#                       labSize = 2.0,
#                       axisLabSize = 10,
#                       selectLab = hgene,
#                       drawConnectors = TRUE,
#                       widthConnectors = 0.75,
#                       legendPosition = 'top',
#                       legendLabSize = 8,
#                       legendIconSize = 2.0,
#                       border = 'full',
#                       boxedLabels = TRUE)
h2 <- EnhancedVolcano(AroPLUSGFPvsWTGFP,
                      lab = AroPLUSGFPvsWTGFP$GeneSymbol,
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = 'AroPlus versus WTGFP',
                      subtitle = NULL,
                      caption = NULL,
                      titleLabSize = 10,
                      pCutoff = 0.01,
                      FCcutoff = 2,
                      pointSize = 1.0,
                      labSize = 2.0,
                      axisLabSize = 10,
                      selectLab = hgene,
                      drawConnectors = TRUE,
                      widthConnectors = 0.75,
                      legendPosition = 'top',
                      legendLabSize = 8,
                      legendIconSize = 2.0,
                      border = 'full',
                      boxedLabels = TRUE)
h3 <- EnhancedVolcano(AroPERFECTGFPvsWTGFP,
                      lab = AroPERFECTGFPvsWTGFP$GeneSymbol,
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = 'AroPerfect versus WTGFP',
                      subtitle = NULL,
                      caption = NULL,
                      titleLabSize = 10,
                      axisLabSize = 10,
                      pCutoff = 0.01,
                      FCcutoff = 2,
                      pointSize = 1.0,
                      labSize = 2.0,
                      selectLab = hgene,
                      drawConnectors = TRUE,
                      widthConnectors = 0.75,
                      legendPosition = 'top',
                      legendLabSize = 8,
                      legendIconSize = 2.0,
                      border = 'full',
                      boxedLabels = TRUE)



