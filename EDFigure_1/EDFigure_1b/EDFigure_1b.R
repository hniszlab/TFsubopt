library(tidyverse)
library(ggthemes)
library(ggpubr)
library(bbplot)
library(circlize)


pblocks_df <- read_csv("Periodic_blocks_HG38_TFs_0-3_4-10_11-20_21-30_threshold_4_20231009.csv")

h <- hist(pblocks_df1$FYW, breaks=10)
text(h$mids,h$counts,labels=h$counts, adj=c(0.5, -0.5))

pblocks_df$FYW<-as.factor(pblocks_df$FYW)

df <- pblocks_df %>%
  group_by(FYW) %>%
  summarise(counts = n())
df

gghistogram(pblocks_df, x = "FYW",
            fill = "#acd7e6", color="black",
            xlab = "Number of aromatic residues per block", 
            ylab = "Number of periodic blocks",
            stat="count", rug = F ,
            position = "identity") + 
  geom_text(aes(x = factor(FYW),
                y = 420,
                label = paste("n = ", df$counts,"\n")),
            aggregate(. ~ FYW ,pblocks_df,length),
            position = position_dodge(.8))

ggsave("ED1b.pdf")

