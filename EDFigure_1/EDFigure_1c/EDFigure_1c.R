``{r}
library(drawProteins)
library(Biostrings)
library(tidyverse)


df_TAD <- read_csv("TAD_database.csv")
df_PB <- read.csv("Coord_block.csv")

df <- read_csv("domains.csv")

df$begin <- df$begin * -1
df$end <- df$end * -1

df_chain <- df %>% dplyr::filter(type == "CHAIN")
df_DNA_BIND <- df %>% dplyr::filter(type == "DOMAIN" | type == "DNA_BIND")
df_Peri <- df %>% dplyr::filter(type == "Periodic")


df_chain <- df_chain %>% arrange(-Order)

orderlist <- df_chain$entryName




df_DNA_BIND <- df_DNA_BIND %>% arrange(match(entryName, orderlist))

df_Peri <- df_Peri %>% arrange(match(entryName, orderlist))

df_chain$entryName <- factor(df_chain$entryName,levels = orderlist)
df_DNA_BIND$entryName <- factor(df_DNA_BIND$entryName,levels = orderlist)
df_Peri$entryName <- factor(df_Peri$entryName,levels = orderlist)


df_TAD <- read_csv("TAD_database_f.csv")

df_TAD$begin <- df_TAD$begin * -1
df_TAD$end <- df_TAD$end * -1


df_TAD <- df_TAD %>% arrange(match(entryName, orderlist))

df_TAD$entryName <- factor(df_TAD$entryName,levels = orderlist)





{r fig.height=4, fig.width=4}
p <- ggplot(df_chain, aes(y = entryName)) +
  geom_rect(data = df_chain, aes(xmin = begin, xmax = end, ymin=entryName, ymax=entryName), color ="#F6FFF6", alpha=1, size = 3) +
  geom_rect(data = df_DNA_BIND , aes( xmin = begin, xmax = end, ymax=entryName, ymin=entryName), colour = "#f6921e", alpha=1, size = 3) +
  geom_rect(data = df_Peri, aes(xmin = begin, xmax = end, ymin=entryName, ymax=entryName), colour ="#be1e2d", alpha=1, size = 3) +
  geom_rect(data = df_TAD, aes(xmin = begin, xmax = end, ymin=entryName, ymax=entryName), colour ="#008d36", alpha=1, size = 3) +
  #geom_point(data = AromaticContentDF %>% dplyr::filter(Aromatic == "Aromatic") ,aes(x = Position, y = GeneName) ,shape = 21, fill = "#f29121", colour = "#000000", size  = 3, stroke = 0.5) +
  #geom_point(data = AromaticContentDF %>% dplyr::filter(Periodic == "Periodic") ,aes(x = Position, y = GeneName) ,shape = 21, fill = "#be1e2d", colour = "#000000", size  = 3, stroke = 0.5) +
  
  theme_classic() +
  theme(axis.text=element_text(size=6)) + scale_y_discrete(position = "right")# + scale_x_continuous(trans = "reverse")

p 
ggsave("Periodic_vs_TAD_TOP80genes.pdf")

{r fig.height=4, fig.width=4}
p <- ggplot(df_chain, aes(y = entryName)) +
  geom_rect(data = df_chain, aes(xmin = begin, xmax = end, ymin=entryName, ymax=entryName), color ="#F6FFF6", alpha=1, size = 3) +
  #geom_rect(data = df_DNA_BIND , aes( xmin = begin, xmax = end, ymax=entryName, ymin=entryName), colour = "#f6921e", alpha=1, size = 3) +
  #geom_rect(data = df_Peri, aes(xmin = begin, xmax = end, ymin=entryName, ymax=entryName), colour ="#be1e2d", alpha=1, size = 3) +
  #geom_rect(data = df_TAD, aes(xmin = begin, xmax = end, ymin=entryName, ymax=entryName), colour ="#008d36", alpha=1, size = 3) +
  #geom_point(data = AromaticContentDF %>% dplyr::filter(Aromatic == "Aromatic") ,aes(x = Position, y = GeneName) ,shape = 21, fill = "#f29121", colour = "#000000", size  = 3, stroke = 0.5) +
  #geom_point(data = AromaticContentDF %>% dplyr::filter(Periodic == "Periodic") ,aes(x = Position, y = GeneName) ,shape = 21, fill = "#be1e2d", colour = "#000000", size  = 3, stroke = 0.5) +
  
  theme_classic() +
  theme(axis.text=element_text(size=6)) + scale_y_discrete(position = "right")# + scale_x_continuous(trans = "reverse")

p 
#ggsave("Periodic_vs_TAD_TOP80genes.pdf")

