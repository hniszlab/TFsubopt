library(universalmotif)


data <- read_meme("glam2_chargedBlock.meme", readsites = TRUE, readsites.meta = TRUE)
 
view_motifs(data$motifs[[1]])
view_motifs(data$motifs[[2]])
view_motifs(data$motifs[[3]])


data2 <- read_meme("glam2_PD.meme", readsites = TRUE, readsites.meta = TRUE)

view_motifs(data2$motifs[[1]])
view_motifs(data2$motifs[[2]])
view_motifs(data2$motifs[[3]])