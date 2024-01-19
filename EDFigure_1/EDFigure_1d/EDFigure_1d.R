```{r}
library(Biostrings)
library(tidyverse)


HG38_pep_ref <- readAAStringSet("HG38_pep_ref.fasta")

TF_PB_4aablock_sequences <- readAAStringSet("TF_PB_4aablock_sequences.fasta")
TF_PB_8aablock_sequences <- readAAStringSet("TF_PB_8aablock_sequences.fasta")
TF_PB_sequences <- readAAStringSet("TF_PB_sequences.fasta")

Charged_blocks_sabari_top500 <- readAAStringSet( "Charged_blocks_sabari_top500.fasta")
Charged_blocks_sabari <- readAAStringSet( "Charged_blocks_sabari.fasta")


HG38_pep_ref_freq
HG38_pep_ref_IDRs_freq
HumanTF_canonical_freq
TF_HG38_pep_IDR_freq
HG38_IDRs_NONperiodic_freq
HG38_IDRs_periodic_freq


HG38_pep_ref_freq <- HG38_pep_ref_freq
HG38_pep_ref_IDRs_freq <- HG38_pep_ref_IDRs_freq
HumanTF_canonical_freq <- HumanTF_canonical_freq
TF_HG38_pep_IDR_freq <- TF_HG38_pep_IDR_freq
HG38_IDRs_NONperiodic_freq <- HG38_IDRs_NONperiodic_freq
HG38_IDRs_periodic_freq <- HG38_IDRs_periodic_freq

, collapse=TRUE, 


amino_letters <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "P", "S", "T", "V", "F", "Y", "W")
#amino_letters <-  "RKEDQNGASTPHYCLVIMFW"

HG38_pep_ref_P <- as.matrix(alphabetFrequency(HG38_pep_ref, collapse=TRUE, as.prob = T))

TF_PB_4aablock_sequences_P <- as.matrix(alphabetFrequency(TF_PB_4aablock_sequences, collapse=TRUE, as.prob = T))
TF_PB_8aablock_sequences_P <- as.matrix(alphabetFrequency(TF_PB_8aablock_sequences, collapse=TRUE, as.prob = T))
TF_PB_sequences_P <- as.matrix(alphabetFrequency(TF_PB_sequences, collapse=TRUE, as.prob = T))

Charged_blocks_sabari_top500_P <- as.matrix(alphabetFrequency( Charged_blocks_sabari_top500, collapse=TRUE, as.prob = T))
Charged_blocks_sabari_P <- as.matrix(alphabetFrequency( Charged_blocks_sabari, collapse=TRUE, as.prob = T))


HG38_pep_ref_freq <- as.matrix(alphabetFrequency(HG38_pep_ref, collapse=TRUE, as.prob = F))
TF_PB_4aablock_sequences_freq <- as.matrix(alphabetFrequency(TF_PB_4aablock_sequences, collapse=TRUE, as.prob = F))
TF_PB_8aablock_sequences_freq <- as.matrix(alphabetFrequency(TF_PB_8aablock_sequences, collapse=TRUE, as.prob = F))
TF_PB_sequences_freq <- as.matrix(alphabetFrequency(TF_PB_sequences, collapse=TRUE, as.prob = F))

Charged_blocks_sabari_top500_freq <- as.matrix(alphabetFrequency( Charged_blocks_sabari_top500, collapse=TRUE, as.prob = F))
Charged_blocks_sabari_freq <- as.matrix(alphabetFrequency( Charged_blocks_sabari, collapse=TRUE, as.prob = F))

BIG_p <- as.data.frame(cbind(HG38_pep_ref_P, TF_PB_4aablock_sequences_P, TF_PB_8aablock_sequences_P, TF_PB_sequences_P, Charged_blocks_sabari_top500_P, Charged_blocks_sabari_P))
BIG_f <- as.data.frame(cbind(HG38_pep_ref_freq, TF_PB_4aablock_sequences_freq, TF_PB_8aablock_sequences_freq, TF_PB_sequences_freq, Charged_blocks_sabari_top500_freq, Charged_blocks_sabari_freq))


colnames(BIG_p) <- c("HG38_pep_ref",	"TF_PB_4aablock",	"TF_PB_8aablock",	"TF_PB", "Charged_blocks_sabari_top500", "Charged_blocks_sabari")
colnames(BIG_f) <- c("HG38_pep_ref",	"TF_PB_4aablock",	"TF_PB_8aablock",	"TF_PB", "Charged_blocks_sabari_top500", "Charged_blocks_sabari")




write.csv(BIG_p, "BIG_p.csv") #This file creates AA_pro_normHG38
write.csv(BIG_f, "BIG_F.csv")
 
AA_pro_normHG38 <- read_csv("AA_pro_normHG38.csv")
 
Amino	Freq	Sample

 
amino_letters <- c("P", "M", "C", "A", "G", "I", "V", "L", "N", "Q", "R", "K", "D", "E", "S", "T", "H", "F", "W", "Y")

AA_pro_normHG38 <- AA_pro_normHG38 %>% arrange(factor(Amino, levels = amino_letters))


 
AA_pro_normHG38 %>% arrange(factor(Amino, levels = amino_letters)) %>% ggplot(
  aes(x=Amino, y=Freq)) +
  geom_bar(stat="identity") +  facet_grid(rows = vars(Sample))+ theme_classic() 


ggsave("rplot.pdf")
 