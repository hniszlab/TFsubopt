library(drawProteins)
library(Biostrings)
library(tidyverse)



AromaticContent <- function (sequence,  Aromatic = c("F", "Y", "W"), Other = c("P", "E", "S", "Q", 
                                                                               "K", "A", "G", "D", "T", "R", "M", "N", "V", "H", "L", "I", "C")) 
{
  datalist = list()
  for (i in 1:length(sequence)){
    seqCharacter <- sequence[[i]][["sequence"]]
    seqCharacter <- as.character(seqCharacter)
    seqCharacter <- strsplit(seqCharacter, "")
    seqCharacterVector <- unlist(seqCharacter)
    sequenceLength <- length(seqCharacterVector)
    NameVector <- rep(NA, sequenceLength)
    AromaticVector <- rep(NA, sequenceLength)
    AromaticResidues <- seqCharacterVector %in% Aromatic
    AromaticVector[AromaticResidues] <- "Aromatic"
    OtherResidues <- seqCharacterVector %in% Other
    AromaticVector[OtherResidues] <- ""
    df <- data.frame(GeneName = sequence[[i]][["entryName"]], Position = seq_len(sequenceLength), 
                     AA = seqCharacterVector, Aromatic = AromaticVector)
    df$AA <- as.character(df$AA)
    df$Aromatic <- as.character(df$Aromatic)
    datalist[[i]] <- df
  }
  AromaticContentDF = do.call(rbind, datalist)
  return(AromaticContentDF)
}



IDR_df <- read_csv("/Users/magalhae/Desktop/patchiness/figure1IDR_lolipop.csv")




rel_json <- get_features("A6NCS4 A6NHT5 A6NJ46 O00358 O00409 O00482 O00570 O14503 O14627 O14867 O15119 O15350 O15353 O43186 O43248 O43365 O43435 O43638 O43711 O60248 O60381 O60393 O60479 O75603 O75626 O94916 O94993 O95076 O95096 O95231 O95343 O95644 O95936 O95947 P04150 P04198 P05549 P08151 P09016 P09017 P09067 P09630 P10070 P10071 P10242 P11161 P11308 P12524 P12980 P14316 P14651 P14653 P15172 P15336 P15407 P15822 P15863 P15884 P15923 P16989 P17481 P17482 P17483 P18146 P18850 P20719 P22670 P23759 P23760 P23769 P28356 P28358 P28698 P31249 P31260 P31267")

rel_json2 <- get_features("P31268 P31269 P31270 P31273 P31274 P31277 P31629 P32242 P32243 P35712 P35869 P36956 P39880 P41161 P41182 P41225 P43694 P43699 P47902 P48378 P48380 P48431 P48436 P49640 P49715 P50222 P50549 P51449 P52945 P52951 P53539 P56178 P56693 P57071 P57682 P58012 P78411 P81133 Q00056 Q00444 Q01538 Q01543 Q01860 Q01954 Q02447 Q02548 Q02930 Q03014 Q03052 Q03112 Q03828 Q04864 Q06416 Q06710 Q07687 Q10586 Q10587 Q12772 Q12800 Q12947 Q12948 Q12950 Q12951 Q12968 Q12986 Q13129 Q13207 Q13285 Q13351 Q13562 Q13887 Q14494 Q14526 Q14865 Q14934 Q14938 Q15270 Q15562 Q15784 Q15788") 

rel_json3 <- get_features("Q15797 Q15853 Q16236 Q16534 Q16650 Q16665 Q33E94 Q3C1V8 Q3SYB3 Q4G112 Q5JT82 Q5T1R4 Q5VV16 Q5VZB9 Q68CJ9 Q6VB84 Q70SY1 Q7Z353 Q7Z6R9 Q86V15 Q86YP4 Q8HWS3 Q8IUM7 Q8IXT2 Q8N693 Q8NEA6 Q8NFW5 Q8WXT5 Q8WY36 Q92481 Q92570 Q92754 Q92786 Q92949 Q96AV8 Q96BA8 Q96IS3 Q96JB3 Q96NK8 Q96NZ1 Q96PN7 Q96QS3 Q96RK0 Q96S65 Q99081 Q99459")

rel_json4 <- get_features("Q99593 Q99607 Q99684 Q99697 Q99742 Q99743 Q99853 Q99958 Q9BQW3 Q9BW11 Q9BWX5 Q9BYV9 Q9BZI1 Q9BZM3 Q9BZS1 Q9C009 Q9C0J9 Q9C0K0 Q9GZU2 Q9GZV8 Q9H161 Q9H2P0 Q9H2S9 Q9H3D4 Q9H4Q3 Q9H4W6 Q9H6I2 Q9H9S0 Q9HAK2 Q9HAZ2 Q9HBE1 Q9HBZ2 Q9HCS4 Q9HD90 Q9NP08 Q9NP62 Q9NP71 Q9NQ03 Q9NQB0 Q9NQL9 Q9NQV7 Q9NQX0 Q9NU39 Q9NYD6 Q9P0K8 Q9UBR4 Q9UHF7 Q9UJU2 Q9UKT9 Q9UL17 Q9UPW0 Q9Y2G1 Q9Y2N7 Q9Y4A8 Q9Y5R6 Q9Y5W3")





rel_json <- append(rel_json, rel_json2)
rel_json <- append(rel_json, rel_json3)
rel_json <- append(rel_json, rel_json4)

prot_data <- feature_to_dataframe(rel_json)

AromaticContentDF <- AromaticContent(rel_json)


write.csv(AromaticContentDF, "/Users/magalhae/Desktop/Paper_suboptimitation/AverageComposition/Merge_chargedVSPB/AromaticContentDF.csv")
write.csv(prot_data, "/Users/magalhae/Desktop/Paper_suboptimitation/AverageComposition/Merge_chargedVSPB/prot_data.csv")




df_TAD <- read_csv("/Users/magalhae/Desktop/Paper_suboptimitation/SubOpt_revision/patchiness/TAD_database.csv")
df_PB <- read.csv("/Users/magalhae/Desktop/Paper_suboptimitation/SubOpt_revision/patchiness/Coord_block.csv")








rel_json <- get_features("P48436 Q92858 Q99814 Q9H0D2 O00321")
prot_data <- feature_to_dataframe(rel_json)
#rel_json[[1]][["sequence"]] <- "MDYNRMNSLEYPLCNRGPYSASAHSAFPTSPPSSFAQAVDSAYSEGRGGGYLSSPAQQFNSGPAQQYPPSTLGVFPPSSAPSYGAPAACSYPSGPSQPYLGQSEGDYGGHPSSGFAQLGGLSYDGGAGGAYGPGPPPQYHPPGNEQYTASAPAAFDLLSEDKYETPCPSEYPNTPTARFTDMKVKRWNPPKTAKVSEPGL"
AromaticContentDF <- AromaticContent(rel_json)
write.csv(AromaticContentDF, "/Users/magalhae/Desktop/Julian/AromaticContentDF.csv")
write.csv(prot_data, "/Users/magalhae/Desktop/Julian/prot_data.csv")
#prot_data <- read_csv("/Users/magalhae/Desktop/Julian/prot_data.csv")


AromaticContentDF <- read_csv("/Users/magalhae/Desktop/patchiness/AromaticContentDF.csv")


prot_data <- read_csv("/Users/magalhae/Desktop/Julian/prot_data.csv")



AromaticContentIDR <- function (sequence,  Aromatic = c("F", "Y", "W"), Other = c("P", "E", "S", "Q", 
                                                                                  "K", "A", "G", "D", "T", "R", "M", "N", "V", "H", "L", "I", "C")) 
{
  datalist = list()
  for (i in 1:length(sequence)){
    beggin
    end
    seqCharacter <- sequence[[i]][["sequence"]]
    seqCharacter <- as.character(seqCharacter)
    seqCharacter <- str_sub(seqCharacter, 2, -2)
    seqCharacter <- strsplit(seqCharacter, "")
    seqCharacterVector <- unlist(seqCharacter)
    sequenceLength <- length(seqCharacterVector)
    NameVector <- rep(NA, sequenceLength)
    AromaticVector <- rep(NA, sequenceLength)
    AromaticResidues <- seqCharacterVector %in% Aromatic
    AromaticVector[AromaticResidues] <- "Aromatic"
    OtherResidues <- seqCharacterVector %in% Other
    AromaticVector[OtherResidues] <- ""
    df <- data.frame(GeneName = sequence[[i]][["entryName"]], Position = seq_len(sequenceLength), 
                     AA = seqCharacterVector, Aromatic = AromaticVector)
    df$AA <- as.character(df$AA)
    df$Aromatic <- as.character(df$Aromatic)
    datalist[[i]] <- df
  }
  AromaticContentDF = do.call(rbind, datalist)
  return(AromaticContentDF)
}



#AromaticContentDF$GeneName <-  factor(AromaticContentDF$GeneName, levels = c("HXD4_HUMAN", "HXB1_HUMAN", "TBX5_HUMAN", "NANOG_HUMAN", "IRX5_HUMAN", "HXB3_HUMAN", "EGR1_HUMAN", #"COE1_HUMAN", "H0Y3W9_HUMAN", "OTX1_HUMAN"))
#prot_data$entryName <- factor(prot_data$entryName, levels = c("HXD4_HUMAN", "HXB1_HUMAN", "TBX5_HUMAN", "NANOG_HUMAN", "IRX5_HUMAN", "HXB3_HUMAN", "EGR1_HUMAN", "COE1_HUMAN", #"H0Y3W9_HUMAN", "OTX1_HUMAN"))

p <- ggplot() +
  geom_rect(data = prot_data %>% dplyr::filter(type == "DNA_BIND"), aes( xmin = begin, xmax = end, ymax=entryName, ymin=entryName), colour = "#f29121", alpha=0.2, linewidth = 3) +
  geom_rect(data = prot_data %>% dplyr::filter(description == "IDR") ,aes(xmin = begin, xmax = end, ymax = entryName, ymin = entryName), colour = "#28316f", alpha=0.2, linewidth = 3) +
  geom_rect(data = prot_data %>% dplyr::filter(description == "TAD") ,aes(xmin = begin, xmax = end, ymax = entryName, ymin = entryName), colour = "#039be5", alpha=0.2, linewidth = 3) +
  geom_rect(data = prot_data %>% dplyr::filter(description == "PBLOCK") ,aes(xmin = begin, xmax = end, ymax = entryName, ymin = entryName), colour = "#be1e2d", alpha=0.2, linewidth = 3) +
  geom_segment(data = prot_data %>% dplyr::filter(type == "CHAIN") ,aes(x = begin, y = entryName, xend = end, yend = entryName), colour = "#000000") +
  geom_point(data = AromaticContentDF %>% dplyr::filter(Aromatic == "Aromatic") ,aes(x = Position, y = GeneName) ,shape = 21, fill = "#f29121", colour = "#000000", size  = 3, stroke = 0.5) +
  #geom_point(data = AromaticContentDF %>% dplyr::filter(Periodic == "Periodic") ,aes(x = Position, y = GeneName) ,shape = 21, fill = "#be1e2d", colour = "#000000", size  = 3, stroke = 0.5) +
  
  theme_classic()

p 


ggsave( "Top100_periodic.pdf",height=15, width=12, units = "in" )


