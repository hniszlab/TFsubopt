library(tidyverse)

df <- read_csv("grid_profile_2ablocks_Long.csv")

ggplot(df, aes(fill=Type, y=AA_amounts, x=Position)) + 
  geom_bar(position="fill", stat="identity")
ggsave("Disorder_freq.pdf")

ggplot(df[df$AA_amounts > 100,], aes(fill=AA, y=AA_amounts, x=Position)) + 
  geom_bar(position="fill", stat="identity")
ggsave("Leter_freq.pdf")