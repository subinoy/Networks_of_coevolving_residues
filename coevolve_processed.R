setwd("~/suny/RESEARCH_DOC/A_THESIS_WORK/CHAPTER_5_COEVOLV_RESIDUE/eukarya")

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)


fen1_pcna <- read_csv("FEN1_INTER_COEV.csv")
fen1_pcna
fen1_splitted <- fen1_pcna %>% separate(AA1, "AA1_algn1", "AA1", sep="\\(")
