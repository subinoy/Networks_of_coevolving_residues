setwd("~/suny/RESEARCH_DOC/A_THESIS_WORK/CHAPTER_5_COEVOLV_RESIDUE/eukarya")

## system('pandoc  -t latex --standalone --smart --number-sections --template=report.tex -f markdown -o whitepaper.pdf "whitepaper.Rmd"')

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)


fen1_pcna <- read_csv("FEN1_INTER_COEV.csv")
fen1_pcna
fen1_splitted <- fen1_pcna %>% 
    separate(., AA1, into=c("AA1_algn", "AA1"), sep="\\(", extra="drop")

fen1_splitted <- fen1_splitted %>% 
    separate(., AA1, into=c("AA1"), sep="\\)", extra="drop")


fen1_splitted <- fen1_splitted %>% 
    separate(., AA2, into=c("AA2_algn", "AA2"), sep="\\(", extra="drop")

fen1_splitted <- fen1_splitted %>% 
    separate(., AA2, into=c("AA2"), sep="\\)", extra="drop")
fen1_splitted
##Sys.getenv("PATH")
