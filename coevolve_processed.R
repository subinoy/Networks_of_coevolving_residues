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
fen1_splitted <- fen1_splitted %>% 
    mutate(AA1_algn=as.numeric(AA1_algn), AA1=as.numeric(AA1),
           AA2_algn=as.numeric(AA2_algn), AA2=as.numeric(AA2)
           )
fen1_freq <- fen1_splitted %>% 
    filter(Correlation>=0.8) %>% group_by(AA1) %>% 
    summarise(freq=n()) %>% mutate(protein="Fen1")
fen1_freq

p <- ggplot(fen1_freq, aes(AA1, y=freq, fill=protein)) + geom_bar(stat="identity")
p+ labs(fill="Protein")
##Sys.getenv("PATH")

data_frame_process <- function(data){
    
require(dplyr)
require(tidyr)

data_splitted <- data %>% 
    separate(., AA1, into=c("AA1_algn", "AA1"), sep="\\(", extra="drop")

data_splitted <- data_splitted %>% 
    separate(., AA1, into=c("AA1"), sep="\\)", extra="drop")


data_splitted <- data_splitted %>% 
    separate(., AA2, into=c("AA2_algn", "AA2"), sep="\\(", extra="drop")

data_splitted <- data_splitted %>% 
    separate(., AA2, into=c("AA2"), sep="\\)", extra="drop")

data_splitted <- data_splitted %>% 
    mutate(AA1_algn=as.numeric(AA1_algn), AA1=as.numeric(AA1),
           AA2_algn=as.numeric(AA2_algn), AA2=as.numeric(AA2)
    )

return(data_splitted)
}

protein_coevolve_frequency <- function(splitted_df, protein="") {
    require(dplyr)
    
    freq_df <- splitted_df %>% 
    filter(Correlation>=0.8) %>% group_by(AA1) %>% 
    summarise(freq=n()) %>% mutate(protein="protein")
    
    return(freq_df)

}

fen1 <- data_frame_process(fen1_pcna)
fen1_co_freq <- protein_coevolve_frequency(fen1,fen1_prot)

draw_bar_plot <- function(frequncy_df){
    require(ggplot2)
    p <- ggplot(frequncy_df, aes(AA1, y=freq, fill=protein)) + geom_bar(stat="identity")
    p+ labs(fill="Protein")
}

draw_bar_plot(fen1_co_freq)
