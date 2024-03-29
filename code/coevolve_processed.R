#setwd("~/suny/RESEARCH_DOC/A_THESIS_WORK/CHAPTER_5_COEVOLV_RESIDUE/eukarya")

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(VennDiagram)

## Testing ggplot sequencially
# p <- ggplot(fen1_freq, aes(AA1, y=freq, fill=protein)) + geom_bar(stat="identity") 
# p <- p + theme_bw()
# p <- p + scale_x_continuous(name="N to C terminal residue position of PCNA ",
#                             breaks=seq(0,263,15)) +
#     scale_y_continuous(name="Number of coevolving residues pair", limits=c(0, 30))
# p+ labs(title="PCNA inter co-evolving residues >0.8 Correlation", fill="Protein")


## Function for processing the raw coevolving residue pairs of inter coevolution from CAPS server.

data_frame_process <- function(data){
## Taking the data frame it will split and remove parenthesis and 
## convert the column type to numeric type.
    

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

    # Charater to numeric type conversion
    data_splitted <- data_splitted %>% 
    mutate(AA1_algn=as.numeric(AA1_algn), AA1=as.numeric(AA1),
           AA2_algn=as.numeric(AA2_algn), AA2=as.numeric(AA2)
            )

    data <- data_splitted
    return(data)
}

## Pair of residues frequency calculation for each position of PCNA residue.
## it takes the processed splitted data frame and the protein name whos data is,
## and after filtering write a filtered dataframe to disk >0.8 correlation as 
## .csv. Finally create a frequency table to use it in next satage.


protein_coevolve_frequency <- function(splitted_df, protein="") {
    require(dplyr)
    
    filtered_df <- splitted_df %>% 
    filter(Correlation>=0.8) 
    
    ## writing filtered coevolved residue pairs to a csv file
    folder="data/"
    filename=paste0(protein,"_inter.csv")
    write.csv(filtered_df, paste(folder, filename), row.names=FALSE)
    
    ## writing filtered coevolved network residue pairs to a csv file
    network_file=paste0(protein, "_inter_network.csv")
    filtered_df %>% select(AA1, AA2, Correlation) %>% 
        write.csv(., paste(folder, network_file), row.names=FALSE)
    
    filtered_df_freq <- filtered_df %>% group_by(AA1) %>% 
    summarise(freq=n()) %>% mutate(protein=protein)
    
    return(filtered_df_freq)

}

## Barplot plotting function with ggplot


draw_bar_plot <- function(frequncy_df){
    require(ggplot2)
    p <- ggplot(frequncy_df, aes(AA1, y=freq, fill=protein)) + 
        geom_bar(stat="identity") 
    p <- p + theme_bw()
    
    p <- p + labs(title="PCNA inter co-evolving residues ( Correlation >0.8 )",
                  fill="Protein") +
        scale_x_continuous(name="N to C terminal residue position of PCNA",
                           breaks=seq(0,263,15)) +
        scale_y_continuous(name="Pairs of coevolving residues",
                           breaks=seq(0,30,5))
    p
    
}

## FEN1 
fen1_pcna <- read_csv("data/FEN1_INTER_COEV.csv")
fen1_pcna
## Applying function to process the data frame
fen1 <- data_frame_process(fen1_pcna)
## Generating frequency table
fen1_co_freq <- protein_coevolve_frequency(fen1,"Fen1")
## Plotting bar plot
draw_bar_plot(fen1_co_freq)
## saving to disk
ggsave("fen1/fen1_inter_plot.tiff", width=9, height=6, units= "in", dpi = 300 )
ggsave("plots/fen1_inter_plot.png", width=9, height=6, units= "in", dpi = 300 )

## POLD

pold_pcna <- read_csv("data/POLD_INTER_COEV.csv")
pold_pcna
## Applying function to process the data frame
pold <- data_frame_process(pold_pcna)
## Generating frequency table
pold_co_freq <- protein_coevolve_frequency(pold,"Pold")
## Plotting bar plot
draw_bar_plot(pold_co_freq)


## saving to disk
ggsave("pold/pold_inter_plot.tiff", width=9, height=6, units= "in", dpi = 300 )
ggsave("plots/pold_inter_plot.png", width=9, height=6, units= "in", dpi = 300 )


## RFC3

rfc3_pcna <- read_csv("data/RFC3_INTER_COEV.csv")
rfc3_pcna
## Applying function to process the data frame
rfc3 <- data_frame_process(rfc3_pcna)
## Generating frequency table
rfc3_co_freq <- protein_coevolve_frequency(rfc3,"RFC3")
## Plotting bar plot
draw_bar_plot(rf3_co_freq)
## saving to disk
ggsave("rfc3/rfc3_inter_plot.tiff", width=9, height=6, units= "in", dpi = 300 )
ggsave("plots/rfc3_inter_plot.png", width=9, height=6, units= "in", dpi = 300 )

## PCNA INTRA ********
# ***********************************************
pcna_intra <- read_csv("data/PCNA_INTRA_COEV.csv")
pcna_intra
## Applying function to process the data frame
pcna <- data_frame_process(pcna_intra)
## Generating frequency table
pcna_co_freq <- protein_coevolve_frequency(pcna,"PCNA")
## Plotting bar plot
draw_bar_plot(pcna_co_freq)

p <- ggplot(pcna_co_freq, aes(AA1, y=freq, fill=protein)) + geom_bar(stat="identity")
#p <- p + facet_wrap(~protein, ncol=1)
p <- p + theme_bw()

p <- p + labs(title="PCNA INTRA co-evolving residues ( Correlation >0.8 )",
              fill="Protein") +
    scale_x_continuous(name="N to C terminal residue position of PCNA",
                       breaks=seq(0,263,15)) +
    scale_y_continuous(name="Pairs of coevolving residues",
                       breaks=seq(0,30,5))
p

## saving to disk
ggsave("pcna/pcna_intra_plot.tiff", width=9, height=6, units= "in", dpi = 300 )
ggsave("plots/pcna_intra_plot.png", width=9, height=6, units= "in", dpi = 300 )

# *********************************************

# WORKING ON Combined data frame

total_df=bind_rows(fen1_co_freq, pold_co_freq, rf3_co_freq )

## facet wrap
p <- ggplot(total_df, aes(AA1, y=freq, fill=protein)) + geom_bar(stat="identity")
p <- p + facet_wrap(~protein, ncol=1)
p <- p + theme_bw()

p <- p + labs(title="PCNA inter co-evolving residues ( Correlation >0.8 )",
              fill="Protein") +
    scale_x_continuous(name="N to C terminal residue position of PCNA",
                       breaks=seq(0,263,15)) +
    scale_y_continuous(name="Pairs of coevolving residues",
                       breaks=seq(0,30,5))
p

ggsave("plots/combined_inter_plot.png", width=9, height=6, units= "in", dpi = 300 )

# ************* PCNA and 3 interacting partners together

PCNA_n_INTER_df=bind_rows(fen1_co_freq, pold_co_freq, rf3_co_freq, pcna_co_freq )
## Reordering plots
#temp$size_f = factor(temp$size, levels=c('50%','100%','150%','200%'))
PCNA_n_INTER_df$protein_f = factor(PCNA_n_INTER_df$protein, levels=c('PCNA','Fen1','Pold','RFC3'))
p <- ggplot(PCNA_n_INTER_df, aes(AA1, y=freq, fill=protein_f)) + geom_bar(stat="identity")
p <- p + facet_wrap(~protein_f, ncol=1)
p <- p + theme_bw()

p <- p + labs(title="Co-evolving residues of PCNA \"intra\" and with three interacting partners  ( Correlation >0.8 )",
              fill="Protein") +
    scale_x_continuous(name="N to C terminal residue position of PCNA",
                       breaks=seq(0,263,15)) +
    scale_y_continuous(name="Pairs of coevolving residues",
                       breaks=seq(0,30,10))
p

ggsave("plots/combined_inter_n_PCNA_INTRA_plot.png", width=9, height=6, units= "in", dpi = 300 )

# *******************************

# mylist <- list(Vector1=pcna_co_freq$AA1, Vector2=fen1_co_freq$AA1)
# venn.diagram(mylist, filename="plots/venn_plot.tiff")


venn.plot <- venn.diagram(
    x = list(
        PCNA = pcna_co_freq$AA1,
        Fen1 = fen1_co_freq$AA1,
        Pold = pold_co_freq$AA1,
        RFC3 = rfc3_co_freq$AA1
    ),
    filename = "plots/Venn_4set_PCNA_pretty.png",
    col = "transparent",
    fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
    alpha = 0.50,
    label.col = c("orange", "white", "darkorchid4", "white", 
                  "white", "white", "white", "white", "darkblue", "white", 
                  "white", "white", "white", "darkgreen", "white"),
    cex = 1.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
    cat.cex = 1.5,
    cat.pos = 0,
    cat.dist = 0.07,
    cat.fontfamily = "serif",
    rotation.degree = 0,
    margin = 0.2
);

# **************
## fucntion for creating inner join data frame
creating_inner_df <- function(co_freq_df){
    inner_df <- as.data.frame(co_freq_df$AA1)
    names(inner_df) <- "AA1"
    return(inner_df)
}

fen1_df <- as.data.frame(fen1_co_freq$AA1)
names(fen1_df) <- "AA1"
fen1_df

pcna_df <- as.data.frame(pcna_co_freq$AA1)
names(pcna_df) <- "AA1"
pcna_df
pcna_fen1_inner <- inner_join(pcna_df, fen1_df, by="AA1")

pold_df <- creating_inner_df(pold_co_freq)

pcna_pold_inner <- inner_join(pcna_df, pold_df, by="AA1")
pcna_pold_inner
pold_fen1_inner <- inner_join(pold_df, fen1_df, by="AA1")
pold_fen1_inner
rfc3_df <- creating_inner_df(rfc3_co_freq)
pold_rfc3_inner <- inner_join(pold_df, rfc3_df, by="AA1")
pold_rfc3_inner
