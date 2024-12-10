setwd("~/suny/RESEARCH_DOC/ortholog_eukarya")
#dir.create("src")
dir.create("data")
dir.create("plots")

library(tidyverse)
#library(VennDiagram)

df_intra_process <- function(data){
    ## Taking the data frame it will split and remove parenthesis and
    ## convert the column type to numeric type.

    require(tidyverse)

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
        ) %>% select("AA1", "AA2", "Correlation", "Bootstrap")

    data <- data_splitted
    return(data)
}

## Pair of residues frequency calculation for each position of PCNA residue.
## it takes the processed splitted data frame and the protein name whos data is,
## and after filtering write a filtered dataframe to disk >0.8 correlation as
## .csv. Finally create a frequency table to use it in next satage.


prot_intra_cov_freq <- function(splitted_df, protein="") {
    require(dplyr)

    filtered_df <- splitted_df %>% filter(Correlation>=0.8)

    ## writing filtered coevolved residue pairs to a csv file
    folder="data/"
    filename=paste0(protein,"_intra.csv")
    write.csv(filtered_df, paste(folder, filename), row.names=FALSE)

    ## writing filtered coevolved network residue pairs to a csv file
    network_file=paste0(protein, "_intra_network.csv")
    filtered_df %>% select(AA1, AA2, Correlation) %>%
        write.csv(., paste(folder, network_file), row.names=FALSE)

    residue_list1= filtered_df %>% select(AA1)
    residue_list2= filtered_df %>% select(AA2)
    names(residue_list2) <- "AA1"
    total_residue= rbind(residue_list1, residue_list2)
    total_uniq_res <- unique(total_residue)
    # total_uniq is utilise in 4cosmic_match.R
    file=paste0("data/total_uniq_", protein)
    write.csv(total_uniq_res, file=file)

    filtered_df_freq <- filtered_df %>% group_by(AA1) %>%
        summarise(freq=n()) %>% mutate(protein=protein)

    return(filtered_df_freq)

}

## Barplot plotting function with ggplot


intra_draw_bar_plot <- function(frequncy_df, protein_name){
    require(ggplot2)
    p <- ggplot(frequncy_df, aes(AA1, y=freq, fill=protein)) +
        geom_bar(stat="identity")
    p <- p + theme_bw()

    intra_title= paste0(protein_name,  " intra co-evolving residues ( Correlation >0.8 )")
    p <- p + labs(title=intra_title,
                  fill="Protein") +
        scale_x_continuous(name=paste("N to C terminal residue position of", protein_name),
                           breaks=seq(0,1400,50)) +
        scale_y_continuous(name="Pairs of coevolving residues",
                           breaks=seq(0,100,10))
    p

}


