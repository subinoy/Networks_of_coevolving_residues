
## ****************************************************************************
## Upto complete network topology plot
## Then move to 8Network_function.R and copy bigclust
## code %>%  move to 11Network_process.R 
## ****************************************************************************

setwd("~/suny/RESEARCH_DOC/ortholog_eukarya")
library(igraph)
library(dplyr)
library(igraphdata)
library(ggplot2)
library(RColorBrewer)
#### ************************
#   PCNA-FEN1
# ***************************

fen1_nnet <- read.csv("data/ Fen1_inter_network.csv")

# without filter
fen1_net <-fen1_nnet %>% group_by(Col1) %>% 
    summarise(freq=n()) 

fen1_node_net <- fen1_net %>% mutate(protein="pcna")
fen1_node_net$Col1 <- as.character(fen1_node_net$Col1)

# data frame with all degree node
fen1_clean_net <- left_join(fen1_net, fen1_nnet)
fen1_clean_net <- fen1_clean_net %>% mutate(t2_node=paste0("F", Col2))
fen1_clean_net <- fen1_clean_net %>% select(Col1, t2_node, Correlation)

## Getting the degree value of target Fen1 protein 
fen1_target_node <-fen1_clean_net %>% group_by(t2_node) %>% 
    summarise(freq=n())

fen1_target_node$protein <- "fen1"
# renaming variable name in node file
names(fen1_target_node) <- c("Col1", "freq", "protein")

# Combining source and target node by row binding to total node list
fen1_tot_node <- bind_rows(fen1_node_net, fen1_target_node)

## Edgelist and network generation 
## *************************************************
# edgelist file without freq(degree) value
edgelist_fen1 <- fen1_clean_net

node_list_fen1 <- fen1_tot_node

# Creating graph from data frame and simplyfying
ig_fen1 = graph_from_data_frame(edgelist_fen1, directed = F, vertices = node_list_fen1)

ig_fen1 <- simplify(ig_fen1, remove.multiple = TRUE, remove.loops = TRUE,
                    edge.attr.comb = igraph_opt("edge.attr.comb"))


### Show all the colour schemes available
#display.brewer.all()

# R has 657 built in color names. To see
# colors()

#plot(ig_fen1, vertex.color="blueviolet", vertex.label=NA)

fen1_label_node <- sort(betweenness(ig_fen1), decreasing = T)[1:10]

# Colouring vertex according to protein name ************
V(ig_fen1)$color= "orange"
# Colouring vertex according to protein name ************
V(ig_fen1)$color=ifelse(V(ig_fen1)$protein =="pcna",
                        "#CC0099", "cyan2")


fen1_full_ns <- ifelse(betweenness(ig_fen1) > 10.0, betweenness(ig_fen1)/(max(betweenness(ig_fen1) * .075)), 3)

plot(ig_fen1, layout=l_fen1a,
     vertex.label=NA,
     #vertex.size=betweenness(ig_fen1)/(max(betweenness(ig_fen1) * .045)),
     vertex.size = fen1_full_ns,
     edge.arrow.size=0, vertex.label.color="black")


## *********************************************
#### ************************
#   PCNA-POLD2
# ***************************
pold2_net <- read.csv("data/ Pold2_inter_network.csv")

## Getting the degree value of target Pold protein 
pold2_target_node_all <- pold2_net %>% group_by(Col2) %>% 
    summarise(freq=n())


#without filtering
 pold2_net_gt4 <-pold2_net %>% group_by(Col1) %>% 
     summarise(freq=n())

pold2_net_gt4 <- pold2_net_gt4 %>% mutate(protein="pcna")
#pold_net_gt4$AA1 <- as.character(pold_net_gt4$AA1)

pold2_gt4_clean_net <- left_join(pold2_net_gt4, pold2_net)
pold2_gt4_clean_net <- pold2_gt4_clean_net %>% mutate(t2_node=paste0("P2-", Col2))
pold2_gt4_clean_net <- pold2_gt4_clean_net %>% select(Col1, t2_node, Correlation)

## Getting the degree value of target Pold protein 
pold2_gt4_target_node <- pold2_gt4_clean_net %>% group_by(t2_node) %>% 
    summarise(freq=n())

pold2_gt4_target_node$protein <- "pold2"
# renaming variable name in node file
names(pold2_gt4_target_node) <- c("Col1", "freq", "protein")


pold2_gt4_node_net <- pold2_net_gt4
pold2_gt4_node_net$Col1 <- as.character(pold2_gt4_node_net$Col1)
# Combining source and target node by row binding to total node list
tot_node_gt4 <- bind_rows(pold2_gt4_node_net, pold2_gt4_target_node)


# edgelist file without freq(degree) value
edgelist_pold2 <- pold2_gt4_clean_net

node_list_pold2 <- tot_node_gt4

ig_pold2 = graph_from_data_frame(edgelist_pold2, directed = F, vertices = node_list_pold2)
ig_pold2 <- simplify(ig_pold2, remove.multiple = TRUE, remove.loops = TRUE,
                    edge.attr.comb = igraph_opt("edge.attr.comb"))


sort(betweenness(ig_pold2), decreasing = T)[1:10]

# Coloring the node based on name of protein node
V(ig_pold2)$color=ifelse(V(ig_pold2)$protein =="pcna",
                        "#CC0099", "skyblue2")


# ++++++++    POLD2 Full net node size definition   ++++++++++

# POLD2 Full net node size definition
pold2_full_ns <- ifelse(betweenness(ig_pold2) > 10.0, 
                        betweenness(ig_pold2)/
                            (max(betweenness(ig_pold2) * .075)), 3)



## *********************************************
#### ************************
#   PCNA-RFC3
# ***************************
rfc3_net_init <- read.csv("data/ RFC3_inter_network.csv")

# Removing 0 
rfc3_net <- rfc3_net_init %>% filter(Col2 != 0)
## Getting the degree value of target RFC3 protein 
rfc3_target_node_all <- rfc3_net %>% group_by(Col2) %>% 
    summarise(freq=n())
# greater than 4 degree vertices are filtered and nodes
# are grouped with unique node name
rfc3_net_gt4 <-rfc3_net %>% group_by(Col1) %>% 
    summarise(freq=n()) 

# without filtering
rfc3_net_gt4 <-rfc3_net %>% group_by(Col1) %>% 
    summarise(freq=n())


rfc3_net_gt4 <- rfc3_net_gt4 %>% mutate(protein="pcna")
rfc3_net_gt4$Col1 <- as.character(rfc3_net_gt4$Col1)
rfc3_net$Col1 <- as.character(rfc3_net$Col1)

rfc3_gt4_clean_net <- left_join(rfc3_net_gt4, rfc3_net)
rfc3_gt4_clean_net <- rfc3_gt4_clean_net %>% mutate(t2_node=paste0("R3-", Col2))
rfc3_gt4_clean_net <- rfc3_gt4_clean_net %>% select(Col1, t2_node, Correlation)

## Getting the degree value of target RFC3 protein 
rfc3_gt4_target_node <- rfc3_gt4_clean_net %>% group_by(t2_node) %>% 
    summarise(freq=n())

rfc3_gt4_target_node$protein <- "rfc3"
# renaming variable name in node file
names(rfc3_gt4_target_node) <- c("Col1", "freq", "protein")


rfc3_gt4_node_net <- rfc3_net_gt4
rfc3_gt4_node_net$Col1 <- as.character(rfc3_gt4_node_net$Col1)
# Combining source and target node by row binding to total node list
tot_node_gt4 <- bind_rows(rfc3_gt4_node_net, rfc3_gt4_target_node)


# edgelist file without freq(degree) value
edgelist_rfc3 <- rfc3_gt4_clean_net

node_list_rfc3 <- tot_node_gt4

# Creating graph from data frame and simplyfying
ig_rfc3 = graph_from_data_frame(edgelist_rfc3, directed = F, 
                                vertices = node_list_rfc3)
ig_rfc3 <- simplify(ig_rfc3, remove.multiple = TRUE, remove.loops = TRUE,
         edge.attr.comb = igraph_opt("edge.attr.comb"))

# Coloring the node based on name of protein node
V(ig_rfc3)$color=ifelse(V(ig_rfc3)$protein =="pcna",
                        "#CC0099", "#FF6600")


sort(betweenness(ig_rfc3), decreasing = T)[1:5]


# ++++++++    RFC3 Full net node size definition   ++++++++++
rfc3_full_ns <- ifelse(betweenness(ig_rfc3) > 10.0, 
                        betweenness(ig_rfc3)/
                            (max(betweenness(ig_rfc3) * .06)), 2.5)
                      
plot(ig_rfc3, layout=layout.fruchterman.reingold, 
     vertex.label=NA,
     vertex.label.cex= 0.5,
     vertex.label.dist= 0.9,
     #vertex.size=betweenness(ig_rfc3)/(max(betweenness(ig_rfc3) * .055)),
     vertex.size=rfc3_full_ns,
     edge.arrow.size=0, vertex.label.color="black")



#### ************************
#   PCNA-RFC1
# ***************************
network_draw <- function(inter_network_file, partner_protein){
    
    
    
    protein_net_init <- 
}

rfc1_net_init <- read.csv("data/ RFC1_inter_network.csv")

# Removing 0 
rfc1_net <- rfc1_net_init %>% filter(Col2 != 0)
## Getting the degree value of target RFC3 protein 
rfc1_target_node_all <- rfc1_net %>% group_by(Col2) %>% 
    summarise(freq=n())

# without filtering
rfc1_net_gt4 <-rfc1_net %>% group_by(Col1) %>% 
    summarise(freq=n())


rfc1_net_gt4 <- rfc1_net_gt4 %>% mutate(protein="pcna")
rfc1_net_gt4$Col1 <- as.character(rfc1_net_gt4$Col1)
rfc1_net$Col1 <- as.character(rfc1_net$Col1)

rfc1_gt4_clean_net <- left_join(rfc1_net_gt4, rfc1_net)
rfc1_gt4_clean_net <- rfc1_gt4_clean_net %>% mutate(t2_node=paste0("R1-", Col2))
rfc1_gt4_clean_net <- rfc1_gt4_clean_net %>% select(Col1, t2_node, Correlation)

## Getting the degree value of target RFC1 protein 
rfc1_gt4_target_node <- rfc1_gt4_clean_net %>% group_by(t2_node) %>% 
    summarise(freq=n())
rfc1_gt4_target_node$protein <- "rfc1"

# renaming variable name in node file
names(rfc1_gt4_target_node) <- c("Col1", "freq", "protein")
rfc1_gt4_node_net <- rfc1_net_gt4
rfc1_gt4_node_net$Col1 <- as.character(rfc1_gt4_node_net$Col1)

# Combining source and target node by row binding to total node list
tot_node_gt4 <- bind_rows(rfc1_gt4_node_net, rfc1_gt4_target_node)


# edgelist file without freq(degree) value
edgelist_rfc1 <- rfc1_gt4_clean_net

node_list_rfc1 <- tot_node_gt4

# Creating graph from data frame and simplyfying
ig_rfc1 = graph_from_data_frame(edgelist_rfc1, directed = F, 
                                vertices = node_list_rfc1)
ig_rfc1 <- simplify(ig_rfc1, remove.multiple = TRUE, remove.loops = TRUE,
                    edge.attr.comb = igraph_opt("edge.attr.comb"))

# Coloring the node based on name of protein node
V(ig_rfc1)$color=ifelse(V(ig_rfc1)$protein =="pcna",
                        "#CC0099", "#FF6600")

sort(betweenness(ig_rfc1), decreasing = T)[1:5]

# ++++++++    RFC1 Full net node size definition   ++++++++++
rfc1_full_ns <- ifelse(betweenness(ig_rfc1) > 10.0, 
                       betweenness(ig_rfc1)/
                           (max(betweenness(ig_rfc1) * .06)), 2.5)

plot(ig_rfc1, layout=layout.fruchterman.reingold, 
     vertex.label=NA,
     vertex.label.cex= 0.8,
     #vertex.label.dist= 0.9,
     vertex.size=rfc1_full_ns,
     edge.arrow.size=0 ,vertex.label.color="black")



#### ************************
#   PCNA-RFC4
# ***************************

rfc4_net_init <- read.csv("data/ RFC4_inter_network.csv")

# Removing 0 
rfc4_net <- rfc4_net_init %>% filter(Col2 != 0)
## Getting the degree value of target RFC3 protein 
rfc4_target_node_all <- rfc4_net %>% group_by(Col2) %>% 
    summarise(freq=n())

# without filtering
rfc4_net_gt4 <-rfc4_net %>% group_by(Col1) %>% 
    summarise(freq=n())


rfc4_net_gt4 <- rfc4_net_gt4 %>% mutate(protein="pcna")
rfc4_net_gt4$Col1 <- as.character(rfc4_net_gt4$Col1)
rfc4_net$Col1 <- as.character(rfc4_net$Col1)

rfc4_gt4_clean_net <- left_join(rfc4_net_gt4, rfc4_net)
rfc4_gt4_clean_net <- rfc4_gt4_clean_net %>% mutate(t2_node=paste0("R4-", Col2))
rfc4_gt4_clean_net <- rfc4_gt4_clean_net %>% select(Col1, t2_node, Correlation)

## Getting the degree value of target rfc4 protein 
rfc4_gt4_target_node <- rfc4_gt4_clean_net %>% group_by(t2_node) %>% 
    summarise(freq=n())
rfc4_gt4_target_node$protein <- "rfc4"

# renaming variable name in node file
names(rfc4_gt4_target_node) <- c("Col1", "freq", "protein")
rfc4_gt4_node_net <- rfc4_net_gt4
rfc4_gt4_node_net$Col1 <- as.character(rfc4_gt4_node_net$Col1)

# Combining source and target node by row binding to total node list
tot_node_gt4 <- bind_rows(rfc4_gt4_node_net, rfc4_gt4_target_node)


# edgelist file without freq(degree) value
edgelist_rfc4 <- rfc4_gt4_clean_net

node_list_rfc4 <- tot_node_gt4

# Creating graph from data frame and simplyfying
ig_rfc4 = graph_from_data_frame(edgelist_rfc4, directed = F, 
                                vertices = node_list_rfc4)
ig_rfc4 <- simplify(ig_rfc4, remove.multiple = TRUE, remove.loops = TRUE,
                    edge.attr.comb = igraph_opt("edge.attr.comb"))

# Coloring the node based on name of protein node
V(ig_rfc4)$color=ifelse(V(ig_rfc4)$protein =="pcna",
                        "#CC0099", "steelblue2")

sort(betweenness(ig_rfc4), decreasing = T)[1:5]

# ++++++++    rfc4 Full net node size definition   ++++++++++
rfc4_full_ns <- ifelse(betweenness(ig_rfc4) > 10.0, 
                       betweenness(ig_rfc4)/
                           (max(betweenness(ig_rfc4) * .06)), 2)

plot(ig_rfc4, layout=layout.fruchterman.reingold, 
     vertex.label=NA,
     vertex.label.cex= 0.8,
     #vertex.label.dist= 0.9,
     vertex.size=rfc4_full_ns,
     edge.arrow.size=0 ,vertex.label.color="black")



#### ************************
#   PCNA-POLD1
# ***************************

pold1_net <- read.csv("data/ POLD1_inter_network.csv")

## Getting the degree value of target Pold protein 
pold1_target_node_all <- pold1_net %>% group_by(Col2) %>% 
    summarise(freq=n())


#without filtering
pold1_net_gt4 <-pold1_net %>% group_by(Col1) %>% 
    summarise(freq=n())

pold1_net_gt4 <- pold1_net_gt4 %>% mutate(protein="pcna")
#pold_net_gt4$AA1 <- as.character(pold_net_gt4$AA1)

pold1_gt4_clean_net <- left_join(pold1_net_gt4, pold1_net)
pold1_gt4_clean_net <- pold1_gt4_clean_net %>% mutate(t2_node=paste0("P1-", Col2))
pold1_gt4_clean_net <- pold1_gt4_clean_net %>% select(Col1, t2_node, Correlation)

## Getting the degree value of target Pold protein 
pold1_gt4_target_node <- pold1_gt4_clean_net %>% group_by(t2_node) %>% 
    summarise(freq=n())

pold1_gt4_target_node$protein <- "pold1"
# renaming variable name in node file
names(pold1_gt4_target_node) <- c("Col1", "freq", "protein")


pold1_gt4_node_net <- pold1_net_gt4
pold1_gt4_node_net$Col1 <- as.character(pold1_gt4_node_net$Col1)
# Combining source and target node by row binding to total node list
tot_node_gt4 <- bind_rows(pold1_gt4_node_net, pold1_gt4_target_node)


# edgelist file without freq(degree) value
edgelist_pold1 <- pold1_gt4_clean_net

node_list_pold1 <- tot_node_gt4

ig_pold1 = graph_from_data_frame(edgelist_pold1, directed = F, vertices = node_list_pold1)
ig_pold1 <- simplify(ig_pold1, remove.multiple = TRUE, remove.loops = TRUE,
                     edge.attr.comb = igraph_opt("edge.attr.comb"))


sort(betweenness(ig_pold1), decreasing = T)[1:10]

# Coloring the node based on name of protein node
V(ig_pold1)$color=ifelse(V(ig_pold1)$protein =="pcna",
                         "#CC0099", "lightslateblue")


# ++++++++    pold1 Full net node size definition   ++++++++++

# pold1 Full net node size definition
pold1_full_ns <- ifelse(betweenness(ig_pold1) > 10.0, 
                        betweenness(ig_pold1)/
                            (max(betweenness(ig_pold1) * .075)), 3)



# ------------------------------------------------------------------------------------------

#### ************************
#   PCNA-FEN1
# ***************************

pcna_nnet <- read.csv("data/ PCNA_intra_network.csv")


# without filter
pcna_net <-pcna_nnet %>% group_by(AA1) %>% 
    summarise(freq=n()) 

pcna_node_net <- pcna_net %>% mutate(protein="pcna")
pcna_node_net$AA1 <- as.character(pcna_node_net$AA1)

# data frame with all degree node
pcna_clean_net <- left_join(pcna_net, pcna_nnet)
pcna_clean_net <- pcna_clean_net %>% mutate(t2_node=paste0("PCNA", AA2))
pcna_clean_net <- pcna_clean_net %>% select(AA1, t2_node, Correlation)

## Getting the degree value of target pcna protein 
pcna_target_node <-pcna_clean_net %>% group_by(t2_node) %>% 
    summarise(freq=n())

pcna_target_node$protein <- "pcna"
# renaming variable name in node file
names(pcna_target_node) <- c("AA1", "freq", "protein")

# Combining source and target node by row binding to total node list
pcna_tot_node <- bind_rows(pcna_node_net, pcna_target_node)

## Edgelist and network generation 
## *************************************************
# edgelist file without freq(degree) value
edgelist_pcna <- pcna_clean_net

node_list_pcna <- pcna_tot_node

# Creating graph from data frame and simplyfying
ig_pcna = graph_from_data_frame(edgelist_pcna, directed = F, vertices = node_list_pcna)

ig_pcna <- simplify(ig_pcna, remove.multiple = TRUE, remove.loops = TRUE,
                    edge.attr.comb = igraph_opt("edge.attr.comb"))


### Show all the colour schemes available
#display.brewer.all()

# R has 657 built in color names. To see
# colors()

#plot(ig_pcna, vertex.color="blueviolet", vertex.label=NA)

pcna_label_node <- sort(betweenness(ig_pcna), decreasing = T)[1:10]

# Colouring vertex according to protein name ************
V(ig_pcna)$color= "orange"
# Colouring vertex according to protein name ************
V(ig_pcna)$color=ifelse(V(ig_pcna)$protein =="pcna",
                        "#CC0099", "cyan2")


pcna_full_ns <- ifelse(betweenness(ig_pcna) > 10.0, betweenness(ig_pcna)/(max(betweenness(ig_pcna) * .075)), 3)

plot(ig_pcna, layout=layout.fruchterman.reingold,
     vertex.label=NA,
     #vertex.size=betweenness(ig_pcna)/(max(betweenness(ig_pcna) * .045)),
     vertex.size = pcna_full_ns,
     edge.arrow.size=0, vertex.label.color="black")


