library(igraph)
library(dplyr)
library(readr)
library(bc3net)
library(tidyr)
setwd("~/suny/RESEARCH_DOC/A_THESIS_WORK/CHAPTER_5_COEVOLV_RESIDUE/eukarya")
total_network <- load("fen1_pold_rfc3_network.RData")
#average.path.length(ig_fen1)
## Function to create Network properties Data frame
## **********************************************************
Network_properties <- function(network_name, protein_complex){
    graph.obj <- network_name
 
    # Average path length
    avg_path_length <- average.path.length(graph.obj)
    # Vertex
    node_len <- length(V(graph.obj))
    # Diameter
    dia <- diameter(graph.obj, directed=F)
    #Edge count
    edge_count <- ecount(graph.obj)
    # No. of cluster
    clust_no <- clusters(graph.obj)$no
    #shortest.paths(g)
    #shortest_paths(g)
    # Density
    edge_dens <- edge_density(graph.obj)
    # Degree
    deg <- degree(graph.obj, mode = "all")
    # Mean degree distribution
    mean_deg <- mean(degree(graph.obj))
    # Degree Centralisation
    deg_centralisation <- centr_degree(graph.obj, mode="all", normalized=T)
    deg_cent <- deg_centralisation$centralization
    # Correlation of attributes across connected nodes
    deg_assort <- assortativity_degree(graph.obj, directed=F)
    # Local cluster coefficient
    clus_coef <- transitivity(graph.obj, type="local")
    # Closeness Centralisation
    # order(closeness(ig_fen1,mode = "all", weights=NA)*100000)
    closeness_centralisation <- centr_clo(graph.obj,mode = "all",
                                          normalized = T)
    close_cent <- closeness_centralisation$centralization
## Data frame output
df <- data.frame(avg_path_length, node_len, dia, edge_count, clust_no, edge_dens,  mean_deg, deg_cent, deg_assort, close_cent)
 
    
    # Closeness (inverse of average dist)
    close <- closeness(graph.obj,mode = "all", weights=NA)
    # Betweenness
    btn <- betweenness(graph.obj)
    # Eigenvector centrality
    eig_vec <- evcent(graph.obj)$vector
    
    
#df_2 <- data.frame(close, btn, eig_vec )

 #colnames(df) <- protein_complex
    return (df)
   # return (df_2)
}

fen1_net_prop <- Network_properties(ig_fen1, PCNA-Fen1)

pold_net_prop <- Network_properties(ig_pold, PCNA-Pold)

rfc3_net_prop <- Network_properties(ig_rfc3, PCNA-Rfc3)

net_df <- data.frame(fen1_net_prop, pold_net_prop, rfc3_net_prop)

# Creating Network Properties table for Manuscript
library(dplyr)
n_df <- bind_rows(fen1_net_prop, pold_net_prop, rfc3_net_prop)
dim(n_df)
head(n_df)

long_df <- data.frame(t(n_df))
long_df %>% 
names(long_df) <- c("Fen1", "Pold", "Rfc3")
long_df %>% select(Fen1) %>% as.numeric
long_df$Fen1 <- as.numeric(long_df$Fen1)
write.csv(long_df, file="network_properties.csv")

head(long_df)
n_df %>% gather(fen1, pold, rfc3, na.rm=F, convert=F)
# ------------------------------------------------
# Shortest paths
distances(graph.obj, weights=NA) # ignore weights
fen1_net_prop <- Network_properties(ig_fen1, "PCNA-Fen1")
eb <- edge.betweenness.community(ig_fen1)
eb
# Degree Centralisation
deg_centralisation <- centr_degree(ig_fen1, mode="all", normalized=T)
deg_centralisation$centralization

# Closeness Centralisation
order(closeness(ig_fen1,mode = "all", weights=NA)*100000)
closeness_centralisation <- centr_clo(ig_fen1,mode = "all", normalized = T)
closeness_centralisation$centralization

# Correlation of attributes across connected nodes
assortativity_degree(ig_fen1, directed=F)
knnnet <- knn(graph = ig_fen1)

plot(knnnet$knn, xlab = "k", ylab = "knn(k)", main = "Degree correlation function lin-lin")

plot(knnnet$knn, xlab = "k", ylab = "knn(k)", log = "xy", main = "Degree correlation function log-log")


myc <- clusters(ig_fen1, mode="strong")
myc$csize
m1 <- myc$membership==1
hs <- hub_score(ig_fen1, weights=NA)$vector
as <- authority_score(ig_fen1, weights=NA)$vector
par(mfrow=c(1,2))
plot(ig_fen1, vertex.size=hs*10, main="Hubs")
plot(ig_fen1, vertex.size=as*10, main="Authorities")
mean_distance(ig_fen1, directed=F)
ebc <- edge.betweenness.community(ig_fen1)
modularity(ig_fen1,membership(ebc))

#Average nearest neighbor degree
av_nnd <- knn(ig_fen1, vids = V(ig_fen1), weights = NULL)
histogram(av_nnd$knnk)
str(av_nnd)
mean(degree(ig_fen1))
#degree distribution

dd2 <- degree.distribution(ig_fen1)
#dd2 is calculated for a node degree ranging from 0 to the max degree found in the network, hence I
#define the following
k_range <- 0:(length(dd2)-1)
#------------
deg <- degree(ig_fen1, mode="all")
deg.dist <- degree_distribution(ig_fen1, cumulative=T, mode="all") 
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange",
                                                                     xlab="Degree", ylab="Cumulative Frequency", main="PCNA-Fen1")

# -----------------------

#Get also the cumulative distribution

dd2_cum <- degree.distribution(ig_fen1, cumulative=TRUE)
plot(k_range,dd2, xlab="Degree",
     ylab="Frequency", pch=3, col=3, type="b",cex.axis=1.4,cex.lab=1.6)

distances(ig_fen1, weights=NA) 

#--------------
big_clus_fen1 <- getgcc(ig_fen1)
plot(big_clus_fen1)

big_clus_pold <- getgcc(ig_pold)
plot(big_clus_pold)

big_clus_rfc3 <- getgcc(ig_rfc3)
plot(big_clus_rfc3)

plot(big_clus_pold, 
     layout=layout.fruchterman.reingold,
     vertex.label.cex=.5,
     vertex.size=betweenness(big_clus_pold)
     /(max(betweenness(big_clus_pold) * .055)),
     edge.arrow.size=0, vertex.label.color="black")

plot(big_clus, 
     layout=layout.fruchterman.reingold,
     vertex.label.cex=.5,
     vertex.size=betweenness(ig_fen1)
     /(max(betweenness(ig_fen1) * .055)),
     edge.arrow.size=0, vertex.label.color="black")

plot(betweenness(big_clus_pold))

# Closeness (inverse of average dist)
close <- closeness(ig_fen1,mode = "all", weights=NA)
plot(close)
dev.off()
# Betweenness
btn <- betweenness(ig_fen1)
plot(btn, log="xy")
