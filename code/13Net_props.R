## *******************************
## NETWORK PROPERTIES
## *******************************

avg_path_length_pold <- average.path.length(ig_pold, directed=F, unconnected=TRUE)
#average.path.length(ig_pold, directed=F)

#/calculates the vertex betweenness, directed = TRUE, weights = NULL, nobigint = TRUE, normalized = FALSE,
betweenness_pold<-betweenness(ig_pold,V(ig_pold),FALSE,NULL,TRUE,FALSE)
which.max(betweenness_pold)

## //calucalating the degree distribution for all nodes in a graph
dd_pold<-degree.distribution(ig_pold)
## //plotting the degree distribution
plot(dd_pold)
dd_pold
## //calucalating the degree distribution for all nodes in a graph
dd_fen1<-degree.distribution(ig_fen1)
plot(dd_fen1)
hist(log2(degree(ig_pold)))
cliques(ig_pold)
largest.cliques(ig_pold)
hist(pold_net$Correlation)
mean(pold_net$Correlation)
sd(pold_net$Correlation)
# edge betweenness
eb <- cluster_edge_betweenness(ig_pold)
eb[7]
dim(eb)

# Centrality Measures
# At the fine grain level, we can look at statistics of individual nodes.  Centrality score measure the social importance of a node in terms of how "central" it is based on a number of measures ...
# Degree centrality gives a higher score to a node that has a high in/out-degree
# Closeness centrality gives a higher score to a node that has short path distance to every other nodes
# Betweenness centrality gives a higher score to a node that sits on many shortest path of other node pairs
# Eigenvector centrality gives a higher score to a node if it connects to many high score nodes
# Local cluster coefficient measures how my neighbors are inter-connected with each other, which means the node becomes less important.
# Degree
degree(ig_pold)>4

# Closeness (inverse of average dist)
closeness(ig_pold)

# Betweenness
betweenness(ig_pold)
# hub score
hub.score(ig_pold)$vector
# Local cluster coefficient
transitivity(ig_pold, type="local")

# Eigenvector centrality
evcent(ig_pold)$vector

# Now rank them
order(degree(ig_pold))

order(closeness(ig_pold))

order(betweenness(ig_pold))

order(evcent(ig_pold)$vector)

# From his studies, Drew Conway has found that people with low Eigenvector centrality but high Betweenness centrality are important gate keepers, while people with high Eigenvector centrality but low Betweenness centrality has direct contact to important persons.  So lets plot Eigenvector centrality against Betweenness centrality.
# Create a graph
g1 <- barabasi.game(100, directed=F)
g2 <- barabasi.game(100, directed=F)
g <- g1 %u% g2
lay <- layout.fruchterman.reingold(g)
# Plot the eigevector and betweenness centrality
plot(evcent(g)$vector, betweenness(g))
text(evcent(g)$vector, betweenness(g), 0:100, 
     cex=0.6, pos=4)
V(g)[12]$color <- 'red'
V(g)[8]$color <- 'green'
plot(g, layout=lay, vertex.size=8, 
     vertex.label.cex=0.6)

# Community strucure via short random walks
# Description
# This function tries to find densely connected subgraphs, also called communities in a graph via random walks. The idea is that short random walks tend to stay in the same community.


comm_strc <- cluster_walktrap(ig_pold, weights = E(ig_pold)$Correlation, steps = 2,
                              merges = TRUE, modularity = TRUE, membership = TRUE)

## COMMUNITY
membership(comm_strc )
sizes(comm_strc)
show_trace(comm_strc )

#components <- clusters(graph, mode = "strong")
components_N <- clusters(ig_pold, mode = "strong")
#print(components$csize)
print(components_N$no)
#qplot(degree(G3), geom = "histogram")
qplot(degree(ig_pold), geom = "histogram")
## ggplot(data=degree(ig_pold)) + geom_histogram()

library("qgraph")
# **** It will draw 3 centrality plot together 
## // very convenient way to plot *******

centralityPlot(Network)
