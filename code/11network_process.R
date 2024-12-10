# Fen1 full network properties definition
# ***************************************
library(bc3net)
library(Cairo)

E(ig_fen1)$width <- E(ig_fen1)$Correlation*2.15
ecol_fen1_0 <- rep("gray40", ecount(ig_fen1))
V(ig_fen1)$frame.color <- "black"

# Fen1 Big cluster properties definition
# ***************************************

big_clus_fen1 <- getgcc(ig_fen1)
plot(big_clus_fen1)

# sort and print out the decreasing betweenness value 
sort(betweenness(big_clus_fen1), decreasing = T)[1:10]

# Fen1 full net node size definition
fen1_fullnet_size <- ifelse(betweenness(ig_fen1)> 10.00, betweenness(ig_fen1)/(max(betweenness(ig_fen1) * .055)), 3)

# fen1_full_ns <- ifelse(betweenness(ig_fen1) > 10.0, betweenness(ig_fen1)/(max(betweenness(ig_fen1) * .075)), 3)

# Fen1 Big cluster node size definition
fen1_net_size=ifelse(betweenness(big_clus_fen1)> 10.0, betweenness(big_clus_fen1)/(max(betweenness(big_clus_fen1)* .055 )), 5)

sort(betweenness(big_clus_fen1), decreasing = T)[1:10]
keep_fen1_label_pic <- as.data.frame(c("L19","S32", "E313", "Q125", "L258", "E275"))
names(keep_fen1_label_pic) <- "PIC_ID"


keep_fen1_label <- as.data.frame(c("19","32", "F313", "125", "F258", "F275"))
names(keep_fen1_label) <- "ID"


V(big_clus_fen1)$label_node_pic <- as.character(keep_fen1_label_pic$PIC_ID
                                            [match(V(big_clus_fen1)$name,
                                                   keep_fen1_label$ID)])

E(big_clus_fen1)$width <- E(big_clus_fen1)$Correlation*2.15
ecol_fen1 <- rep("gray40", ecount(big_clus_fen1))
V(big_clus_fen1)$frame.color <- "white"

## Fen1
Cairo(file="Fen1_big_clus.jpeg", 
      type="jpeg",
      units="in", 
      width=12, 
      height=9, 
      dpi=350)

plot(big_clus_fen1, layout=layout.fruchterman.reingold, 
     vertex.label=V(big_clus_fen1)$label_node_pic,
     vertex.label.cex=0.8,
     edge.color=ecol_fen1,
     edge.curved=F,
     vertex.size=fen1_net_size,
     edge.arrow.size=0, vertex.label.color="black")

dev.off()

# POLD2  ***********************

# ***************************************
# Pold full network properties definition
# ***************************************

E(ig_pold2)$width <- E(ig_pold2)$Correlation*2.15
ecol_pold2_0 <- rep("gray40", ecount(ig_pold2))
V(ig_pold2)$frame.color <- "black"

# Big cluster POLD2
big_clus_pold2 <- getgcc(ig_pold2)
plot(big_clus_pold2)

sort(betweenness(big_clus_pold2), decreasing = T)[1:10]
# 125    P2-302       137        32    P2-268    P2-104    P2-195     P2-23    P2-239    P2-252 
# 68.333333 26.833333 26.166667 14.500000 11.500000  7.666667  0.000000  0.000000  0.000000  0.000000 

keep_pold2_label_pic <- as.data.frame(c("Q125","T302", "V137", "S32", "T268", "Q104"))
names(keep_pold2_label_pic) <- "PIC_ID"

# node_list pold2 have these names
keep_pold2_label <- as.data.frame(c("125","P2-302", "137", "32", "P2-268", "P2-104"))
names(keep_pold2_label) <- "ID"

V(big_clus_pold2)$label_node_pic <- as.character(keep_pold2_label_pic$PIC_ID
                                                [match(V(big_clus_pold2)$name,
                                                       keep_pold2_label$ID)])

E(big_clus_pold2)$width <- E(big_clus_pold2)$Correlation*2.15
ecol_pold2 <- "gray40"
#ecol_pold2 <- rep("gray40", ecount(big_clus_pold2))
V(big_clus_pold2)$frame.color <- "white"

# POLD2 Big cluster node size definition
pold2_bc_ns=ifelse(betweenness(big_clus_pold2)> 10.0, 
                   betweenness(big_clus_pold2)/(
                       max(betweenness(big_clus_pold2)))*25, 3)

# RFC3 ***********************************************************
# ****************************************************************

E(ig_rfc3)$width <- E(ig_rfc3)$Correlation*2.15
ecol_rfc3_0 <- rep("gray40", ecount(ig_rfc3))
V(ig_rfc3)$frame.color <- "black"

# Big cluster RFC3
big_clus_rfc3 <- getgcc(ig_rfc3)
plot(big_clus_rfc3)

sort(betweenness(big_clus_rfc3), decreasing = T)[1:10]
# 19 R3-139    258 R3-125 R3-171 R3-233 R3-243 R3-260 R3-269 R3-284 
# 88     24     13      0      0      0      0      0      0      0 

keep_rfc3_label_pic <- as.data.frame(c("L19","K139", "E258"))
names(keep_rfc3_label_pic) <- "PIC_ID"

# node_list pold2 have these names
keep_rfc3_label <- as.data.frame(c("19","R3-139", "258"))
names(keep_rfc3_label) <- "ID"

V(big_clus_rfc3)$label_node_pic <- as.character(keep_rfc3_label_pic$PIC_ID
                                                 [match(V(big_clus_rfc3)$name,
                                                        keep_rfc3_label$ID)])


E(big_clus_rfc3)$width <- E(big_clus_rfc3)$Correlation*2.15
ecol_rfc3 <- "gray40"
#ecol_pold2 <- rep("gray40", ecount(big_clus_pold2))
V(big_clus_rfc3)$frame.color <- "white"

# RFC3 Big cluster node size definition
rfc3_bc_ns=ifelse(betweenness(big_clus_rfc3)> 10.0, 
                   betweenness(big_clus_rfc3)/(
                       max(betweenness(big_clus_rfc3)))*35, 3)


# RFC1 ***********************************************************
# ****************************************************************
E(ig_rfc1)$width <- E(ig_rfc1)$Correlation*2.15
ecol_rfc1_0 <- rep("gray40", ecount(ig_rfc1))
V(ig_rfc1)$frame.color <- "black"

# Big cluster RFC3
big_clus_rfc1 <- getgcc(ig_rfc1)
plot(big_clus_rfc1)

sort(betweenness(big_clus_rfc1), decreasing = T)[1:10]
# 32       125        19    R1-972       163       256    R1-675    R1-240        71    R1-836 
# 8138.1858 7148.1932 1965.6010  433.2791  351.1077  349.9819  306.2299  282.1956  257.7012  215.4353 
#rfc1_prot[c(972,675, 240, 836 )]
# PVKA
# HLKD
#pcna_prot[c(32, 125, 19, 163, 256, 71)]
#SQLAEN

keep_rfc1_label_pic <- as.data.frame(c("S32","Q125", "L19", "H972", "A163", "E256", "L675","K240", "N71","D836" ))
names(keep_rfc1_label_pic) <- "PIC_ID"

# node_list rfc1 have these names
keep_rfc1_label <- as.data.frame(c("32","125", "19", "R1-972", "163", "256", "R1-675","R1-240", "71","R1-836" ))
names(keep_rfc1_label) <- "ID"

V(big_clus_rfc1)$label_node_pic <- as.character(keep_rfc1_label_pic$PIC_ID
                                                 [match(V(big_clus_rfc1)$name,
                                                        keep_rfc1_label$ID)])

E(big_clus_rfc1)$width <- E(big_clus_rfc1)$Correlation*2.15
ecol_rfc1 <- "gray40"
#ecol_pold2 <- rep("gray40", ecount(big_clus_pold2))
V(big_clus_rfc1)$frame.color <- "white"

# RFC1 Big cluster node size definition
rfc1_bc_ns=ifelse(betweenness(big_clus_rfc1)> 10.0, 
                  betweenness(big_clus_rfc1)/(
                      max(betweenness(big_clus_rfc1)))*20, 2.5)




# RFC4 ***********************************************************
# ****************************************************************
E(ig_rfc4)$width <- E(ig_rfc4)$Correlation*2.15
ecol_rfc4_0 <- rep("gray40", ecount(ig_rfc4))
V(ig_rfc4)$frame.color <- "black"

# Big cluster RFC4
big_clus_rfc4 <- getgcc(ig_rfc4)
plot(big_clus_rfc4)

sort(betweenness(big_clus_rfc4), decreasing = T)[1:10]
# 174    192 R4-102 R4-117 R4-353 R4-356 R4-133 R4-136 R4-188 R4-204 
# 264.0   26.0    9.5    9.5    9.5    9.5    0.0    0.0    0.0    0.0 

keep_rfc4_label_pic <- as.data.frame(c("E174","E192", "R102", "Q117", "A353", "M356"))
names(keep_rfc4_label_pic) <- "PIC_ID"

# node_list rfc1 have these names
keep_rfc4_label <- as.data.frame(c("174","192", "102", "117", "353", "356"))
names(keep_rfc4_label) <- "ID"

V(big_clus_rfc4)$label_node_pic <- as.character(keep_rfc4_label_pic$PIC_ID
                                                [match(V(big_clus_rfc4)$name,
                                                       keep_rfc4_label$ID)])

E(big_clus_rfc4)$width <- E(big_clus_rfc4)$Correlation*2.15
ecol_rfc4 <- "gray40"
#ecol_pold2 <- rep("gray40", ecount(big_clus_pold2))
#V(big_clus_rfc4)$frame.color <- "white"

# rfc4 Big cluster node size definition
rfc4_bc_ns=ifelse(betweenness(big_clus_rfc4)> 9.0, 
                  betweenness(big_clus_rfc4)/(
                      max(betweenness(big_clus_rfc4)))* 50,2.5)




# POLD1  ***********************

# ***************************************
# Pold1 full network properties definition
# ***************************************
E(ig_pold1)$width <- E(ig_pold1)$Correlation*2.15
ecol_pold1_0 <- rep("gray40", ecount(ig_pold1))
V(ig_pold1)$frame.color <- "black"

# Big cluster pold1
big_clus_pold1 <- getgcc(ig_pold1)
plot(big_clus_pold1)

sort(betweenness(big_clus_pold1), decreasing = T)[1:10]
# 49 P1-232 P1-573 P1-649   <NA>   <NA>   <NA>   <NA>   <NA>   <NA> 
# 3      0      0      0     NA     NA     NA     NA     NA     NA 

keep_pold1_label_pic <- as.data.frame(c("Q49"))
names(keep_pold1_label_pic) <- "PIC_ID"

# node_list pold1 have these names
keep_pold1_label <- as.data.frame(c("49"))
names(keep_pold1_label) <- "ID"

V(big_clus_pold1)$label_node_pic <- as.character(keep_pold1_label_pic$PIC_ID
                                                 [match(V(big_clus_pold1)$name,
                                                        keep_pold1_label$ID)])

E(big_clus_pold1)$width <- E(big_clus_pold1)$Correlation*2.15
ecol_pold1 <- "gray40"
#ecol_pold1 <- rep("gray40", ecount(big_clus_pold1))
V(big_clus_pold1)$frame.color <- "white"

# pold1 Big cluster node size definition
pold1_bc_ns=ifelse(betweenness(big_clus_pold1)> 10.0, 
                   betweenness(big_clus_pold1)/(
                       max(betweenness(big_clus_pold1)))*25, 4)



# -------------------------------------------------------------------------------------
# ***************************************
# Pold full network properties definition
# ***************************************

E(ig_pcna)$width <- E(ig_pcna)$Correlation*2.15
ecol_pcna_0 <- rep("gray40", ecount(ig_pcna))
V(ig_pcna)$frame.color <- "black"

# Big cluster pcna
big_clus_pcna <- getgcc(ig_pcna)
plot(big_clus_pcna)

sort(betweenness(big_clus_pcna), decreasing = T)[1:10]
# 125    P2-302       137        32    P2-268    P2-104    P2-195     P2-23    P2-239    P2-252 
# 68.333333 26.833333 26.166667 14.500000 11.500000  7.666667  0.000000  0.000000  0.000000  0.000000 

keep_pcna_label_pic <- as.data.frame(c("Q125","T302", "V137", "S32", "T268", "Q104"))
names(keep_pcna_label_pic) <- "PIC_ID"

# node_list pcna have these names
keep_pcna_label <- as.data.frame(c("125","P2-302", "137", "32", "P2-268", "P2-104"))
names(keep_pcna_label) <- "ID"

V(big_clus_pcna)$label_node_pic <- as.character(keep_pcna_label_pic$PIC_ID
                                                [match(V(big_clus_pcna)$name,
                                                       keep_pcna_label$ID)])

E(big_clus_pcna)$width <- E(big_clus_pcna)$Correlation*2.15
ecol_pcna <- "gray40"
#ecol_pcna <- rep("gray40", ecount(big_clus_pcna))
V(big_clus_pcna)$frame.color <- "white"

# pcna Big cluster node size definition
pcna_bc_ns=ifelse(betweenness(big_clus_pcna)> 10.0, 
                  betweenness(big_clus_pcna)/(
                      max(betweenness(big_clus_pcna)))*25, 3)