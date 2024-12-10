# layout(matrix(1:10, 2, 5, byrow = TRUE))
# par(mar=c(3,5,3,1))
# set.seed(2017)
# The layout fruchterman.reingold.grid is similar to fruchterman.reingold, but faster.
#l <- layout.fruchterman.reingold(net, repulserad=vcount(net)^3,
#                                 area=vcount(net)^2.4)

# ----------------------------------------------------------------------------
library(Cairo)


## Fen1
# ***************************************
# Full Network
Cairo(file="Fen1_full_350_1.jpeg", 
      type="jpeg",
      units="in", 
      width=12, 
      height=9, 
      dpi=350)

plot(ig_fen1, layout=layout.fruchterman.reingold, 
     vertex.label=NA,
     edge.color=ecol_fen1_0,
     vertex.size=fen1_full_ns,
     edge.arrow.size=0, vertex.label.color="black")

dev.off()

# finalised Big cluster
Cairo(file="Fen1_big_clus_350_1.jpeg", 
      type="jpeg",
      units="in", 
      width=12, 
      height=9, 
      dpi=350)

plot(big_clus_fen1, layout=layout.fruchterman.reingold, 
     vertex.label=V(big_clus_fen1)$label_node,
     vertex.label.cex=0.9,
     edge.color=ecol_fen1,
     edge.curved=F,
     vertex.size=fen1_net_size,
     edge.arrow.size=0, vertex.label.color="black")


dev.off()

# ___________________________________________________________________

# POLD2 
# ***************************************

Cairo(file="Pold2_full_350_1.jpeg", 
      type="jpeg",
      units="in", 
      width=9, 
      height=12, 
      dpi=350)

plot(ig_pold2, layout=layout.fruchterman.reingold,
     vertex.label=NA,
     edge.color=ecol_pold2,
     vertex.size=pold2_full_ns,
     edge.arrow.size=0, vertex.label.color="#CCCCCC")


dev.off()

# ---------------------------------
Cairo(file="Pold2_big_clus_1.jpeg", 
      type="jpeg",
      units="in", 
      width=9, 
      height=12, 
      dpi=350)


plot(big_clus_pold2, layout=layout.fruchterman.reingold, 
     vertex.label.cex=0.9,
     vertex.label=V(big_clus_pold2)$label_node_pic,
     edge.curved=F,
     edge.color=ecol_pold2,
     vertex.label.degree= (-pi/2),
     vertex.size=pold2_bc_ns,
     edge.arrow.size=0, vertex.label.color="black")


dev.off()


# -----------------------------------------------------------

# **************************************
# RFC3 full network properties definition
# ***************************************

E(ig_rfc3)$width <- E(ig_rfc3)$Correlation*2.15
ecol_rfc3_0 <- rep("gray40", ecount(ig_rfc3))
V(ig_rfc3)$frame.color <- "black"

# **************************************
# RFC3 Plotting full network
# ***************************************


Cairo(file="rfc3_full_350_1.jpeg", 
      type="jpeg",
      units="in", 
      width=12, 
      height=9, 
      dpi=350)


plot(ig_rfc3, layout=layout.fruchterman.reingold, 
     vertex.label=NA, 
     edge.color=ecol_rfc3_0,
     vertex.size=rfc3_full_ns,
     edge.arrow.size=0, vertex.label.color="black")

dev.off()
# ---------------------------

Cairo(file="rfc3_350_big_clus_1.jpeg", 
      type="jpeg",
      units="in", 
      width=12, 
      height=9, 
      dpi=350)


plot(big_clus_rfc3, layout=layout.fruchterman.reingold, 
     vertex.label=V(big_clus_rfc3)$label_node_pic,
     vertex.label.cex=0.9,
     edge.color=ecol_rfc3,
     vertex.size=rfc3_bc_ns,
     edge.arrow.size=0, vertex.label.color="black")


dev.off()
# -------------------------------------------------

# **************************************
# RFC1 full network properties definition
# ***************************************
Cairo(file="rfc1_full_350_2.jpeg", 
      type="jpeg",
      units="in", 
      width=12, 
      height=9, 
      dpi=350)


plot(ig_rfc1, layout=layout.fruchterman.reingold, 
     vertex.label=NA, 
     edge.color=ecol_rfc1_0,
     vertex.size=rfc1_full_ns,
     edge.arrow.size=0, vertex.label.color="black")

dev.off()

# RFC1 Big clust plot
Cairo(file="rfc1_350_big_clus_3.jpeg", 
      type="jpeg",
      units="in", 
      width=12, 
      height=9, 
      dpi=350)


plot(big_clus_rfc1, layout=layout.fruchterman.reingold, 
     vertex.label=V(big_clus_rfc1)$label_node_pic,
     vertex.label.cex=0.9,
     edge.color=ecol_rfc1,
     vertex.size=rfc1_bc_ns,
     edge.arrow.size=0, vertex.label.color="black")


dev.off()


# **************************************
# RFC4 full network properties definition
# ***************************************

Cairo(file="rfc4_full_350_2.jpeg", 
      type="jpeg",
      units="in", 
      width=12, 
      height=9, 
      dpi=350)


plot(ig_rfc4, layout=layout.fruchterman.reingold, 
     vertex.label=NA, 
     edge.color=ecol_rfc4_0,
     vertex.size=rfc4_full_ns,
     edge.arrow.size=0, vertex.label.color="black")

dev.off()

# rfc4 Big clust plot
Cairo(file="rfc4_350_big_clus_4.jpeg", 
      type="jpeg",
      units="in", 
      width=12, 
      height=9, 
      dpi=350)


plot(big_clus_rfc4, layout=layout.fruchterman.reingold, 
     vertex.label=V(big_clus_rfc4)$label_node_pic,
     vertex.label.cex=0.9,
     edge.color=ecol_rfc4,
     vertex.size=rfc4_bc_ns,
     edge.arrow.size=0, vertex.label.color="black")


dev.off()





# ___________________________________________________________________

# POLD1
# ***************************************
Cairo(file="pold1_full_350_1.jpeg", 
      type="jpeg",
      units="in", 
      width=9, 
      height=12, 
      dpi=350)

plot(ig_pold1, layout=layout.fruchterman.reingold,
     vertex.label=NA,
     edge.color=ecol_pold1_0,
     vertex.size=pold1_full_ns,
     edge.arrow.size=0, vertex.label.color="#CCCCCC")


dev.off()

# ---------------------------------
Cairo(file="pold1_big_clus_1.jpeg", 
      type="jpeg",
      units="in", 
      width=9, 
      height=12, 
      dpi=350)


plot(big_clus_pold1, layout=layout.fruchterman.reingold, 
     vertex.label.cex=0.9,
     vertex.label=V(big_clus_pold1)$label_node_pic,
     edge.curved=F,
     edge.color=ecol_pold1,
     vertex.label.degree= (-pi/2),
     vertex.size=pold1_bc_ns,
     edge.arrow.size=0, vertex.label.color="black")


dev.off()



# -----------------------------------------------------------
