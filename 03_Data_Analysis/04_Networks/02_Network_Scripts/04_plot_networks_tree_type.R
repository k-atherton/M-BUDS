library(stringr)
library(igraph)

plot_network <- function(network, site){
  # Isolated = which(degree(network)==0)
  # network = delete_vertices(network, Isolated)
  layout <- layout.fruchterman.reingold(network)
  png(filename = paste0("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/",
                        site,".png"), 
      width = 6, height = 4, units = "in", res = 600)
  par(mar=c(1, 1, 2, 1))
  plot(network, 
       layout=layout, 
       vertex.color=V(network)$color, 
       vertex.size=degree(network), 
       vertex.size=3,
       vertex.label=NA, 
       vertex.shape=V(network)$shape,
       edge.width=E(network)$weight^4,
       main=paste0(site))
  legend("topright",bty = "n",
         legend=levels(as.factor(sapply(str_wrap(V(network)$guild,28),unlist))),
         fill=unique(V(network)[order(V(network)$guild)]$color), border=NA, cex = 0.35)
  # theme(legend.key.height = unit(1, "cm"))
  invisible(dev.off())
}

leaf_network_atts <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/leaf_network_atts_treetype.RDS")
root_network_atts <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/root_network_atts.RDS")
m_network_atts <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/m_network_atts.RDS")
o_network_atts <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/o_network_atts.RDS")



sites_leaf <- c("leaf_pit", "leaf_grass", "leaf_urban_edge", "leaf_urban_interior", "leaf_rural_edge", "leaf_rural_interior")
for(i in 1:length(leaf_network_atts)){
  plot_network(leaf_network_atts[[i]], sites_leaf[[i]])
}

sites_root <- c("root_pit", "root_grass", "root_urban_edge", "root_urban_interior", "root_rural_edge", "root_rural_interior")
for(i in 1:length(root_network_atts)){
  plot_network(root_network_atts[[i]], sites_root[[i]])
}

sites_m <- c("m_pit", "m_grass", "m_urban_edge", "m_urban_interior", "m_rural_edge", "m_rural_interior")
for(i in 1:length(m_network_atts)){
  plot_network(m_network_atts[[i]], sites_m[[i]])
}

sites_o <- c("o_pit", "o_grass", "o_urban_edge", "o_urban_interior", "o_rural_edge", "o_rural_interior")
for(i in 1:length(o_network_atts)){
  plot_network(o_network_atts[[i]], sites_o[[i]])
}

