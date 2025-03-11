library(stringr)
library(igraph)

plot_network <- function(network, site){
  Isolated = which(degree(network)==0)
  network = delete_vertices(network, Isolated)
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

leaf_network_atts <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/leaf_network_atts.RDS")
root_network_atts <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/root_network_atts.RDS")
m_network_atts <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/m_network_atts.RDS")
o_network_atts <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/o_network_atts.RDS")



sites_leaf <- c("leaf_aa01", "leaf_bm03", "leaf_hf04", "leaf_hf06", "leaf_hf069c", "leaf_hf069d", "leaf_hw07", "leaf_jp", "leaf_sw08")
for(i in 1:length(leaf_network_atts)){
  plot_network(leaf_network_atts[[i]], sites_leaf[[i]])
}

sites_root <- c("root_aa01", "root_bm03", "root_dosb", "root_hf069c", "root_hw07", "root_jp", "root_mcse", "root_sw08", "root_wr")
for(i in 1:length(root_network_atts)){
  plot_network(root_network_atts[[i]], sites_root[[i]])
}

sites_m <- c("m_aa01", "m_abbb", "m_bm03", "m_dosb", "m_hf04", "m_hf06", "m_hf069_c", "m_hf069d", "m_hf069e", "m_hw07", "m_jp", "m_mcse", "m_mp", "m_wr")
for(i in 1:length(m_network_atts)){
  plot_network(m_network_atts[[i]], sites_m[[i]])
}

sites_o <- c("o_aa01", "o_abbb", "o_bm03", "o_dosb", "o_hf04", "o_hf069c", "o_hf069d", "o_hf069_e", "o_hw07", "o_jp", "o_mcse", "o_mp", "o_sw08", "o_wr")
for(i in 1:length(o_network_atts)){
  plot_network(o_network_atts[[i]], sites_o[[i]])
}

