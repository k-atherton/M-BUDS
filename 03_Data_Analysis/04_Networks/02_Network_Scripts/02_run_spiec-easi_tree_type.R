library(SpiecEasi)

leaf_its <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/leaf_its_treetype.RDS")
root_its <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/root_its_treetype.RDS")
m_its <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/m_its_treetype.RDS")
o_its <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/o_its_treetype.RDS")

leaf_16s <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/leaf_16s_treetype.RDS")
root_16s <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/root_16s_treetype.RDS")
m_16s <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/m_16s_treetype.RDS")
o_16s <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/o_16s_treetype.RDS")

source('/projectnb2/talbot-lab-data/Katies_data/UNENetworks/UNE-Networks/Functions/04_Run_SPIEC-EASI.R')

leaf_networks <- vector("list", 6)
for(i in 1:length(leaf_networks)){
  print(paste0("Starting network ", i, ":"))
  a <- leaf_its[[i]]
  b <- leaf_16s[[i]]
  samples <- intersect(colnames(a), colnames(b))
    a <- a[,colnames(a) %in% samples]
    b <- b[,colnames(b) %in% samples]
    a <- t(a)
    b <- t(b)
    leaf_networks[[i]] <- network_cross(a, b, 0.1, 500)
    print(paste0("Network ", i, " complete.")) 
}
saveRDS(leaf_networks, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/leaf_networks.RDS")

root_networks <- vector("list", 6)
for(i in 1:length(root_networks)){
  print(paste0("Starting network ", i, ":"))
  a <- root_its[[i]]
  b <- root_16s[[i]]
  samples <- intersect(colnames(a), colnames(b))
    a <- a[,colnames(a) %in% samples]
    b <- b[,colnames(b) %in% samples]
    a <- t(a)
    b <- t(b)
    root_networks[[i]] <- network_cross(a, b, 0.1, 700)
    print(paste0("Network ", i, " complete.")) 
}

saveRDS(root_networks, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/root_networks.RDS")

m_networks <- vector("list", 6)
for(i in 1:length(m_networks)){
  print(paste0("Starting network ", i, ":"))
  a <- m_its[[i]]
  b <- m_16s[[i]]
  samples <- intersect(colnames(a), colnames(b))
    a <- a[,colnames(a) %in% samples]
    b <- b[,colnames(b) %in% samples]
    a <- t(a)
    b <- t(b)
    m_networks[[i]] <- network_cross(a, b, 0.1, 500)
    print(paste0("Network ", i, " complete.")) 
}
saveRDS(m_networks, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/m_networks.RDS")

o_networks <- vector("list", 6)
for(i in 1:length(o_networks)){
  print(paste0("Starting network ", i, ":"))
  a <- o_its[[i]]
  b <- o_16s[[i]]
  a <- a[,order(colnames(a))]
  b <- b[,order(colnames(b))]
  samples <- intersect(colnames(a), colnames(b))
  a <- a[,colnames(a) %in% samples]
  b <- b[,colnames(b) %in% samples]
  a <- t(a)
  b <- t(b)
  o_networks[[i]] <- network_cross(a, b, 0.1, 500)
  print(paste0("Network ", i, " complete.")) 
}
saveRDS(o_networks, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/o_networks.RDS")
