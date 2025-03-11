library(dplyr)

# load data 
leaf_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/leaf_its_genus.csv", row.names = 1)
root_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/root_its_genus.csv", row.names = 1)
m_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/m_its_genus.csv", row.names = 1)
o_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/o_its_genus.csv", row.names = 1)

leaf_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/leaf_16s_genus.csv", row.names = 1)
root_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/root_16s_genus.csv", row.names = 1)
m_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/m_16s_genus.csv", row.names = 1)
o_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/o_16s_genus.csv", row.names = 1)

metadata_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Data/metadata_its.csv")
metadata_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Data/metadata_16s_200.csv")

metadata_its <- metadata_its[!is.na(metadata_its$sample_name),]
metadata_16s <- metadata_16s[!is.na(metadata_16s$sample_name),]

metadata_its <- metadata_its[which(metadata_its$sequences_dropped == "No"),]
metadata_16s <- metadata_16s[which(metadata_16s$sequences_dropped == "No"),]

# get names for each tree type and sample type
leaf_its_pit_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$tree_pit_type == "Pit")]
leaf_its_grass_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$tree_pit_type == "Grass")]
leaf_its_urban_edge_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$tree_pit_type == "Urban Edge")]
leaf_its_urban_int_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$tree_pit_type == "Urban Interior")]
leaf_its_rural_edge_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$tree_pit_type == "Rural Edge")]
leaf_its_rural_int_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$tree_pit_type == "Rural Interior")]

root_its_pit_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$tree_pit_type == "Pit")]
root_its_grass_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$tree_pit_type == "Grass")]
root_its_urban_edge_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$tree_pit_type == "Urban Edge")]
root_its_urban_int_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$tree_pit_type == "Urban Interior")]
root_its_rural_edge_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$tree_pit_type == "Rural Edge")]
root_its_rural_int_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$tree_pit_type == "Rural Interior")]

m_its_pit_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$tree_pit_type == "Pit")]
m_its_grass_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$tree_pit_type == "Grass")]
m_its_urban_edge_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$tree_pit_type == "Urban Edge")]
m_its_urban_int_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$tree_pit_type == "Urban Interior")]
m_its_rural_edge_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$tree_pit_type == "Rural Edge")]
m_its_rural_int_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$tree_pit_type == "Rural Interior")]

o_its_pit_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$tree_pit_type == "Pit")]
o_its_grass_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$tree_pit_type == "Grass")]
o_its_urban_edge_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$tree_pit_type == "Urban Edge")]
o_its_urban_int_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$tree_pit_type == "Urban Interior")]
o_its_rural_edge_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$tree_pit_type == "Rural Edge")]
o_its_rural_int_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$tree_pit_type == "Rural Interior")]

leaf_16s_pit_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$tree_pit_type == "Pit")]
leaf_16s_grass_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$tree_pit_type == "Grass")]
leaf_16s_urban_edge_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$tree_pit_type == "Urban Edge")]
leaf_16s_urban_int_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$tree_pit_type == "Urban Interior")]
leaf_16s_rural_edge_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$tree_pit_type == "Rural Edge")]
leaf_16s_rural_int_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$tree_pit_type == "Rural Interior")]

root_16s_pit_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$tree_pit_type == "Pit")]
root_16s_grass_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$tree_pit_type == "Grass")]
root_16s_urban_edge_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$tree_pit_type == "Urban Edge")]
root_16s_urban_int_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$tree_pit_type == "Urban Interior")]
root_16s_rural_edge_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$tree_pit_type == "Rural Edge")]
root_16s_rural_int_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$tree_pit_type == "Rural Interior")]

m_16s_pit_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$tree_pit_type == "Pit")]
m_16s_grass_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$tree_pit_type == "Grass")]
m_16s_urban_edge_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$tree_pit_type == "Urban Edge")]
m_16s_urban_int_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$tree_pit_type == "Urban Interior")]
m_16s_rural_edge_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$tree_pit_type == "Rural Edge")]
m_16s_rural_int_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$tree_pit_type == "Rural Interior")]

o_16s_pit_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$tree_pit_type == "Pit")]
o_16s_grass_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$tree_pit_type == "Grass")]
o_16s_urban_edge_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$tree_pit_type == "Urban Edge")]
o_16s_urban_int_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$tree_pit_type == "Urban Interior")]
o_16s_rural_edge_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$tree_pit_type == "Rural Edge")]
o_16s_rural_int_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$tree_pit_type == "Rural Interior")]

# make matrices
leaf_its_pit <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_pit_names)]
leaf_its_grass <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_grass_names)]
leaf_its_urban_edge <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_urban_edge_names)]
leaf_its_urban_int <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_urban_int_names)]
leaf_its_rural_edge <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_rural_edge_names,)]
leaf_its_rural_int <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_rural_int_names)]

root_its_pit <- root_its[,which(colnames(root_its) %in% root_its_pit_names)]
root_its_grass <- root_its[,which(colnames(root_its) %in% root_its_grass_names)]
root_its_urban_edge <- root_its[,which(colnames(root_its) %in% root_its_urban_edge_names)]
root_its_urban_int <- root_its[,which(colnames(root_its) %in% root_its_urban_int_names)]
root_its_rural_edge <- root_its[,which(colnames(root_its) %in% root_its_rural_edge_names,)]
root_its_rural_int <- root_its[,which(colnames(root_its) %in% root_its_rural_int_names)]

m_its_pit <- m_its[,which(colnames(m_its) %in% m_its_pit_names)]
m_its_grass <- m_its[,which(colnames(m_its) %in% m_its_grass_names)]
m_its_urban_edge <- m_its[,which(colnames(m_its) %in% m_its_urban_edge_names)]
m_its_urban_int <- m_its[,which(colnames(m_its) %in% m_its_urban_int_names)]
m_its_rural_edge <- m_its[,which(colnames(m_its) %in% m_its_rural_edge_names,)]
m_its_rural_int <- m_its[,which(colnames(m_its) %in% m_its_rural_int_names)]

o_its_pit <- o_its[,which(colnames(o_its) %in% o_its_pit_names)]
o_its_grass <- o_its[,which(colnames(o_its) %in% o_its_grass_names)]
o_its_urban_edge <- o_its[,which(colnames(o_its) %in% o_its_urban_edge_names)]
o_its_urban_int <- o_its[,which(colnames(o_its) %in% o_its_urban_int_names)]
o_its_rural_edge <- o_its[,which(colnames(o_its) %in% o_its_rural_edge_names,)]
o_its_rural_int <- o_its[,which(colnames(o_its) %in% o_its_rural_int_names)]

leaf_16s_pit <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_pit_names)]
leaf_16s_grass <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_grass_names)]
leaf_16s_urban_edge <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_urban_edge_names)]
leaf_16s_urban_int <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_urban_int_names)]
leaf_16s_rural_edge <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_rural_edge_names,)]
leaf_16s_rural_int <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_rural_int_names)]

root_16s_pit <- root_16s[,which(colnames(root_16s) %in% root_16s_pit_names)]
root_16s_grass <- root_16s[,which(colnames(root_16s) %in% root_16s_grass_names)]
root_16s_urban_edge <- root_16s[,which(colnames(root_16s) %in% root_16s_urban_edge_names)]
root_16s_urban_int <- root_16s[,which(colnames(root_16s) %in% root_16s_urban_int_names)]
root_16s_rural_edge <- root_16s[,which(colnames(root_16s) %in% root_16s_rural_edge_names,)]
root_16s_rural_int <- root_16s[,which(colnames(root_16s) %in% root_16s_rural_int_names)]

m_16s_pit <- m_16s[,which(colnames(m_16s) %in% m_16s_pit_names)]
m_16s_grass <- m_16s[,which(colnames(m_16s) %in% m_16s_grass_names)]
m_16s_urban_edge <- m_16s[,which(colnames(m_16s) %in% m_16s_urban_edge_names)]
m_16s_urban_int <- m_16s[,which(colnames(m_16s) %in% m_16s_urban_int_names)]
m_16s_rural_edge <- m_16s[,which(colnames(m_16s) %in% m_16s_rural_edge_names,)]
m_16s_rural_int <- m_16s[,which(colnames(m_16s) %in% m_16s_rural_int_names)]

o_16s_pit <- o_16s[,which(colnames(o_16s) %in% o_16s_pit_names)]
o_16s_grass <- o_16s[,which(colnames(o_16s) %in% o_16s_grass_names)]
o_16s_urban_edge <- o_16s[,which(colnames(o_16s) %in% o_16s_urban_edge_names)]
o_16s_urban_int <- o_16s[,which(colnames(o_16s) %in% o_16s_urban_int_names)]
o_16s_rural_edge <- o_16s[,which(colnames(o_16s) %in% o_16s_rural_edge_names,)]
o_16s_rural_int <- o_16s[,which(colnames(o_16s) %in% o_16s_rural_int_names)]

# make lists
leaf_its <- list(leaf_its_pit, leaf_its_grass, leaf_its_urban_edge, 
                 leaf_its_urban_int, leaf_its_rural_edge, leaf_its_rural_int)

root_its <- list(root_its_pit, root_its_grass, root_its_urban_edge, 
                 root_its_urban_int, root_its_rural_edge, root_its_rural_int)

m_its <- list(m_its_pit, m_its_grass, m_its_urban_edge, 
                 m_its_urban_int, m_its_rural_edge, m_its_rural_int)

o_its <- list(o_its_pit, o_its_grass, o_its_urban_edge, 
                 o_its_urban_int, o_its_rural_edge, o_its_rural_int)

leaf_16s <- list(leaf_16s_pit, leaf_16s_grass, leaf_16s_urban_edge, 
                 leaf_16s_urban_int, leaf_16s_rural_edge, leaf_16s_rural_int)

root_16s <- list(root_16s_pit, root_16s_grass, root_16s_urban_edge, 
                 root_16s_urban_int, root_16s_rural_edge, root_16s_rural_int)

m_16s <- list(m_16s_pit, m_16s_grass, m_16s_urban_edge, 
              m_16s_urban_int, m_16s_rural_edge, m_16s_rural_int)

o_16s <- list(o_16s_pit, o_16s_grass, o_16s_urban_edge, 
              o_16s_urban_int, o_16s_rural_edge, o_16s_rural_int)

# filter out genera with < 25 reads
for(i in 1:length(leaf_its)){
  x <- leaf_its[[i]]
  x <- x[which(rowSums(x) > 25),]
  leaf_its[[i]] <- x
}

for(i in 1:length(root_its)){
  x <- root_its[[i]]
  x <- x[which(rowSums(x) > 25),]
  root_its[[i]] <- x
}

for(i in 1:length(m_its)){
  x <- m_its[[i]]
  x <- x[which(rowSums(x) > 25),]
  m_its[[i]] <- x
}

for(i in 1:length(o_its)){
  x <- o_its[[i]]
  x <- x[which(rowSums(x) > 25),]
  o_its[[i]] <- x
}

for(i in 1:length(leaf_16s)){
  x <- leaf_16s[[i]]
  x <- x[which(rowSums(x) > 25),]
  leaf_16s[[i]] <- x
}

for(i in 1:length(root_16s)){
  x <- root_16s[[i]]
  x <- x[which(rowSums(x) > 25),]
  root_16s[[i]] <- x
}

for(i in 1:length(m_16s)){
  x <- m_16s[[i]]
  x <- x[which(rowSums(x) > 25),]
  m_16s[[i]] <- x
}

for(i in 1:length(o_16s)){
  x <- o_16s[[i]]
  x <- x[which(rowSums(x) > 25),]
  o_16s[[i]] <- x
}

# filter out genera present in < 20% of samples in a network
for(i in 1:length(leaf_its)){
  x <- leaf_its[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  leaf_its[[i]] <- x
}

for(i in 1:length(root_its)){
  x <- root_its[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  root_its[[i]] <- x
}

for(i in 1:length(m_its)){
  x <- m_its[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  m_its[[i]] <- x
}

for(i in 1:length(o_its)){
  x <- o_its[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  o_its[[i]] <- x
}

for(i in 1:length(leaf_16s)){
  x <- leaf_16s[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  leaf_16s[[i]] <- x
}

for(i in 1:length(root_16s)){
  x <- root_16s[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  root_16s[[i]] <- x
}

for(i in 1:length(m_16s)){
  x <- m_16s[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  m_16s[[i]] <- x
}

for(i in 1:length(o_16s)){
  x <- o_16s[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  o_16s[[i]] <- x
}

# filter out genera < 1% abundance
for(i in 1:length(leaf_its)){
  x <- leaf_its[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  leaf_its[[i]] <- x
}

for(i in 1:length(root_its)){
  x <- root_its[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  root_its[[i]] <- x
}

for(i in 1:length(m_its)){
  x <- m_its[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  m_its[[i]] <- x
}

for(i in 1:length(o_its)){
  x <- o_its[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  o_its[[i]] <- x
}

for(i in 1:length(leaf_16s)){
  x <- leaf_16s[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  leaf_16s[[i]] <- x
}

for(i in 1:length(root_16s)){
  x <- root_16s[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  root_16s[[i]] <- x
}

for(i in 1:length(m_16s)){
  x <- m_16s[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  m_16s[[i]] <- x
}

for(i in 1:length(o_16s)){
  x <- o_16s[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  o_16s[[i]] <- x
}

# save tables
for(i in 1:length(leaf_its)){
  rownames(leaf_its[[i]]) <- gsub("g__", "f__", rownames(leaf_its[[i]]))
  colnames(leaf_its[[i]]) <- gsub("_ITS.*", "", colnames(leaf_its[[i]]))
}

for(i in 1:length(root_its)){
  rownames(root_its[[i]]) <- gsub("g__", "f__", rownames(root_its[[i]]))
  colnames(root_its[[i]]) <- gsub("_ITS.*", "", colnames(root_its[[i]]))
}

for(i in 1:length(m_its)){
  rownames(m_its[[i]]) <- gsub("g__", "f__", rownames(m_its[[i]]))
  colnames(m_its[[i]]) <- gsub("_ITS.*", "", colnames(m_its[[i]]))
}

for(i in 1:length(o_its)){
  rownames(o_its[[i]]) <- gsub("g__", "f__", rownames(o_its[[i]]))
  colnames(o_its[[i]]) <- gsub("_ITS.*", "", colnames(o_its[[i]]))
}

for(i in 1:length(leaf_16s)){
  rownames(leaf_16s[[i]]) <- gsub("g__", "b__", rownames(leaf_16s[[i]]))
  colnames(leaf_16s[[i]]) <- gsub("_16S.*", "", colnames(leaf_16s[[i]]))
}

for(i in 1:length(root_16s)){
  rownames(root_16s[[i]]) <- gsub("g__", "b__", rownames(root_16s[[i]]))
  colnames(root_16s[[i]]) <- gsub("_16S.*", "", colnames(root_16s[[i]]))
}

for(i in 1:length(m_16s)){
  rownames(m_16s[[i]]) <- gsub("g__", "b__", rownames(m_16s[[i]]))
  colnames(m_16s[[i]]) <- gsub("_16S.*", "", colnames(m_16s[[i]]))
}

for(i in 1:length(o_16s)){
  rownames(o_16s[[i]]) <- gsub("g__", "b__", rownames(o_16s[[i]]))
  colnames(o_16s[[i]]) <- gsub("_16S.*", "", colnames(o_16s[[i]]))
}

saveRDS(leaf_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/leaf_its_treetype.RDS")
saveRDS(root_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/root_its_treetype.RDS")
saveRDS(m_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/m_its_treetype.RDS")
saveRDS(o_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/o_its_treetype.RDS")

saveRDS(leaf_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/leaf_16s_treetype.RDS")
saveRDS(root_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/root_16s_treetype.RDS")
saveRDS(m_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/m_16s_treetype.RDS")
saveRDS(o_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/o_16s_treetype.RDS")