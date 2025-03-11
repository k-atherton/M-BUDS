library(ggplot2)
library(ggrepel)
library(vegan)
library(lmerTest)
library(tidyr)
library(MuMIn)
library(dplyr)
library(multcompView)
library(plotrix)
library(ggpubr)
library(Hmisc)
library(reshape2)
library(wesanderson)
library(patchwork)
library(indicspecies)
library(stringr)
library(Biostrings)
library(seqinr)
library(ggtree)
library(ape)
library(ggnewscale)

# read in metadata
metadata_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Data/metadata_its.csv")

# read in ASV data
all_clr_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/CLR/ITS/all_its_asv_clr_concat.csv", row.names = 1)
# all_rare <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/all_its_asv_rarefied.csv", row.names = 1)
leaf_rare_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Rarefaction/ITS/leaf_its_asv_rare_wtax_guild.csv", row.names = 1)
m_rare_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Rarefaction/ITS/m_its_asv_rare_wtax_guild.csv", row.names = 1)
o_rare_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Rarefaction/ITS/o_its_asv_rare_wtax_guild.csv", row.names = 1)
root_rare_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Rarefaction/ITS/root_its_asv_rare_wtax_guild.csv", row.names = 1)

leaf_rare_its <- leaf_rare_its[,10:ncol(leaf_rare_its)]
root_rare_its <- root_rare_its[,10:ncol(root_rare_its)]
m_rare_its <- m_rare_its[,10:ncol(m_rare_its)]
o_rare_its <- o_rare_its[,10:ncol(o_rare_its)]

leaf_rare_its <- aggregate(leaf_rare_its, . ~ Species, FUN = "sum")
root_rare_its <- aggregate(root_rare_its, . ~ Species, FUN = "sum")
m_rare_its <- aggregate(m_rare_its, . ~ Species, FUN = "sum")
o_rare_its <- aggregate(o_rare_its, . ~ Species, FUN = "sum")

bg_rare_its <- merge(root_rare_its, m_rare_its, by = "Species", all = T)
bg_rare_its <- merge(bg_rare_its, o_rare_its, by = "Species", all = T)

metadata_its$sample_type<- paste0(metadata_its$soil_horizon, " ", metadata_its$sample_type)
metadata_its$sample_type <- gsub("NA ", "", metadata_its$sample_type)

metadata_its$sample_type <- gsub("Soil", "Belowground", metadata_its$sample_type)
metadata_its$sample_type <- gsub("Root", "Belowground", metadata_its$sample_type)


metadata_its$sample_type <- factor(metadata_its$sample_type, levels = c("Leaf", "Belowground"))


metadata_its$tree_pit_type <- gsub("Pit", "Street Tree", metadata_its$tree_pit_type)
metadata_its$tree_pit_type <- gsub("Grass", "Street Tree", metadata_its$tree_pit_type)


metadata_its <- metadata_its[which(metadata_its$sample_name %in% colnames(all_clr_its)),]
metadata_its$tree_pit_type <- factor(metadata_its$tree_pit_type, 
                                     levels = c("Street Tree", "Urban Edge", 
                                                "Urban Interior", "Rural Edge", 
                                                "Rural Interior"))

fungal_trait_names <- c("perc_algal_parasites", "perc_animal_parasites", "perc_amf", "perc_arthropod_associated", "perc_dung_saprotrophs", "perc_ecm", "perc_epiphytes", "perc_foliar_endophytes", "perc_lichen_parasites", "perc_lichenized", "perc_litter_saprotrophs", "perc_moss_symbionts", "perc_mycoparasites", "perc_nectar_sap_saprotrophs", "perc_plant_pathogenic_capacity", "perc_pollen_saprotrophs", "perc_root_endophytes", "perc_soil_saprotrophs", "perc_sooty_molds", "perc_wood_saprotrophs", "perc_animal_associated_pathogen", "perc_human_pathogen", "perc_opportunistic_pathogen", "perc_root_associated_pathogen", "perc_leaf_fruit_seed_pathogen", "perc_wood_pathogen", "shannon_diversity")

# format data
metadata_its_leaf <- metadata_its[which(metadata_its$sample_type == "Leaf"),]
metadata_its_bg <- metadata_its[which(metadata_its$sample_type == "Belowground"),]

asv_leaf_its <- leaf_rare_its[,11:(ncol(leaf_rare_its)-1)]
asv_o_its <- o_rare_its[,11:(ncol(o_rare_its)-1)]
asv_m_its <- m_rare_its[,11:(ncol(m_rare_its)-1)]
asv_root_its <- root_rare_its[,11:(ncol(root_rare_its)-1)]

rownames(asv_leaf_its) <- leaf_rare_its$asv
rownames(asv_o_its) <- o_rare_its$asv
rownames(asv_m_its) <- m_rare_its$asv
rownames(asv_root_its) <- root_rare_its$asv

asv_leaf.t_its <- t(asv_leaf_its)
asv_o.t_its <- t(asv_o_its)
asv_m.t_its <- t(asv_m_its)
asv_root.t_its <- t(asv_root_its)

asv_leaf.t_its <- asv_leaf.t_its[order(rownames(asv_leaf.t_its)),]
asv_m.t_its <- asv_m.t_its[order(rownames(asv_m.t_its)),]
asv_o.t_its <- asv_o.t_its[order(rownames(asv_o.t_its)),]
asv_root.t_its <- asv_root.t_its[order(rownames(asv_root.t_its)),]

percent_leaf_its <- asv_leaf_its/mean(colSums(asv_leaf_its)) * 100
percent_o_its <- asv_leaf_its/mean(colSums(asv_o_its)) * 100
percent_m_its <- asv_leaf_its/mean(colSums(asv_m_its)) * 100
percent_root_its <- asv_leaf_its/mean(colSums(asv_root_its)) * 100

# run indicspecies
inv_leaf_its <- multipatt(asv_leaf.t_its, cluster=factor(metadata_its_leaf$tree_pit_type), func = "r.g", duleg=TRUE)
inv_o_its <- multipatt(asv_o.t_its, cluster=factor(metadata_its_o$tree_pit_type), func = "r.g", duleg=TRUE)
inv_m_its <- multipatt(asv_m.t_its, cluster=factor(metadata_its_m$tree_pit_type), func = "r.g", duleg=TRUE)
inv_root_its <- multipatt(asv_root.t_its, cluster=factor(metadata_its_root$tree_pit_type), func = "r.g", duleg=TRUE)

# save information
options(max.print=100000)
sink("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_leaf_its.csv")
summary(inv_leaf_its)
sink()

sink("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_o_its.csv")
summary(inv_o_its)
sink()

sink("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_m_its.csv")
summary(inv_m_its)
sink()

sink("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_root_its.csv")
summary(inv_root_its)
sink()

str_leaf_its <- inv_leaf_its$str
write.csv(str_leaf_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/correlation.inv_leaf_its.csv")

str_o_its <- inv_o_its$str
write.csv(str_o_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/correlation.inv_o_its.csv")

str_m_its <- inv_m_its$str
write.csv(str_m_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/correlation.inv_m_its.csv")

str_root_its <- inv_root_its$str
write.csv(str_root_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/correlation.inv_root_its.csv")

# get indicator asvs of each location group
summary_leaf_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_leaf_its.csv")
rows_leaf_its <- grep(" Group", summary_leaf_its$X..summary.inv_leaf_its.)

summary_o_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_o_its.csv")
rows_o_its <- grep(" Group", summary_o_its$X..summary.inv_o_its.)

summary_m_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_m_its.csv")
rows_m_its <- grep(" Group", summary_m_its$X..summary.inv_m_its.)

summary_root_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_root_its.csv")
rows_root_its <- grep(" Group", summary_root_its$X..summary.inv_root_its.)

leaf_street_tree_its <- data.frame(summary_leaf_its[(rows_leaf_its[1]+1):(rows_leaf_its[2]-1),])
leaf_urban_edge_its <- data.frame(summary_leaf_its[(rows_leaf_its[2]+1):(rows_leaf_its[3]-1),])
leaf_urban_interior_its <- data.frame(summary_leaf_its[(rows_leaf_its[3]+1):(rows_leaf_its[4]-1),])
leaf_rural_edge_its <- data.frame(summary_leaf_its[(rows_leaf_its[4]+1):(rows_leaf_its[5]-1),])
leaf_rural_interior_its <- data.frame(summary_leaf_its[(rows_leaf_its[5]+1):(nrow(summary_leaf_its)-2),])

o_street_tree_its <- data.frame(summary_o_its[(rows_o_its[1]+1):(rows_o_its[2]-1),])
o_urban_edge_its <- data.frame(summary_o_its[(rows_o_its[2]+1):(rows_o_its[3]-1),])
o_urban_interior_its <- data.frame(summary_o_its[(rows_o_its[3]+1):(rows_o_its[4]-1),])
o_rural_edge_its <- data.frame(summary_o_its[(rows_o_its[4]+1):(rows_o_its[5]-1),])
o_rural_interior_its <- data.frame(summary_o_its[(rows_o_its[5]+1):(nrow(summary_o_its)-2),])

m_street_tree_its <- data.frame(summary_m_its[(rows_m_its[1]+1):(rows_m_its[2]-1),])
m_urban_edge_its <- data.frame(summary_m_its[(rows_m_its[2]+1):(rows_m_its[3]-1),])
m_urban_interior_its <- data.frame(summary_m_its[(rows_m_its[3]+1):(rows_m_its[4]-1),])
m_rural_edge_its <- data.frame(summary_m_its[(rows_m_its[4]+1):(rows_m_its[5]-1),])
m_rural_interior_its <- data.frame(summary_m_its[(rows_m_its[5]+1):(nrow(summary_m_its)-2),])

root_street_tree_its <- data.frame(summary_root_its[(rows_root_its[1]+1):(rows_root_its[2]-1),])
root_urban_edge_its <- data.frame(summary_root_its[(rows_root_its[2]+1):(rows_root_its[3]-1),])
root_urban_interior_its <- data.frame(summary_root_its[(rows_root_its[3]+1):(rows_root_its[4]-1),])
root_rural_edge_its <- data.frame(summary_root_its[(rows_root_its[4]+1):(rows_root_its[5]-1),])
root_rural_interior_its <- data.frame(summary_root_its[(rows_root_its[5]+1):(nrow(summary_root_its)-2),])

name_its <- c("Street.Tree", "Urban.Edge", "Urban.Interior", "Rural.Edge", "Rural.Interior")

sample.list_leaf_its <- list()
sample.list_leaf_its[[1]] <- leaf_street_tree_its
sample.list_leaf_its[[2]] <- leaf_urban_edge_its
sample.list_leaf_its[[3]] <- leaf_urban_interior_its
sample.list_leaf_its[[4]] <- leaf_rural_edge_its
sample.list_leaf_its[[5]] <- leaf_street_tree_its

sample.list_o_its <- list()
sample.list_o_its[[1]] <- o_street_tree_its
sample.list_o_its[[2]] <- o_urban_edge_its
sample.list_o_its[[3]] <- o_urban_interior_its
sample.list_o_its[[4]] <- o_rural_edge_its
sample.list_o_its[[5]] <- o_street_tree_its

sample.list_m_its <- list()
sample.list_m_its[[1]] <- m_street_tree_its
sample.list_m_its[[2]] <- m_urban_edge_its
sample.list_m_its[[3]] <- m_urban_interior_its
sample.list_m_its[[4]] <- m_rural_edge_its
sample.list_m_its[[5]] <- m_street_tree_its

sample.list_root_its <- list()
sample.list_root_its[[1]] <- root_street_tree_its
sample.list_root_its[[2]] <- root_urban_edge_its
sample.list_root_its[[3]] <- root_urban_interior_its
sample.list_root_its[[4]] <- root_rural_edge_its
sample.list_root_its[[5]] <- root_street_tree_its

# save
for(i in 1:5){
  sample.list_leaf_its[[i]] <- str_split_fixed(sample.list_leaf_its[[i]][2:nrow(sample.list_leaf_its[[i]]),], " 0.",3)
  # rownames(sample.list_leaf_its[[i]]) <- sample.list_leaf_its[[i]][,1]
  # rownames(sample.list_leaf_its[[i]]) <- gsub(" ", "", rownames(sample.list_leaf_its[[i]]))
  # sample.list_leaf_its[[i]][,1] <- rownames(sample.list_leaf_its[[i]])
  colnames(sample.list_leaf_its[[i]]) <- c("ASV_ID", "stat", "p.value")
  write.csv(sample.list_leaf_its[[i]], paste0("/projectnb/talbot-lab-data/Katies_data/Street_Trees/",name_its[i],"_indic_leaf_its.csv"), row.names = F)
}

for(i in 1:5){
  sample.list_o_its[[i]] <- str_split_fixed(sample.list_o_its[[i]][2:nrow(sample.list_o_its[[i]]),], " 0.",3)
  # rownames(sample.list_o_its[[i]]) <- sample.list_o_its[[i]][,1]
  # rownames(sample.list_o_its[[i]]) <- gsub(" ", "", rownames(sample.list_o_its[[i]]))
  # sample.list_o_its[[i]][,1] <- rownames(sample.list_o_its[[i]])
  colnames(sample.list_o_its[[i]]) <- c("ASV_ID", "stat", "p.value")
  write.csv(sample.list_o_its[[i]], paste0("/projectnb/talbot-lab-data/Katies_data/Street_Trees/",name_its[i],"_indic_o_its.csv"), row.names = F)
}

for(i in 1:5){
  sample.list_m_its[[i]] <- str_split_fixed(sample.list_m_its[[i]][2:nrow(sample.list_m_its[[i]]),], " 0.",3)
  # rownames(sample.list_m_its[[i]]) <- sample.list_m_its[[i]][,1]
  # rownames(sample.list_m_its[[i]]) <- gsub(" ", "", rownames(sample.list_m_its[[i]]))
  # sample.list_m_its[[i]][,1] <- rownames(sample.list_m_its[[i]])
  colnames(sample.list_m_its[[i]]) <- c("ASV_ID", "stat", "p.value")
  write.csv(sample.list_m_its[[i]], paste0("/projectnb/talbot-lab-data/Katies_data/Street_Trees/",name_its[i],"_indic_m_its.csv"), row.names = F)
}

for(i in 1:5){
  sample.list_root_its[[i]] <- str_split_fixed(sample.list_root_its[[i]][2:nrow(sample.list_root_its[[i]]),], " 0.",3)
  # rownames(sample.list_root_its[[i]]) <- sample.list_root_its[[i]][,1]
  # rownames(sample.list_root_its[[i]]) <- gsub(" ", "", rownames(sample.list_root_its[[i]]))
  # sample.list_root_its[[i]][,1] <- rownames(sample.list_root_its[[i]])
  colnames(sample.list_root_its[[i]]) <- c("ASV_ID", "stat", "p.value")
  write.csv(sample.list_root_its[[i]], paste0("/projectnb/talbot-lab-data/Katies_data/Street_Trees/",name_its[i],"_indic_root_its.csv"), row.names = F)
}

cols_in_database_its <- c("Kingdom", "Phylum", "Class", "Order", "Family", "GENUS", "Species", "asv", "primary_lifestyle")
tax_leaf_its <- leaf_rare_its[, colnames(leaf_rare_its) %in% cols_in_database_its]
tax_o_its <- o_rare_its[, colnames(o_rare_its) %in% cols_in_database_its]
tax_m_its <- m_rare_its[, colnames(m_rare_its) %in% cols_in_database_its]
tax_root_its <- root_rare_its[, colnames(root_rare_its) %in% cols_in_database_its]

leaf_street_tree_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Street.Tree_indic_leaf_its.csv", header = T)
leaf_urban_edge_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Edge_indic_leaf_its.csv", header = T)
leaf_urban_interior_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Interior_indic_leaf_its.csv", header = T)
leaf_rural_edge_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Edge_indic_leaf_its.csv", header = T)
leaf_rural_interior_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Interior_indic_leaf_its.csv", header = T)

o_street_tree_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Street.Tree_indic_o_its.csv", header = T)
o_urban_edge_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Edge_indic_o_its.csv", header = T)
o_urban_interior_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Interior_indic_o_its.csv", header = T)
o_rural_edge_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Edge_indic_o_its.csv", header = T)
o_rural_interior_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Interior_indic_o_its.csv", header = T)

m_street_tree_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Street.Tree_indic_m_its.csv", header = T)
m_urban_edge_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Edge_indic_m_its.csv", header = T)
m_urban_interior_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Interior_indic_m_its.csv", header = T)
m_rural_edge_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Edge_indic_m_its.csv", header = T)
m_rural_interior_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Interior_indic_m_its.csv", header = T)

root_street_tree_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Street.Tree_indic_root_its.csv", header = T)
root_urban_edge_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Edge_indic_root_its.csv", header = T)
root_urban_interior_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Interior_indic_root_its.csv", header = T)
root_rural_edge_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Edge_indic_root_its.csv", header = T)
root_rural_interior_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Interior_indic_root_its.csv", header = T)

leaf_street_tree_its_merged <- merge(leaf_street_tree_its, tax_leaf_its, by.x = "ASV_ID", by.y = "asv")
leaf_urban_edge_its_merged <- merge(leaf_urban_edge_its, tax_leaf_its, by.x = "ASV_ID", by.y = "asv")
leaf_urban_interior_its_merged <- merge(leaf_urban_interior_its, tax_leaf_its, by.x = "ASV_ID", by.y = "asv")
leaf_rural_edge_its_merged <- merge(leaf_rural_edge_its, tax_leaf_its, by.x = "ASV_ID", by.y = "asv")
leaf_rural_interior_its_merged <- merge(leaf_rural_interior_its, tax_leaf_its, by.x = "ASV_ID", by.y = "asv")

o_street_tree_its_merged <- merge(o_street_tree_its, tax_o_its, by.x = "ASV_ID", by.y = "asv")
o_urban_edge_its_merged <- merge(o_urban_edge_its, tax_o_its, by.x = "ASV_ID", by.y = "asv")
o_urban_interior_its_merged <- merge(o_urban_interior_its, tax_o_its, by.x = "ASV_ID", by.y = "asv")
o_rural_edge_its_merged <- merge(o_rural_edge_its, tax_o_its, by.x = "ASV_ID", by.y = "asv")
o_rural_interior_its_merged <- merge(o_rural_interior_its, tax_o_its, by.x = "ASV_ID", by.y = "asv")

m_street_tree_its_merged <- merge(m_street_tree_its, tax_m_its, by.x = "ASV_ID", by.y = "asv")
m_urban_edge_its_merged <- merge(m_urban_edge_its, tax_m_its, by.x = "ASV_ID", by.y = "asv")
m_urban_interior_its_merged <- merge(m_urban_interior_its, tax_m_its, by.x = "ASV_ID", by.y = "asv")
m_rural_edge_its_merged <- merge(m_rural_edge_its, tax_m_its, by.x = "ASV_ID", by.y = "asv")
m_rural_interior_its_merged <- merge(m_rural_interior_its, tax_m_its, by.x = "ASV_ID", by.y = "asv")

root_street_tree_its_merged <- merge(root_street_tree_its, tax_root_its, by.x = "ASV_ID", by.y = "asv")
root_urban_edge_its_merged <- merge(root_urban_edge_its, tax_root_its, by.x = "ASV_ID", by.y = "asv")
root_urban_interior_its_merged <- merge(root_urban_interior_its, tax_root_its, by.x = "ASV_ID", by.y = "asv")
root_rural_edge_its_merged <- merge(root_rural_edge_its, tax_root_its, by.x = "ASV_ID", by.y = "asv")
root_rural_interior_its_merged <- merge(root_rural_interior_its, tax_root_its, by.x = "ASV_ID", by.y = "asv")

rownames(leaf_street_tree_its_merged) <- leaf_street_tree_its_merged[,1]
rownames(leaf_urban_edge_its_merged) <- leaf_urban_edge_its_merged[,1]
rownames(leaf_urban_interior_its_merged) <- leaf_urban_interior_its_merged[,1]
rownames(leaf_rural_edge_its_merged) <- leaf_rural_edge_its_merged[,1]
rownames(leaf_rural_interior_its_merged) <- leaf_rural_interior_its_merged[,1]

rownames(o_street_tree_its_merged) <- o_street_tree_its_merged[,1]
rownames(o_urban_edge_its_merged) <- o_urban_edge_its_merged[,1]
rownames(o_urban_interior_its_merged) <- o_urban_interior_its_merged[,1]
rownames(o_rural_edge_its_merged) <- o_rural_edge_its_merged[,1]
rownames(o_rural_interior_its_merged) <- o_rural_interior_its_merged[,1]

rownames(m_street_tree_its_merged) <- m_street_tree_its_merged[,1]
rownames(m_urban_edge_its_merged) <- m_urban_edge_its_merged[,1]
rownames(m_urban_interior_its_merged) <- m_urban_interior_its_merged[,1]
rownames(m_rural_edge_its_merged) <- m_rural_edge_its_merged[,1]
rownames(m_rural_interior_its_merged) <- m_rural_interior_its_merged[,1]

rownames(root_street_tree_its_merged) <- root_street_tree_its_merged[,1]
rownames(root_urban_edge_its_merged) <- root_urban_edge_its_merged[,1]
rownames(root_urban_interior_its_merged) <- root_urban_interior_its_merged[,1]
rownames(root_rural_edge_its_merged) <- root_rural_edge_its_merged[,1]
rownames(root_rural_interior_its_merged) <- root_rural_interior_its_merged[,1]

leaf_street_tree_its_merged$Group <- "Street.tree"
leaf_urban_edge_its_merged$Group <- "Urban.edge"
leaf_urban_interior_its_merged$Group <- "Urban.interior"
leaf_rural_edge_its_merged$Group <- "Rural.edge"
leaf_rural_interior_its_merged$Group <- "Rural.interior"

o_street_tree_its_merged$Group <- "Street.tree"
o_urban_edge_its_merged$Group <- "Urban.edge"
o_urban_interior_its_merged$Group <- "Urban.interior"
o_rural_edge_its_merged$Group <- "Rural.edge"
o_rural_interior_its_merged$Group <- "Rural.interior"

m_street_tree_its_merged$Group <- "Street.tree"
m_urban_edge_its_merged$Group <- "Urban.edge"
m_urban_interior_its_merged$Group <- "Urban.interior"
m_rural_edge_its_merged$Group <- "Rural.edge"
m_rural_interior_its_merged$Group <- "Rural.interior"

root_street_tree_its_merged$Group <- "Street.tree"
root_urban_edge_its_merged$Group <- "Urban.edge"
root_urban_interior_its_merged$Group <- "Urban.interior"
root_rural_edge_its_merged$Group <- "Rural.edge"
root_rural_interior_its_merged$Group <- "Rural.interior"

leaf_street_tree_its_merged$Sample.Type <- "Leaf"
leaf_urban_edge_its_merged$Sample.Type <- "Leaf"
leaf_urban_interior_its_merged$Sample.Type <- "Leaf"
leaf_rural_edge_its_merged$Sample.Type <- "Leaf"
leaf_rural_interior_its_merged$Sample.Type <- "Leaf"

o_street_tree_its_merged$Sample.Type <- "O"
o_urban_edge_its_merged$Sample.Type <- "O"
o_urban_interior_its_merged$Sample.Type <- "O"
o_rural_edge_its_merged$Sample.Type <- "O"
o_rural_interior_its_merged$Sample.Type <- "O"

m_street_tree_its_merged$Sample.Type <- "M"
m_urban_edge_its_merged$Sample.Type <- "M"
m_urban_interior_its_merged$Sample.Type <- "M"
m_rural_edge_its_merged$Sample.Type <- "M"
m_rural_interior_its_merged$Sample.Type <- "M"

root_street_tree_its_merged$Sample.Type <- "Root"
root_urban_edge_its_merged$Sample.Type <- "Root"
root_urban_interior_its_merged$Sample.Type <- "Root"
root_rural_edge_its_merged$Sample.Type <- "Root"
root_rural_interior_its_merged$Sample.Type <- "Root"

leaf_street_tree_its <- leaf_street_tree_its_merged[,-c(1:3)]
leaf_urban_edge_its <- leaf_urban_edge_its_merged[,-c(1:3)]
leaf_urban_interior_its <- leaf_urban_interior_its_merged[,-c(1:3)]
leaf_rural_edge_its <- leaf_rural_edge_its_merged[,-c(1:3)]
leaf_rural_interior_its <- leaf_rural_interior_its_merged[,-c(1:3)]

o_street_tree_its <- o_street_tree_its_merged[,-c(1:3)]
o_urban_edge_its <- o_urban_edge_its_merged[,-c(1:3)]
o_urban_interior_its <- o_urban_interior_its_merged[,-c(1:3)]
o_rural_edge_its <- o_rural_edge_its_merged[,-c(1:3)]
o_rural_interior_its <- o_rural_interior_its_merged[,-c(1:3)]

m_street_tree_its <- m_street_tree_its_merged[,-c(1:3)]
m_urban_edge_its <- m_urban_edge_its_merged[,-c(1:3)]
m_urban_interior_its <- m_urban_interior_its_merged[,-c(1:3)]
m_rural_edge_its <- m_rural_edge_its_merged[,-c(1:3)]
m_rural_interior_its <- m_rural_interior_its_merged[,-c(1:3)]

root_street_tree_its <- root_street_tree_its_merged[,-c(1:3)]
root_urban_edge_its <- root_urban_edge_its_merged[,-c(1:3)]
root_urban_interior_its <- root_urban_interior_its_merged[,-c(1:3)]
root_rural_edge_its <- root_rural_edge_its_merged[,-c(1:3)]
root_rural_interior_its <- root_rural_interior_its_merged[,-c(1:3)]

# save sequences
rep_seqs_its <- readDNAStringSet("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/its_representative_seqs_unique.fasta")
seq_name_its <- names(rep_seqs_its)
sequence_its <- paste(rep_seqs_its)
rep_seqs_its <- data.frame(seq_name_its, sequence_its)
all_leaf_its <- rbind(leaf_street_tree_its, leaf_urban_edge_its, leaf_urban_interior_its, leaf_rural_edge_its, leaf_rural_interior_its)
all_o_its <- rbind(o_street_tree_its, o_urban_edge_its, o_urban_interior_its, o_rural_edge_its, o_rural_interior_its)
all_m_its <- rbind(m_street_tree_its, m_urban_edge_its, m_urban_interior_its, m_rural_edge_its, m_rural_interior_its)
all_root_its <- rbind(root_street_tree_its, root_urban_edge_its, root_urban_interior_its, root_rural_edge_its, root_rural_interior_its)

write.csv(all_leaf_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_indic_phylogenetic_info.csv")
write.csv(all_o_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_indic_phylogenetic_info.csv")
write.csv(all_m_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_indic_phylogenetic_info.csv")
write.csv(all_root_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_indic_phylogenetic_info.csv")

all_leaf_fasta_its <- merge(all_leaf_its, rep_seqs_its, by.x = 0, by.y = "seq_name_its")
all_o_fasta_its <- merge(all_o_its, rep_seqs_its, by.x = 0, by.y = "seq_name_its")
all_m_fasta_its <- merge(all_m_its, rep_seqs_its, by.x = 0, by.y = "seq_name_its")
all_root_fasta_its <- merge(all_root_its, rep_seqs_its, by.x = 0, by.y = "seq_name_its")

write.fasta(sequences = as.list(all_leaf_fasta_its$sequence_its), names = all_leaf_fasta_its$Row.names, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_indic_seqs.fasta")
write.fasta(sequences = as.list(all_o_fasta_its$sequence_its), names = all_o_fasta_its$Row.names, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_indic_seqs.fasta")
write.fasta(sequences = as.list(all_m_fasta_its$sequence_its), names = all_m_fasta_its$Row.names, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_indic_seqs.fasta")
write.fasta(sequences = as.list(all_root_fasta_its$sequence_its), names = all_root_fasta_its$Row.names, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_indic_seqs.fasta")

# get 28S sequence
sample_leaf_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_indic_phylogenetic_info.csv", header = T)
sample_o_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_indic_phylogenetic_info.csv", header = T)
sample_m_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_indic_phylogenetic_info.csv", header = T)
sample_root_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_indic_phylogenetic_info.csv", header = T)

sample_leaf_its <- data.frame(lapply(sample_leaf_its, function(x){gsub(pattern="o__", replacement = "", x)}), stringsAsFactors = FALSE)
sample_o_its <- data.frame(lapply(sample_o_its, function(x){gsub(pattern="o__", replacement = "", x)}), stringsAsFactors = FALSE)
sample_m_its <- data.frame(lapply(sample_m_its, function(x){gsub(pattern="o__", replacement = "", x)}), stringsAsFactors = FALSE)
sample_root_its <- data.frame(lapply(sample_root_its, function(x){gsub(pattern="o__", replacement = "", x)}), stringsAsFactors = FALSE)

colnames(sample_leaf_its)[1] <- "ASV_ID"
colnames(sample_o_its)[1] <- "ASV_ID"
colnames(sample_m_its)[1] <- "ASV_ID"
colnames(sample_root_its)[1] <- "ASV_ID"

database_its <- read.fasta("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/SILVA_138.1_LSURef_NR99_tax_silva.fasta", as.string=TRUE)
at_its <- data.frame(attributes(database_its))

# format database
tax.list_its <- NULL
id.list_its <- NULL
seq.list_its <- NULL
for(i in 1:nrow(at_its)){
  description <- attributes(database_its[[at_its[i,]]])[["Annot"]]
  annot <- strsplit(description, " +")[[1]][2]
  if(grepl("Fungi", annot, fixed = TRUE)){
    id <- at_its[i,]
    id.list_its <- c(id.list_its, id)
    seq <- database_its[[at_its[i,]]]
    seq.list_its <- rbind(seq.list_its, seq)
    tax.name <- strsplit(annot, ";")[[1]]
    tax.name <- c(tax.name, rep(NA, (14-length(tax.name))))
    tax.list_its <- rbind(tax.list_its, tax.name)}
}
tax.data_its <- data.frame(tax.list_its)
seq.data_its <- data.frame(seq.list_its)

database_its <- tax.data_its[,6:14]
database_its$ID <- id.list_its
database_its$seq <- seq.data_its$seq.list_its
database_its <- data.frame(database_its)

colnames(database_its) <- c("Kingdom","Subkingdom","Phylum","Subphylum","Class","Order","Family","Genus","Species","ID","Sequence")

# assign 28S sequences to ASVs referring to Genus
database_its.genus <- na.omit(cbind(database_its$Genus, database_its$Sequence))
colnames(database_its.genus) <- c("Genus", "Sequence")
database_its.genus <- data.frame(database_its.genus)
db_its.genus <- database_its.genus[!duplicated(database_its.genus$Genus),]

merge_leaf_its.genus <- merge(sample_leaf_its, db_its.genus, by.x = "GENUS", by.y = "Genus", all.x = TRUE)
merge_leaf_its.genus <- merge_leaf_its.genus[order(merge_leaf_its.genus$ASV_ID, decreasing = FALSE),]
merge_o_its.genus <- merge(sample_o_its, db_its.genus, by.x = "GENUS", by.y = "Genus", all.x = TRUE)
merge_o_its.genus <- merge_o_its.genus[order(merge_o_its.genus$ASV_ID, decreasing = FALSE),]
merge_m_its.genus <- merge(sample_m_its, db_its.genus, by.x = "GENUS", by.y = "Genus", all.x = TRUE)
merge_m_its.genus <- merge_m_its.genus[order(merge_m_its.genus$ASV_ID, decreasing = FALSE),]
merge_root_its.genus <- merge(sample_root_its, db_its.genus, by.x = "GENUS", by.y = "Genus", all.x = TRUE)
merge_root_its.genus <- merge_root_its.genus[order(merge_root_its.genus$ASV_ID, decreasing = FALSE),]

# assign 28S sequences to ASVs referring to Family
database_its.family <- na.omit(cbind(database_its$Family, database_its$Sequence))
colnames(database_its.family) <- c("Family", "Sequence")
database_its.family <- data.frame(database_its.family)
db_its.family <- database_its.family[!duplicated(database_its.family$Family),]

merge_leaf_its.family <- merge(sample_leaf_its, db_its.family, by.x = "Family", by.y = "Family", all.x = TRUE)
merge_leaf_its.family <- merge_leaf_its.family[order(merge_leaf_its.family$ASV_ID, decreasing = FALSE),]
merge_o_its.family <- merge(sample_o_its, db_its.family, by.x = "Family", by.y = "Family", all.x = TRUE)
merge_o_its.family <- merge_o_its.family[order(merge_o_its.family$ASV_ID, decreasing = FALSE),]
merge_m_its.family <- merge(sample_m_its, db_its.family, by.x = "Family", by.y = "Family", all.x = TRUE)
merge_m_its.family <- merge_m_its.family[order(merge_m_its.family$ASV_ID, decreasing = FALSE),]
merge_root_its.family <- merge(sample_root_its, db_its.family, by.x = "Family", by.y = "Family", all.x = TRUE)
merge_root_its.family <- merge_root_its.family[order(merge_root_its.family$ASV_ID, decreasing = FALSE),]

# assign 28S sequences to ASVs referring to Order
database_its.order <- na.omit(cbind(database_its$Order, database_its$Sequence))
colnames(database_its.order) <- c("Order", "Sequence")
database_its.order <- data.frame(database_its.order)
db_its.order <- database_its.order[!duplicated(database_its.order$Order),]

merge_leaf_its.order <- merge(sample_leaf_its, db_its.order, by.x = "Order", by.y = "Order", all.x = TRUE)
merge_leaf_its.order <- merge_leaf_its.order[order(merge_leaf_its.order$ASV_ID, decreasing = FALSE),]
merge_o_its.order <- merge(sample_o_its, db_its.order, by.x = "Order", by.y = "Order", all.x = TRUE)
merge_o_its.order <- merge_o_its.order[order(merge_o_its.order$ASV_ID, decreasing = FALSE),]
merge_m_its.order <- merge(sample_m_its, db_its.order, by.x = "Order", by.y = "Order", all.x = TRUE)
merge_m_its.order <- merge_m_its.order[order(merge_m_its.order$ASV_ID, decreasing = FALSE),]
merge_root_its.order <- merge(sample_root_its, db_its.order, by.x = "Order", by.y = "Order", all.x = TRUE)
merge_root_its.order <- merge_root_its.order[order(merge_root_its.order$ASV_ID, decreasing = FALSE),]

# combine taxonomy information
data_leaf_its <- sample_leaf_its[order(sample_leaf_its$ASV_ID, decreasing = FALSE),]
data_leaf_its$seqs <- merge_leaf_its.genus$Sequence
data_o_its <- sample_o_its[order(sample_o_its$ASV_ID, decreasing = FALSE),]
data_o_its$seqs <- merge_o_its.genus$Sequence
data_m_its <- sample_m_its[order(sample_m_its$ASV_ID, decreasing = FALSE),]
data_m_its$seqs <- merge_m_its.genus$Sequence
data_root_its <- sample_root_its[order(sample_root_its$ASV_ID, decreasing = FALSE),]
data_root_its$seqs <- merge_root_its.genus$Sequence

for(i in 1:nrow(data_leaf_its)){
  if(is.na(data_leaf_its$seqs[i])){data_leaf_its$seqs[i] <- merge_leaf_its.family$Sequence[i]}
  if(is.na(data_leaf_its$seqs[i])){data_leaf_its$seqs[i] <- merge_leaf_its.order$Sequence[i]}}
for(i in 1:nrow(data_o_its)){
  if(is.na(data_o_its$seqs[i])){data_o_its$seqs[i] <- merge_o_its.family$Sequence[i]}
  if(is.na(data_o_its$seqs[i])){data_o_its$seqs[i] <- merge_o_its.order$Sequence[i]}}
for(i in 1:nrow(data_m_its)){
  if(is.na(data_m_its$seqs[i])){data_m_its$seqs[i] <- merge_m_its.family$Sequence[i]}
  if(is.na(data_m_its$seqs[i])){data_m_its$seqs[i] <- merge_m_its.order$Sequence[i]}}
for(i in 1:nrow(data_root_its)){
  if(is.na(data_root_its$seqs[i])){data_root_its$seqs[i] <- merge_root_its.family$Sequence[i]}
  if(is.na(data_root_its$seqs[i])){data_root_its$seqs[i] <- merge_root_its.order$Sequence[i]}}

data_leaf_its <- data_leaf_its[!is.na(data_leaf_its$seqs),]
data_o_its <- data_o_its[!is.na(data_o_its$seqs),]
data_m_its <- data_m_its[!is.na(data_m_its$seqs),]
data_root_its <- data_root_its[!is.na(data_root_its$seqs),]

rownames(data_leaf_its) <- data_leaf_its$ASV_ID
rownames(data_o_its) <- data_o_its$ASV_ID
rownames(data_m_its) <- data_m_its$ASV_ID
rownames(data_root_its) <- data_root_its$ASV_ID

write.csv(data_leaf_its[,-1], "/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_phylogenetic_info_28S.csv")
write.csv(data_o_its[,-1], "/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_phylogenetic_info_28S.csv")
write.csv(data_m_its[,-1], "/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_phylogenetic_info_28S.csv")
write.csv(data_root_its[,-1], "/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_phylogenetic_info_28S.csv")

write.fasta(sequences = as.list(data_leaf_its$seqs), names = data_leaf_its$ASV_ID, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_seqs_28s.fasta")
write.fasta(sequences = as.list(data_o_its$seqs), names = data_o_its$ASV_ID, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_seqs_28s.fasta")
write.fasta(sequences = as.list(data_m_its$seqs), names = data_m_its$ASV_ID, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_seqs_28s.fasta")
write.fasta(sequences = as.list(data_root_its$seqs), names = data_root_its$ASV_ID, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_seqs_28s.fasta")

# to visualize in iTOL
# make a binary file to show their presence or absence
data_leaf_its$Street.tree <- 0
data_leaf_its$Urban.edge <- 0
data_leaf_its$Urban.interior <- 0
data_leaf_its$Rural.edge <- 0
data_leaf_its$Rural.interior <- 0

data_o_its$Street.tree <- 0
data_o_its$Urban.edge <- 0
data_o_its$Urban.interior <- 0
data_o_its$Rural.edge <- 0
data_o_its$Rural.interior <- 0

data_m_its$Street.tree <- 0
data_m_its$Urban.edge <- 0
data_m_its$Urban.interior <- 0
data_m_its$Rural.edge <- 0
data_m_its$Rural.interior <- 0

data_root_its$Street.tree <- 0
data_root_its$Urban.edge <- 0
data_root_its$Urban.interior <- 0
data_root_its$Rural.edge <- 0
data_root_its$Rural.interior <- 0

for(i in 1:nrow(data_leaf_its)){
  if(data_leaf_its$Group[i] == "Street.tree"){data_leaf_its$Street.tree[i] <- 1}else{data_leaf_its$Street.tree[i] <- 0}
  if(data_leaf_its$Group[i] == "Urban.edge"){data_leaf_its$Urban.edge[i] <- 1}else{data_leaf_its$Urban.edge[i] <- 0}
  if(data_leaf_its$Group[i] == "Urban.interior"){data_leaf_its$Urban.interior[i] <- 1}else{data_leaf_its$Urban.interior[i] <- 0}
  if(data_leaf_its$Group[i] == "Rural.edge"){data_leaf_its$Rural.edge[i] <- 1}else{data_leaf_its$Rural.edge[i] <- 0}
  if(data_leaf_its$Group[i] == "Rural.interior"){data_leaf_its$Rural.interior[i] <- 1}else{data_leaf_its$Rural.interior[i] <- 0}
}

for(i in 1:nrow(data_o_its)){
  if(data_o_its$Group[i] == "Street.tree"){data_o_its$Street.tree[i] <- 1}else{data_o_its$Street.tree[i] <- 0}
  if(data_o_its$Group[i] == "Urban.edge"){data_o_its$Urban.edge[i] <- 1}else{data_o_its$Urban.edge[i] <- 0}
  if(data_o_its$Group[i] == "Urban.interior"){data_o_its$Urban.interior[i] <- 1}else{data_o_its$Urban.interior[i] <- 0}
  if(data_o_its$Group[i] == "Rural.edge"){data_o_its$Rural.edge[i] <- 1}else{data_o_its$Rural.edge[i] <- 0}
  if(data_o_its$Group[i] == "Rural.interior"){data_o_its$Rural.interior[i] <- 1}else{data_o_its$Rural.interior[i] <- 0}
}

for(i in 1:nrow(data_m_its)){
  if(data_m_its$Group[i] == "Street.tree"){data_m_its$Street.tree[i] <- 1}else{data_m_its$Street.tree[i] <- 0}
  if(data_m_its$Group[i] == "Urban.edge"){data_m_its$Urban.edge[i] <- 1}else{data_m_its$Urban.edge[i] <- 0}
  if(data_m_its$Group[i] == "Urban.interior"){data_m_its$Urban.interior[i] <- 1}else{data_m_its$Urban.interior[i] <- 0}
  if(data_m_its$Group[i] == "Rural.edge"){data_m_its$Rural.edge[i] <- 1}else{data_m_its$Rural.edge[i] <- 0}
  if(data_m_its$Group[i] == "Rural.interior"){data_m_its$Rural.interior[i] <- 1}else{data_m_its$Rural.interior[i] <- 0}
}

for(i in 1:nrow(data_root_its)){
  if(data_root_its$Group[i] == "Street.tree"){data_root_its$Street.tree[i] <- 1}else{data_root_its$Street.tree[i] <- 0}
  if(data_root_its$Group[i] == "Urban.edge"){data_root_its$Urban.edge[i] <- 1}else{data_root_its$Urban.edge[i] <- 0}
  if(data_root_its$Group[i] == "Urban.interior"){data_root_its$Urban.interior[i] <- 1}else{data_root_its$Urban.interior[i] <- 0}
  if(data_root_its$Group[i] == "Rural.edge"){data_root_its$Rural.edge[i] <- 1}else{data_root_its$Rural.edge[i] <- 0}
  if(data_root_its$Group[i] == "Rural.interior"){data_root_its$Rural.interior[i] <- 1}else{data_root_its$Rural.interior[i] <- 0}
}

heatmap_leaf_its <- cbind(data_leaf_its$ASV_ID, data_leaf_its$Street.tree, data_leaf_its$Urban.edge, data_leaf_its$Urban.interior, data_leaf_its$Rural.edge, data_leaf_its$Rural.interior)
heatmap_o_its <- cbind(data_o_its$ASV_ID, data_o_its$Street.tree, data_o_its$Urban.edge, data_o_its$Urban.interior, data_o_its$Rural.edge, data_o_its$Rural.interior)
heatmap_m_its <- cbind(data_m_its$ASV_ID, data_m_its$Street.tree, data_m_its$Urban.edge, data_m_its$Urban.interior, data_m_its$Rural.edge, data_m_its$Rural.interior)
heatmap_root_its <- cbind(data_root_its$ASV_ID, data_root_its$Street.tree, data_root_its$Urban.edge, data_root_its$Urban.interior, data_root_its$Rural.edge, data_root_its$Rural.interior)

write.csv(heatmap_leaf_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/heatmap_leaf_its_iTOL_28S.csv")
write.csv(heatmap_o_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/heatmap_o_its_iTOL_28S.csv")
write.csv(heatmap_m_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/heatmap_m_its_iTOL_28S.csv")
write.csv(heatmap_root_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/heatmap_root_its_iTOL_28S.csv")

# tree_its <- read.tree("Phylogenetic_tree_28S_Newick")
# tree_its
# 
# data_leaf_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_phylogenetic_info_28S.csv")
# data_o_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_phylogenetic_info_28S.csv")
# data_m_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_phylogenetic_info_28S.csv")
# data_root_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_phylogenetic_info_28S.csv")
# 
# data_leaf_its <- data.frame(lapply(data_leaf_its, function(x){gsub(pattern = "p__", replacement = "", x)}), stringsAsFactors = FALSE)
# data_o_its <- data.frame(lapply(data_o_its, function(x){gsub(pattern = "p__", replacement = "", x)}), stringsAsFactors = FALSE)
# data_m_its <- data.frame(lapply(data_m_its, function(x){gsub(pattern = "p__", replacement = "", x)}), stringsAsFactors = FALSE)
# data_root_its <- data.frame(lapply(data_root_its, function(x){gsub(pattern = "p__", replacement = "", x)}), stringsAsFactors = FALSE)
# 
# # to change color by phylum
# groupInfo_leaf_its <- split(data_leaf_its$X, data_leaf_its$Phylum)
# groupInfo_o_its <- split(data_o_its$X, data_o_its$Phylum)
# groupInfo_m_its <- split(data_m_its$X, data_m_its$Phylum)
# groupInfo_root_its <- split(data_root_its$X, data_root_its$Phylum)
# 
# tree_leaf_its <- groupOTU(tree_its, groupInfo_leaf_its, group_name = "Phylum")
# tree_o_its <- groupOTU(tree_its, groupInfo_o_its, group_name = "Phylum")
# tree_m_its <- groupOTU(tree_its, groupInfo_m_its, group_name = "Phylum")
# tree_root_its <- groupOTU(tree_its, groupInfo_root_its, group_name = "Phylum")
# 
# # save
# ggtree(tree_leaf_its, aes(color = Phylum), layout = 'circular')
# ggtree(tree_o_its, aes(color = Phylum), layout = 'circular')
# ggtree(tree_m_its, aes(color = Phylum), layout = 'circular')
# ggtree(tree_root_its, aes(color = Phylum), layout = 'circular')