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
metadata_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Data/metadata_16s_200.csv")

# read in ASV data
# all_aitch_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/all_16s_sample_aitch.csv", row.names = 1)
# all_clr_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/all_16s_asv_clr_wtax.csv", row.names = 1)
# # all_rare <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/all_16s_asv_rarefied.csv", row.names = 1)
leaf_rare_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Rarefaction/16S/leaf_16s_asv_rare_wtax_guild_200.csv", row.names = 1)
m_rare_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Rarefaction/16S/m_16s_asv_rare_wtax_guild.csv", row.names = 1)
o_rare_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Rarefaction/16S/o_16s_asv_rare_wtax_guild.csv", row.names = 1)
root_rare_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Rarefaction/16S/root_16s_asv_rare_wtax_guild.csv", row.names = 1)

colnames(leaf_rare_16s)[1] <- "asv"
colnames(root_rare_16s)[1] <- "asv"
colnames(m_rare_16s)[1] <- "asv"
colnames(o_rare_16s)[1] <- "asv"

metadata_16s$sample_type<- paste0(metadata_16s$soil_horizon, " ", metadata_16s$sample_type)
metadata_16s$sample_type <- gsub("NA ", "", metadata_16s$sample_type)

metadata_16s$sample_type <- gsub("O Soil", "Soil, 0-15 cm", metadata_16s$sample_type)
metadata_16s$sample_type <- gsub("M Soil", "Soil, 15-30 cm", metadata_16s$sample_type)

metadata_16s$sample_type <- gsub("O Root", "Roots, 0-15 cm", metadata_16s$sample_type)
metadata_16s$sample_type <- gsub("M Root", "Roots, 15-30 cm", metadata_16s$sample_type)

metadata_16s$sample_type <- factor(metadata_16s$sample_type, levels = c("Leaf", "Soil, 0-15 cm", "Soil, 15-30 cm", "Roots, 0-15 cm", "Roots, 15-30 cm"))


metadata_16s$tree_pit_type <- gsub("Pit", "Street Tree", metadata_16s$tree_pit_type)
metadata_16s$tree_pit_type <- gsub("Grass", "Street Tree", metadata_16s$tree_pit_type)

samples <- unique(c(colnames(leaf_rare_16s), colnames(root_rare_16s), colnames(m_rare_16s), colnames(o_rare_16s)))
metadata_16s <- metadata_16s[which(metadata_16s$sample_name %in% samples),]
metadata_16s$tree_pit_type <- factor(metadata_16s$tree_pit_type, 
                                     levels = c("Street Tree", "Urban Edge", 
                                                "Urban Interior", "Rural Edge", 
                                                "Rural Interior"))

bacterial_trait_names <- c("perc_c_fixation", "perc_cellulolytic", 
                           "perc_chitinolytic", "perc_lignolytic",
                           "perc_carbon_monoxide_oxidation", "perc_n_fixation",
                           "perc_copiotroph", "perc_denitrification",
                           "perc_dissim_nitrate_reduction", 
                           "perc_hydrocarbon_degradation",
                           "perc_sulfonate_desulfurization",
                           "perc_other_p_cycling",
                           "perc_oxidize_reduced_sulfur", "perc_methanotroph",
                           "perc_oligotroph", "perc_partial_nitrification", 
                           "perc_archaea", "perc_plant_pathogen", 
                           "perc_animal_pathogen", "perc_zoonotic_pathogen", 
                           "shannon", "shannon_nopath")
# format data
metadata_16s_leaf <- metadata_16s[which(metadata_16s$sample_type == "Leaf"),]
metadata_16s_o <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]
metadata_16s_m <- metadata_16s[which(metadata_16s$sample_type == "Soil, 15-30 cm"),]
metadata_16s_root <- metadata_16s[which(metadata_16s$sample_type %in% c("Roots, 0-15 cm", "Roots, 15-30 cm")),]

asv_leaf_16s <- leaf_rare_16s[,11:(ncol(leaf_rare_16s)-1)]
asv_o_16s <- o_rare_16s[,11:(ncol(o_rare_16s)-1)]
asv_m_16s <- m_rare_16s[,11:(ncol(m_rare_16s)-1)]
asv_root_16s <- root_rare_16s[,11:(ncol(root_rare_16s)-1)]

asv_leaf_16s <- asv_leaf_16s[,colnames(asv_leaf_16s) %in% metadata_16s_leaf$sample_name]
asv_root_16s <- asv_root_16s[,colnames(asv_root_16s) %in% metadata_16s_root$sample_name]
asv_m_16s <- asv_m_16s[,colnames(asv_m_16s) %in% metadata_16s_m$sample_name]
asv_o_16s <- asv_o_16s[,colnames(asv_o_16s) %in% metadata_16s_o$sample_name]

metadata_16s_leaf <- metadata_16s_leaf[which(metadata_16s_leaf$sample_name %in% colnames(asv_leaf_16s)),]
metadata_16s_o <- metadata_16s_o[which(metadata_16s_o$sample_name %in% colnames(asv_o_16s)),]
metadata_16s_m <- metadata_16s_m[which(metadata_16s_m$sample_name %in% colnames(asv_m_16s)),]
metadata_16s_root <- metadata_16s_root[which(metadata_16s_root$sample_name %in% colnames(asv_root_16s)),]

metadata_16s_leaf <- distinct(metadata_16s_leaf)
metadata_16s_o <- distinct(metadata_16s_o)
metadata_16s_m <- distinct(metadata_16s_m)
metadata_16s_root <- distinct(metadata_16s_root)

rownames(asv_leaf_16s) <- leaf_rare_16s$asv
rownames(asv_o_16s) <- o_rare_16s$asv
rownames(asv_m_16s) <- m_rare_16s$asv
rownames(asv_root_16s) <- root_rare_16s$asv

asv_leaf.t_16s <- t(asv_leaf_16s)
asv_o.t_16s <- t(asv_o_16s)
asv_m.t_16s <- t(asv_m_16s)
asv_root.t_16s <- t(asv_root_16s)

asv_leaf.t_16s <- asv_leaf.t_16s[order(rownames(asv_leaf.t_16s)),]
asv_m.t_16s <- asv_m.t_16s[order(rownames(asv_m.t_16s)),]
asv_o.t_16s <- asv_o.t_16s[order(rownames(asv_o.t_16s)),]
asv_root.t_16s <- asv_root.t_16s[order(rownames(asv_root.t_16s)),]

percent_leaf_16s <- asv_leaf_16s/mean(colSums(asv_leaf_16s)) * 100
percent_o_16s <- asv_leaf_16s/mean(colSums(asv_o_16s)) * 100
percent_m_16s <- asv_leaf_16s/mean(colSums(asv_m_16s)) * 100
percent_root_16s <- asv_leaf_16s/mean(colSums(asv_root_16s)) * 100

# run indicspecies
inv_leaf_16s <- multipatt(asv_leaf.t_16s, cluster=factor(metadata_16s_leaf$tree_pit_type), func = "r.g", duleg=TRUE)
inv_o_16s <- multipatt(asv_o.t_16s, cluster=factor(metadata_16s_o$tree_pit_type), func = "r.g", duleg=TRUE)
inv_m_16s <- multipatt(asv_m.t_16s, cluster=factor(metadata_16s_m$tree_pit_type), func = "r.g", duleg=TRUE)
inv_root_16s <- multipatt(asv_root.t_16s, cluster=factor(metadata_16s_root$tree_pit_type), func = "r.g", duleg=TRUE)

# save information
options(max.print=100000)
sink("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_leaf_16s.csv")
summary(inv_leaf_16s)
sink()

sink("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_o_16s.csv")
summary(inv_o_16s)
sink()

sink("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_m_16s.csv")
summary(inv_m_16s)
sink()

sink("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_root_16s.csv")
summary(inv_root_16s)
sink()

str_leaf_16s <- inv_leaf_16s$str
write.csv(str_leaf_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/correlation.inv_leaf_16s.csv")

str_o_16s <- inv_o_16s$str
write.csv(str_o_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/correlation.inv_o_16s.csv")

str_m_16s <- inv_m_16s$str
write.csv(str_m_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/correlation.inv_m_16s.csv")

str_root_16s <- inv_root_16s$str
write.csv(str_root_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/correlation.inv_root_16s.csv")

# get indicator asvs of each location group
summary_leaf_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_leaf_16s.csv")
rows_leaf_16s <- grep("Group", summary_leaf_16s$Multilevel.pattern.analysis)

summary_o_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_o_16s.csv")
rows_o_16s <- grep("Group", summary_o_16s$Multilevel.pattern.analysis)

summary_m_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_m_16s.csv")
rows_m_16s <- grep("Group", summary_m_16s$Multilevel.pattern.analysis)

summary_root_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/summary.inv_root_16s.csv")
rows_root_16s <- grep("Group", summary_root_16s$Multilevel.pattern.analysis)

leaf_street_tree_16s <- data.frame(summary_leaf_16s[(rows_leaf_16s[1]+1):(rows_leaf_16s[2]-1),])
leaf_urban_edge_16s <- data.frame(summary_leaf_16s[(rows_leaf_16s[2]+1):(rows_leaf_16s[3]-1),])
leaf_urban_interior_16s <- data.frame(summary_leaf_16s[(rows_leaf_16s[3]+1):(rows_leaf_16s[4]-1),])
leaf_rural_edge_16s <- data.frame(summary_leaf_16s[(rows_leaf_16s[4]+1):(rows_leaf_16s[5]-1),])
leaf_rural_interior_16s <- data.frame(summary_leaf_16s[(rows_leaf_16s[5]+1):(nrow(summary_leaf_16s)-2),])

o_street_tree_16s <- data.frame(summary_o_16s[(rows_o_16s[1]+1):(rows_o_16s[2]-1),])
o_urban_edge_16s <- data.frame(summary_o_16s[(rows_o_16s[2]+1):(rows_o_16s[3]-1),])
o_urban_interior_16s <- data.frame(summary_o_16s[(rows_o_16s[3]+1):(rows_o_16s[4]-1),])
o_rural_edge_16s <- data.frame(summary_o_16s[(rows_o_16s[4]+1):(rows_o_16s[5]-1),])
o_rural_interior_16s <- data.frame(summary_o_16s[(rows_o_16s[5]+1):(nrow(summary_o_16s)-2),])

m_street_tree_16s <- data.frame(summary_m_16s[(rows_m_16s[1]+1):(rows_m_16s[2]-1),])
m_urban_edge_16s <- data.frame(summary_m_16s[(rows_m_16s[2]+1):(rows_m_16s[3]-1),])
m_urban_interior_16s <- data.frame(summary_m_16s[(rows_m_16s[3]+1):(rows_m_16s[4]-1),])
m_rural_edge_16s <- data.frame(summary_m_16s[(rows_m_16s[4]+1):(rows_m_16s[5]-1),])
m_rural_interior_16s <- data.frame(summary_m_16s[(rows_m_16s[5]+1):(nrow(summary_m_16s)-2),])

root_street_tree_16s <- data.frame(summary_root_16s[(rows_root_16s[1]+1):(rows_root_16s[2]-1),])
root_urban_edge_16s <- data.frame(summary_root_16s[(rows_root_16s[2]+1):(rows_root_16s[3]-1),])
root_urban_interior_16s <- data.frame(summary_root_16s[(rows_root_16s[3]+1):(rows_root_16s[4]-1),])
root_rural_edge_16s <- data.frame(summary_root_16s[(rows_root_16s[4]+1):(rows_root_16s[5]-1),])
root_rural_interior_16s <- data.frame(summary_root_16s[(rows_root_16s[5]+1):(nrow(summary_root_16s)-2),])

name_16s <- c("Street.Tree", "Urban.Edge", "Urban.Interior", "Rural.Edge", "Rural.Interior")

sample.list_leaf_16s <- list()
sample.list_leaf_16s[[1]] <- leaf_street_tree_16s
sample.list_leaf_16s[[2]] <- leaf_urban_edge_16s
sample.list_leaf_16s[[3]] <- leaf_urban_interior_16s
sample.list_leaf_16s[[4]] <- leaf_rural_edge_16s
sample.list_leaf_16s[[5]] <- leaf_street_tree_16s

sample.list_o_16s <- list()
sample.list_o_16s[[1]] <- o_street_tree_16s
sample.list_o_16s[[2]] <- o_urban_edge_16s
sample.list_o_16s[[3]] <- o_urban_interior_16s
sample.list_o_16s[[4]] <- o_rural_edge_16s
sample.list_o_16s[[5]] <- o_street_tree_16s

sample.list_m_16s <- list()
sample.list_m_16s[[1]] <- m_street_tree_16s
sample.list_m_16s[[2]] <- m_urban_edge_16s
sample.list_m_16s[[3]] <- m_urban_interior_16s
sample.list_m_16s[[4]] <- m_rural_edge_16s
sample.list_m_16s[[5]] <- m_street_tree_16s

sample.list_root_16s <- list()
sample.list_root_16s[[1]] <- root_street_tree_16s
sample.list_root_16s[[2]] <- root_urban_edge_16s
sample.list_root_16s[[3]] <- root_urban_interior_16s
sample.list_root_16s[[4]] <- root_rural_edge_16s
sample.list_root_16s[[5]] <- root_street_tree_16s

# save
for(i in 1:5){
  sample.list_leaf_16s[[i]] <- str_split_fixed(sample.list_leaf_16s[[i]][2:nrow(sample.list_leaf_16s[[i]]),], " 0.",3)
  # rownames(sample.list_leaf_16s[[i]]) <- sample.list_leaf_16s[[i]][,1]
  # rownames(sample.list_leaf_16s[[i]]) <- gsub(" ", "", rownames(sample.list_leaf_16s[[i]]))
  # sample.list_leaf_16s[[i]][,1] <- rownames(sample.list_leaf_16s[[i]])
  colnames(sample.list_leaf_16s[[i]]) <- c("ASV_ID", "stat", "p.value")
  write.csv(sample.list_leaf_16s[[i]], paste0("/projectnb/talbot-lab-data/Katies_data/Street_Trees/",name_16s[i],"_indic_leaf_16s.csv"), row.names = F)
}

for(i in 1:5){
  sample.list_o_16s[[i]] <- str_split_fixed(sample.list_o_16s[[i]][2:nrow(sample.list_o_16s[[i]]),], " 0.",3)
  # rownames(sample.list_o_16s[[i]]) <- sample.list_o_16s[[i]][,1]
  # rownames(sample.list_o_16s[[i]]) <- gsub(" ", "", rownames(sample.list_o_16s[[i]]))
  # sample.list_o_16s[[i]][,1] <- rownames(sample.list_o_16s[[i]])
  colnames(sample.list_o_16s[[i]]) <- c("ASV_ID", "stat", "p.value")
  write.csv(sample.list_o_16s[[i]], paste0("/projectnb/talbot-lab-data/Katies_data/Street_Trees/",name_16s[i],"_indic_o_16s.csv"), row.names = F)
}

for(i in 1:5){
  sample.list_m_16s[[i]] <- str_split_fixed(sample.list_m_16s[[i]][2:nrow(sample.list_m_16s[[i]]),], " 0.",3)
  # rownames(sample.list_m_16s[[i]]) <- sample.list_m_16s[[i]][,1]
  # rownames(sample.list_m_16s[[i]]) <- gsub(" ", "", rownames(sample.list_m_16s[[i]]))
  # sample.list_m_16s[[i]][,1] <- rownames(sample.list_m_16s[[i]])
  colnames(sample.list_m_16s[[i]]) <- c("ASV_ID", "stat", "p.value")
  write.csv(sample.list_m_16s[[i]], paste0("/projectnb/talbot-lab-data/Katies_data/Street_Trees/",name_16s[i],"_indic_m_16s.csv"), row.names = F)
}

for(i in 1:5){
  sample.list_root_16s[[i]] <- str_split_fixed(sample.list_root_16s[[i]][2:nrow(sample.list_root_16s[[i]]),], " 0.",3)
  # rownames(sample.list_root_16s[[i]]) <- sample.list_root_16s[[i]][,1]
  # rownames(sample.list_root_16s[[i]]) <- gsub(" ", "", rownames(sample.list_root_16s[[i]]))
  # sample.list_root_16s[[i]][,1] <- rownames(sample.list_root_16s[[i]])
  colnames(sample.list_root_16s[[i]]) <- c("ASV_ID", "stat", "p.value")
  write.csv(sample.list_root_16s[[i]], paste0("/projectnb/talbot-lab-data/Katies_data/Street_Trees/",name_16s[i],"_indic_root_16s.csv"), row.names = F)
}

cols_in_database_16s <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "asv", "guild")
tax_leaf_16s <- leaf_rare_16s[, colnames(leaf_rare_16s) %in% cols_in_database_16s]
tax_o_16s <- o_rare_16s[, colnames(o_rare_16s) %in% cols_in_database_16s]
tax_m_16s <- m_rare_16s[, colnames(m_rare_16s) %in% cols_in_database_16s]
tax_root_16s <- root_rare_16s[, colnames(root_rare_16s) %in% cols_in_database_16s]

leaf_street_tree_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Street.Tree_indic_leaf_16s.csv", header = T)
leaf_urban_edge_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Edge_indic_leaf_16s.csv", header = T)
leaf_urban_interior_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Interior_indic_leaf_16s.csv", header = T)
leaf_rural_edge_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Edge_indic_leaf_16s.csv", header = T)
leaf_rural_interior_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Interior_indic_leaf_16s.csv", header = T)

o_street_tree_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Street.Tree_indic_o_16s.csv", header = T)
o_urban_edge_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Edge_indic_o_16s.csv", header = T)
o_urban_interior_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Interior_indic_o_16s.csv", header = T)
o_rural_edge_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Edge_indic_o_16s.csv", header = T)
o_rural_interior_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Interior_indic_o_16s.csv", header = T)

m_street_tree_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Street.Tree_indic_m_16s.csv", header = T)
m_urban_edge_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Edge_indic_m_16s.csv", header = T)
m_urban_interior_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Interior_indic_m_16s.csv", header = T)
m_rural_edge_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Edge_indic_m_16s.csv", header = T)
m_rural_interior_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Interior_indic_m_16s.csv", header = T)

root_street_tree_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Street.Tree_indic_root_16s.csv", header = T)
root_urban_edge_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Edge_indic_root_16s.csv", header = T)
root_urban_interior_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Urban.Interior_indic_root_16s.csv", header = T)
root_rural_edge_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Edge_indic_root_16s.csv", header = T)
root_rural_interior_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Rural.Interior_indic_root_16s.csv", header = T)

leaf_street_tree_16s_merged <- merge(leaf_street_tree_16s, tax_leaf_16s, by.x = "ASV_ID", by.y = "asv")
leaf_urban_edge_16s_merged <- merge(leaf_urban_edge_16s, tax_leaf_16s, by.x = "ASV_ID", by.y = "asv")
leaf_urban_interior_16s_merged <- merge(leaf_urban_interior_16s, tax_leaf_16s, by.x = "ASV_ID", by.y = "asv")
leaf_rural_edge_16s_merged <- merge(leaf_rural_edge_16s, tax_leaf_16s, by.x = "ASV_ID", by.y = "asv")
leaf_rural_interior_16s_merged <- merge(leaf_rural_interior_16s, tax_leaf_16s, by.x = "ASV_ID", by.y = "asv")

o_street_tree_16s_merged <- merge(o_street_tree_16s, tax_o_16s, by.x = "ASV_ID", by.y = "asv")
o_urban_edge_16s_merged <- merge(o_urban_edge_16s, tax_o_16s, by.x = "ASV_ID", by.y = "asv")
o_urban_interior_16s_merged <- merge(o_urban_interior_16s, tax_o_16s, by.x = "ASV_ID", by.y = "asv")
o_rural_edge_16s_merged <- merge(o_rural_edge_16s, tax_o_16s, by.x = "ASV_ID", by.y = "asv")
o_rural_interior_16s_merged <- merge(o_rural_interior_16s, tax_o_16s, by.x = "ASV_ID", by.y = "asv")

m_street_tree_16s_merged <- merge(m_street_tree_16s, tax_m_16s, by.x = "ASV_ID", by.y = "asv")
m_urban_edge_16s_merged <- merge(m_urban_edge_16s, tax_m_16s, by.x = "ASV_ID", by.y = "asv")
m_urban_interior_16s_merged <- merge(m_urban_interior_16s, tax_m_16s, by.x = "ASV_ID", by.y = "asv")
m_rural_edge_16s_merged <- merge(m_rural_edge_16s, tax_m_16s, by.x = "ASV_ID", by.y = "asv")
m_rural_interior_16s_merged <- merge(m_rural_interior_16s, tax_m_16s, by.x = "ASV_ID", by.y = "asv")

root_street_tree_16s_merged <- merge(root_street_tree_16s, tax_root_16s, by.x = "ASV_ID", by.y = "asv")
root_urban_edge_16s_merged <- merge(root_urban_edge_16s, tax_root_16s, by.x = "ASV_ID", by.y = "asv")
root_urban_interior_16s_merged <- merge(root_urban_interior_16s, tax_root_16s, by.x = "ASV_ID", by.y = "asv")
root_rural_edge_16s_merged <- merge(root_rural_edge_16s, tax_root_16s, by.x = "ASV_ID", by.y = "asv")
root_rural_interior_16s_merged <- merge(root_rural_interior_16s, tax_root_16s, by.x = "ASV_ID", by.y = "asv")

rownames(leaf_street_tree_16s_merged) <- leaf_street_tree_16s_merged[,1]
rownames(leaf_urban_edge_16s_merged) <- leaf_urban_edge_16s_merged[,1]
rownames(leaf_urban_interior_16s_merged) <- leaf_urban_interior_16s_merged[,1]
rownames(leaf_rural_edge_16s_merged) <- leaf_rural_edge_16s_merged[,1]
rownames(leaf_rural_interior_16s_merged) <- leaf_rural_interior_16s_merged[,1]

rownames(o_street_tree_16s_merged) <- o_street_tree_16s_merged[,1]
rownames(o_urban_edge_16s_merged) <- o_urban_edge_16s_merged[,1]
rownames(o_urban_interior_16s_merged) <- o_urban_interior_16s_merged[,1]
rownames(o_rural_edge_16s_merged) <- o_rural_edge_16s_merged[,1]
rownames(o_rural_interior_16s_merged) <- o_rural_interior_16s_merged[,1]

rownames(m_street_tree_16s_merged) <- m_street_tree_16s_merged[,1]
rownames(m_urban_edge_16s_merged) <- m_urban_edge_16s_merged[,1]
rownames(m_urban_interior_16s_merged) <- m_urban_interior_16s_merged[,1]
rownames(m_rural_edge_16s_merged) <- m_rural_edge_16s_merged[,1]
rownames(m_rural_interior_16s_merged) <- m_rural_interior_16s_merged[,1]

rownames(root_street_tree_16s_merged) <- root_street_tree_16s_merged[,1]
rownames(root_urban_edge_16s_merged) <- root_urban_edge_16s_merged[,1]
rownames(root_urban_interior_16s_merged) <- root_urban_interior_16s_merged[,1]
rownames(root_rural_edge_16s_merged) <- root_rural_edge_16s_merged[,1]
rownames(root_rural_interior_16s_merged) <- root_rural_interior_16s_merged[,1]

leaf_street_tree_16s_merged$Group <- "Street.tree"
leaf_urban_edge_16s_merged$Group <- "Urban.edge"
leaf_urban_interior_16s_merged$Group <- "Urban.interior"
leaf_rural_edge_16s_merged$Group <- "Rural.edge"
leaf_rural_interior_16s_merged$Group <- "Rural.interior"

o_street_tree_16s_merged$Group <- "Street.tree"
o_urban_edge_16s_merged$Group <- "Urban.edge"
o_urban_interior_16s_merged$Group <- "Urban.interior"
o_rural_edge_16s_merged$Group <- "Rural.edge"
o_rural_interior_16s_merged$Group <- "Rural.interior"

m_street_tree_16s_merged$Group <- "Street.tree"
m_urban_edge_16s_merged$Group <- "Urban.edge"
m_urban_interior_16s_merged$Group <- "Urban.interior"
m_rural_edge_16s_merged$Group <- "Rural.edge"
m_rural_interior_16s_merged$Group <- "Rural.interior"

root_street_tree_16s_merged$Group <- "Street.tree"
root_urban_edge_16s_merged$Group <- "Urban.edge"
root_urban_interior_16s_merged$Group <- "Urban.interior"
root_rural_edge_16s_merged$Group <- "Rural.edge"
root_rural_interior_16s_merged$Group <- "Rural.interior"

leaf_street_tree_16s_merged$Sample.Type <- "Leaf"
leaf_urban_edge_16s_merged$Sample.Type <- "Leaf"
leaf_urban_interior_16s_merged$Sample.Type <- "Leaf"
leaf_rural_edge_16s_merged$Sample.Type <- "Leaf"
leaf_rural_interior_16s_merged$Sample.Type <- "Leaf"

o_street_tree_16s_merged$Sample.Type <- "O"
o_urban_edge_16s_merged$Sample.Type <- "O"
o_urban_interior_16s_merged$Sample.Type <- "O"
o_rural_edge_16s_merged$Sample.Type <- "O"
o_rural_interior_16s_merged$Sample.Type <- "O"

m_street_tree_16s_merged$Sample.Type <- "M"
m_urban_edge_16s_merged$Sample.Type <- "M"
m_urban_interior_16s_merged$Sample.Type <- "M"
m_rural_edge_16s_merged$Sample.Type <- "M"
m_rural_interior_16s_merged$Sample.Type <- "M"

root_street_tree_16s_merged$Sample.Type <- "Root"
root_urban_edge_16s_merged$Sample.Type <- "Root"
root_urban_interior_16s_merged$Sample.Type <- "Root"
root_rural_edge_16s_merged$Sample.Type <- "Root"
root_rural_interior_16s_merged$Sample.Type <- "Root"

leaf_street_tree_16s <- leaf_street_tree_16s_merged[,-c(1:3)]
leaf_urban_edge_16s <- leaf_urban_edge_16s_merged[,-c(1:3)]
leaf_urban_interior_16s <- leaf_urban_interior_16s_merged[,-c(1:3)]
leaf_rural_edge_16s <- leaf_rural_edge_16s_merged[,-c(1:3)]
leaf_rural_interior_16s <- leaf_rural_interior_16s_merged[,-c(1:3)]

o_street_tree_16s <- o_street_tree_16s_merged[,-c(1:3)]
o_urban_edge_16s <- o_urban_edge_16s_merged[,-c(1:3)]
o_urban_interior_16s <- o_urban_interior_16s_merged[,-c(1:3)]
o_rural_edge_16s <- o_rural_edge_16s_merged[,-c(1:3)]
o_rural_interior_16s <- o_rural_interior_16s_merged[,-c(1:3)]

m_street_tree_16s <- m_street_tree_16s_merged[,-c(1:3)]
m_urban_edge_16s <- m_urban_edge_16s_merged[,-c(1:3)]
m_urban_interior_16s <- m_urban_interior_16s_merged[,-c(1:3)]
m_rural_edge_16s <- m_rural_edge_16s_merged[,-c(1:3)]
m_rural_interior_16s <- m_rural_interior_16s_merged[,-c(1:3)]

root_street_tree_16s <- root_street_tree_16s_merged[,-c(1:3)]
root_urban_edge_16s <- root_urban_edge_16s_merged[,-c(1:3)]
root_urban_interior_16s <- root_urban_interior_16s_merged[,-c(1:3)]
root_rural_edge_16s <- root_rural_edge_16s_merged[,-c(1:3)]
root_rural_interior_16s <- root_rural_interior_16s_merged[,-c(1:3)]

# save sequences
rep_seqs_16s <- readDNAStringSet("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/dada2_output/16s_rep_unique.fasta")
seq_name_16s <- names(rep_seqs_16s)
sequence_16s <- paste(rep_seqs_16s)
rep_seqs_16s <- data.frame(seq_name_16s, sequence_16s)
all_leaf_16s <- rbind(leaf_street_tree_16s, leaf_urban_edge_16s, leaf_urban_interior_16s, leaf_rural_edge_16s, leaf_rural_interior_16s)
all_o_16s <- rbind(o_street_tree_16s, o_urban_edge_16s, o_urban_interior_16s, o_rural_edge_16s, o_rural_interior_16s)
all_m_16s <- rbind(m_street_tree_16s, m_urban_edge_16s, m_urban_interior_16s, m_rural_edge_16s, m_rural_interior_16s)
all_root_16s <- rbind(root_street_tree_16s, root_urban_edge_16s, root_urban_interior_16s, root_rural_edge_16s, root_rural_interior_16s)

write.csv(all_leaf_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_indic_phylogenetic_info.csv")
write.csv(all_o_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_indic_phylogenetic_info.csv")
write.csv(all_m_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_indic_phylogenetic_info.csv")
write.csv(all_root_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_indic_phylogenetic_info.csv")

all_leaf_fasta_16s <- merge(all_leaf_16s, rep_seqs_16s, by.x = 0, by.y = "seq_name_16s")
all_o_fasta_16s <- merge(all_o_16s, rep_seqs_16s, by.x = 0, by.y = "seq_name_16s")
all_m_fasta_16s <- merge(all_m_16s, rep_seqs_16s, by.x = 0, by.y = "seq_name_16s")
all_root_fasta_16s <- merge(all_root_16s, rep_seqs_16s, by.x = 0, by.y = "seq_name_16s")

write.fasta(sequences = as.list(all_leaf_fasta_16s$sequence_16s), names = all_leaf_fasta_16s$Row.names, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_indic_seqs.fasta")
write.fasta(sequences = as.list(all_o_fasta_16s$sequence_16s), names = all_o_fasta_16s$Row.names, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_indic_seqs.fasta")
write.fasta(sequences = as.list(all_m_fasta_16s$sequence_16s), names = all_m_fasta_16s$Row.names, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_indic_seqs.fasta")
write.fasta(sequences = as.list(all_root_fasta_16s$sequence_16s), names = all_root_fasta_16s$Row.names, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_indic_seqs.fasta")

# pathogen_status 
tax_w_paths <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Raw_Data/16S/tax_16s_allruns_wpathogen.csv")
tax_w_paths <- tax_w_paths[,c(2,10)]

leaf_phylo <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/bacteria_leaf_indic_phylogenetic_info.csv")
root_phylo <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/bacteria_root_indic_phylogenetic_info.csv")
o_phylo <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/bacteria_o_indic_phylogenetic_info.csv")
m_phylo <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/bacteria_m_indic_phylogenetic_info.csv")

colnames(leaf_phylo)[1] <- "asv_id"
colnames(root_phylo)[1] <- "asv_id"
colnames(m_phylo)[1] <- "asv_id"
colnames(o_phylo)[1] <- "asv_id"

leaf_phylo_pathogen <- merge(leaf_phylo, tax_w_paths, by = "asv_id", all.x = TRUE)
root_phylo_pathogen <- merge(root_phylo, tax_w_paths, by = "asv_id", all.x = TRUE)
m_phylo_pathogen <- merge(m_phylo, tax_w_paths, by = "asv_id", all.x = TRUE)
o_phylo_pathogen <- merge(o_phylo, tax_w_paths, by = "asv_id", all.x = TRUE)

write.csv(leaf_phylo_pathogen, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/bacteria_leaf_indic_phylogenetic_info.csv")
write.csv(root_phylo_pathogen, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/bacteria_root_indic_phylogenetic_info.csv")
write.csv(m_phylo_pathogen, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/bacteria_m_indic_phylogenetic_info.csv")
write.csv(o_phylo_pathogen, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/bacteria_o_indic_phylogenetic_info.csv")

# get 28S sequence
sample_leaf_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_indic_phylogenetic_info.csv", header = T)
sample_o_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_indic_phylogenetic_info.csv", header = T)
sample_m_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_indic_phylogenetic_info.csv", header = T)
sample_root_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_indic_phylogenetic_info.csv", header = T)

sample_leaf_16s <- data.frame(lapply(sample_leaf_16s, function(x){gsub(pattern="o__", replacement = "", x)}), stringsAsFactors = FALSE)
sample_o_16s <- data.frame(lapply(sample_o_16s, function(x){gsub(pattern="o__", replacement = "", x)}), stringsAsFactors = FALSE)
sample_m_16s <- data.frame(lapply(sample_m_16s, function(x){gsub(pattern="o__", replacement = "", x)}), stringsAsFactors = FALSE)
sample_root_16s <- data.frame(lapply(sample_root_16s, function(x){gsub(pattern="o__", replacement = "", x)}), stringsAsFactors = FALSE)

colnames(sample_leaf_16s)[1] <- "ASV_ID"
colnames(sample_o_16s)[1] <- "ASV_ID"
colnames(sample_m_16s)[1] <- "ASV_ID"
colnames(sample_root_16s)[1] <- "ASV_ID"

database_16s <- read.fasta("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/SILVA_138.1_LSURef_NR99_tax_silva.fasta", as.string=TRUE)
at_16s <- data.frame(attributes(database_16s))

# format database
tax.list_16s <- NULL
id.list_16s <- NULL
seq.list_16s <- NULL
for(i in 1:nrow(at_16s)){
  description <- attributes(database_16s[[at_16s[i,]]])[["Annot"]]
  annot <- strsplit(description, " +")[[1]][2]
  if(grepl("Fungi", annot, fixed = TRUE)){
    id <- at_16s[i,]
    id.list_16s <- c(id.list_16s, id)
    seq <- database_16s[[at_16s[i,]]]
    seq.list_16s <- rbind(seq.list_16s, seq)
    tax.name <- strsplit(annot, ";")[[1]]
    tax.name <- c(tax.name, rep(NA, (14-length(tax.name))))
    tax.list_16s <- rbind(tax.list_16s, tax.name)}
}
tax.data_16s <- data.frame(tax.list_16s)
seq.data_16s <- data.frame(seq.list_16s)

database_16s <- tax.data_16s[,6:14]
database_16s$ID <- id.list_16s
database_16s$seq <- seq.data_16s$seq.list_16s
database_16s <- data.frame(database_16s)

colnames(database_16s) <- c("Kingdom","Subkingdom","Phylum","Subphylum","Class","Order","Family","Genus","Species","ID","Sequence")

# assign 28S sequences to ASVs referring to Genus
database_16s.genus <- na.omit(cbind(database_16s$Genus, database_16s$Sequence))
colnames(database_16s.genus) <- c("Genus", "Sequence")
database_16s.genus <- data.frame(database_16s.genus)
db_16s.genus <- database_16s.genus[!duplicated(database_16s.genus$Genus),]

merge_leaf_16s.genus <- merge(sample_leaf_16s, db_16s.genus, by.x = "GENUS", by.y = "Genus", all.x = TRUE)
merge_leaf_16s.genus <- merge_leaf_16s.genus[order(merge_leaf_16s.genus$ASV_ID, decreasing = FALSE),]
merge_o_16s.genus <- merge(sample_o_16s, db_16s.genus, by.x = "GENUS", by.y = "Genus", all.x = TRUE)
merge_o_16s.genus <- merge_o_16s.genus[order(merge_o_16s.genus$ASV_ID, decreasing = FALSE),]
merge_m_16s.genus <- merge(sample_m_16s, db_16s.genus, by.x = "GENUS", by.y = "Genus", all.x = TRUE)
merge_m_16s.genus <- merge_m_16s.genus[order(merge_m_16s.genus$ASV_ID, decreasing = FALSE),]
merge_root_16s.genus <- merge(sample_root_16s, db_16s.genus, by.x = "GENUS", by.y = "Genus", all.x = TRUE)
merge_root_16s.genus <- merge_root_16s.genus[order(merge_root_16s.genus$ASV_ID, decreasing = FALSE),]

# assign 28S sequences to ASVs referring to Family
database_16s.family <- na.omit(cbind(database_16s$Family, database_16s$Sequence))
colnames(database_16s.family) <- c("Family", "Sequence")
database_16s.family <- data.frame(database_16s.family)
db_16s.family <- database_16s.family[!duplicated(database_16s.family$Family),]

merge_leaf_16s.family <- merge(sample_leaf_16s, db_16s.family, by.x = "Family", by.y = "Family", all.x = TRUE)
merge_leaf_16s.family <- merge_leaf_16s.family[order(merge_leaf_16s.family$ASV_ID, decreasing = FALSE),]
merge_o_16s.family <- merge(sample_o_16s, db_16s.family, by.x = "Family", by.y = "Family", all.x = TRUE)
merge_o_16s.family <- merge_o_16s.family[order(merge_o_16s.family$ASV_ID, decreasing = FALSE),]
merge_m_16s.family <- merge(sample_m_16s, db_16s.family, by.x = "Family", by.y = "Family", all.x = TRUE)
merge_m_16s.family <- merge_m_16s.family[order(merge_m_16s.family$ASV_ID, decreasing = FALSE),]
merge_root_16s.family <- merge(sample_root_16s, db_16s.family, by.x = "Family", by.y = "Family", all.x = TRUE)
merge_root_16s.family <- merge_root_16s.family[order(merge_root_16s.family$ASV_ID, decreasing = FALSE),]

# assign 28S sequences to ASVs referring to Order
database_16s.order <- na.omit(cbind(database_16s$Order, database_16s$Sequence))
colnames(database_16s.order) <- c("Order", "Sequence")
database_16s.order <- data.frame(database_16s.order)
db_16s.order <- database_16s.order[!duplicated(database_16s.order$Order),]

merge_leaf_16s.order <- merge(sample_leaf_16s, db_16s.order, by.x = "Order", by.y = "Order", all.x = TRUE)
merge_leaf_16s.order <- merge_leaf_16s.order[order(merge_leaf_16s.order$ASV_ID, decreasing = FALSE),]
merge_o_16s.order <- merge(sample_o_16s, db_16s.order, by.x = "Order", by.y = "Order", all.x = TRUE)
merge_o_16s.order <- merge_o_16s.order[order(merge_o_16s.order$ASV_ID, decreasing = FALSE),]
merge_m_16s.order <- merge(sample_m_16s, db_16s.order, by.x = "Order", by.y = "Order", all.x = TRUE)
merge_m_16s.order <- merge_m_16s.order[order(merge_m_16s.order$ASV_ID, decreasing = FALSE),]
merge_root_16s.order <- merge(sample_root_16s, db_16s.order, by.x = "Order", by.y = "Order", all.x = TRUE)
merge_root_16s.order <- merge_root_16s.order[order(merge_root_16s.order$ASV_ID, decreasing = FALSE),]

# combine taxonomy information
data_leaf_16s <- sample_leaf_16s[order(sample_leaf_16s$ASV_ID, decreasing = FALSE),]
data_leaf_16s$seqs <- merge_leaf_16s.genus$Sequence
data_o_16s <- sample_o_16s[order(sample_o_16s$ASV_ID, decreasing = FALSE),]
data_o_16s$seqs <- merge_o_16s.genus$Sequence
data_m_16s <- sample_m_16s[order(sample_m_16s$ASV_ID, decreasing = FALSE),]
data_m_16s$seqs <- merge_m_16s.genus$Sequence
data_root_16s <- sample_root_16s[order(sample_root_16s$ASV_ID, decreasing = FALSE),]
data_root_16s$seqs <- merge_root_16s.genus$Sequence

for(i in 1:nrow(data_leaf_16s)){
  if(is.na(data_leaf_16s$seqs[i])){data_leaf_16s$seqs[i] <- merge_leaf_16s.family$Sequence[i]}
  if(is.na(data_leaf_16s$seqs[i])){data_leaf_16s$seqs[i] <- merge_leaf_16s.order$Sequence[i]}}
for(i in 1:nrow(data_o_16s)){
  if(is.na(data_o_16s$seqs[i])){data_o_16s$seqs[i] <- merge_o_16s.family$Sequence[i]}
  if(is.na(data_o_16s$seqs[i])){data_o_16s$seqs[i] <- merge_o_16s.order$Sequence[i]}}
for(i in 1:nrow(data_m_16s)){
  if(is.na(data_m_16s$seqs[i])){data_m_16s$seqs[i] <- merge_m_16s.family$Sequence[i]}
  if(is.na(data_m_16s$seqs[i])){data_m_16s$seqs[i] <- merge_m_16s.order$Sequence[i]}}
for(i in 1:nrow(data_root_16s)){
  if(is.na(data_root_16s$seqs[i])){data_root_16s$seqs[i] <- merge_root_16s.family$Sequence[i]}
  if(is.na(data_root_16s$seqs[i])){data_root_16s$seqs[i] <- merge_root_16s.order$Sequence[i]}}

data_leaf_16s <- data_leaf_16s[!is.na(data_leaf_16s$seqs),]
data_o_16s <- data_o_16s[!is.na(data_o_16s$seqs),]
data_m_16s <- data_m_16s[!is.na(data_m_16s$seqs),]
data_root_16s <- data_root_16s[!is.na(data_root_16s$seqs),]

rownames(data_leaf_16s) <- data_leaf_16s$ASV_ID
rownames(data_o_16s) <- data_o_16s$ASV_ID
rownames(data_m_16s) <- data_m_16s$ASV_ID
rownames(data_root_16s) <- data_root_16s$ASV_ID

write.csv(data_leaf_16s[,-1], "/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_phylogenetic_info_28S.csv")
write.csv(data_o_16s[,-1], "/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_phylogenetic_info_28S.csv")
write.csv(data_m_16s[,-1], "/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_phylogenetic_info_28S.csv")
write.csv(data_root_16s[,-1], "/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_phylogenetic_info_28S.csv")

write.fasta(sequences = as.list(data_leaf_16s$seqs), names = data_leaf_16s$ASV_ID, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_seqs_28s.fasta")
write.fasta(sequences = as.list(data_o_16s$seqs), names = data_o_16s$ASV_ID, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_seqs_28s.fasta")
write.fasta(sequences = as.list(data_m_16s$seqs), names = data_m_16s$ASV_ID, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_seqs_28s.fasta")
write.fasta(sequences = as.list(data_root_16s$seqs), names = data_root_16s$ASV_ID, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_seqs_28s.fasta")

# to visualize in iTOL
# make a binary file to show their presence or absence
data_leaf_16s$Street.tree <- 0
data_leaf_16s$Urban.edge <- 0
data_leaf_16s$Urban.interior <- 0
data_leaf_16s$Rural.edge <- 0
data_leaf_16s$Rural.interior <- 0

data_o_16s$Street.tree <- 0
data_o_16s$Urban.edge <- 0
data_o_16s$Urban.interior <- 0
data_o_16s$Rural.edge <- 0
data_o_16s$Rural.interior <- 0

data_m_16s$Street.tree <- 0
data_m_16s$Urban.edge <- 0
data_m_16s$Urban.interior <- 0
data_m_16s$Rural.edge <- 0
data_m_16s$Rural.interior <- 0

data_root_16s$Street.tree <- 0
data_root_16s$Urban.edge <- 0
data_root_16s$Urban.interior <- 0
data_root_16s$Rural.edge <- 0
data_root_16s$Rural.interior <- 0

for(i in 1:nrow(data_leaf_16s)){
  if(data_leaf_16s$Group[i] == "Street.tree"){data_leaf_16s$Street.tree[i] <- 1}else{data_leaf_16s$Street.tree[i] <- 0}
  if(data_leaf_16s$Group[i] == "Urban.edge"){data_leaf_16s$Urban.edge[i] <- 1}else{data_leaf_16s$Urban.edge[i] <- 0}
  if(data_leaf_16s$Group[i] == "Urban.interior"){data_leaf_16s$Urban.interior[i] <- 1}else{data_leaf_16s$Urban.interior[i] <- 0}
  if(data_leaf_16s$Group[i] == "Rural.edge"){data_leaf_16s$Rural.edge[i] <- 1}else{data_leaf_16s$Rural.edge[i] <- 0}
  if(data_leaf_16s$Group[i] == "Rural.interior"){data_leaf_16s$Rural.interior[i] <- 1}else{data_leaf_16s$Rural.interior[i] <- 0}
}

for(i in 1:nrow(data_o_16s)){
  if(data_o_16s$Group[i] == "Street.tree"){data_o_16s$Street.tree[i] <- 1}else{data_o_16s$Street.tree[i] <- 0}
  if(data_o_16s$Group[i] == "Urban.edge"){data_o_16s$Urban.edge[i] <- 1}else{data_o_16s$Urban.edge[i] <- 0}
  if(data_o_16s$Group[i] == "Urban.interior"){data_o_16s$Urban.interior[i] <- 1}else{data_o_16s$Urban.interior[i] <- 0}
  if(data_o_16s$Group[i] == "Rural.edge"){data_o_16s$Rural.edge[i] <- 1}else{data_o_16s$Rural.edge[i] <- 0}
  if(data_o_16s$Group[i] == "Rural.interior"){data_o_16s$Rural.interior[i] <- 1}else{data_o_16s$Rural.interior[i] <- 0}
}

for(i in 1:nrow(data_m_16s)){
  if(data_m_16s$Group[i] == "Street.tree"){data_m_16s$Street.tree[i] <- 1}else{data_m_16s$Street.tree[i] <- 0}
  if(data_m_16s$Group[i] == "Urban.edge"){data_m_16s$Urban.edge[i] <- 1}else{data_m_16s$Urban.edge[i] <- 0}
  if(data_m_16s$Group[i] == "Urban.interior"){data_m_16s$Urban.interior[i] <- 1}else{data_m_16s$Urban.interior[i] <- 0}
  if(data_m_16s$Group[i] == "Rural.edge"){data_m_16s$Rural.edge[i] <- 1}else{data_m_16s$Rural.edge[i] <- 0}
  if(data_m_16s$Group[i] == "Rural.interior"){data_m_16s$Rural.interior[i] <- 1}else{data_m_16s$Rural.interior[i] <- 0}
}

for(i in 1:nrow(data_root_16s)){
  if(data_root_16s$Group[i] == "Street.tree"){data_root_16s$Street.tree[i] <- 1}else{data_root_16s$Street.tree[i] <- 0}
  if(data_root_16s$Group[i] == "Urban.edge"){data_root_16s$Urban.edge[i] <- 1}else{data_root_16s$Urban.edge[i] <- 0}
  if(data_root_16s$Group[i] == "Urban.interior"){data_root_16s$Urban.interior[i] <- 1}else{data_root_16s$Urban.interior[i] <- 0}
  if(data_root_16s$Group[i] == "Rural.edge"){data_root_16s$Rural.edge[i] <- 1}else{data_root_16s$Rural.edge[i] <- 0}
  if(data_root_16s$Group[i] == "Rural.interior"){data_root_16s$Rural.interior[i] <- 1}else{data_root_16s$Rural.interior[i] <- 0}
}

heatmap_leaf_16s <- cbind(data_leaf_16s$ASV_ID, data_leaf_16s$Street.tree, data_leaf_16s$Urban.edge, data_leaf_16s$Urban.interior, data_leaf_16s$Rural.edge, data_leaf_16s$Rural.interior)
heatmap_o_16s <- cbind(data_o_16s$ASV_ID, data_o_16s$Street.tree, data_o_16s$Urban.edge, data_o_16s$Urban.interior, data_o_16s$Rural.edge, data_o_16s$Rural.interior)
heatmap_m_16s <- cbind(data_m_16s$ASV_ID, data_m_16s$Street.tree, data_m_16s$Urban.edge, data_m_16s$Urban.interior, data_m_16s$Rural.edge, data_m_16s$Rural.interior)
heatmap_root_16s <- cbind(data_root_16s$ASV_ID, data_root_16s$Street.tree, data_root_16s$Urban.edge, data_root_16s$Urban.interior, data_root_16s$Rural.edge, data_root_16s$Rural.interior)

write.csv(heatmap_leaf_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/heatmap_leaf_16s_iTOL_28S.csv")
write.csv(heatmap_o_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/heatmap_o_16s_iTOL_28S.csv")
write.csv(heatmap_m_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/heatmap_m_16s_iTOL_28S.csv")
write.csv(heatmap_root_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/heatmap_root_16s_iTOL_28S.csv")

# tree_16s <- read.tree("Phylogenetic_tree_28S_Newick")
# tree_16s
# 
# data_leaf_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_phylogenetic_info_28S.csv")
# data_o_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_phylogenetic_info_28S.csv")
# data_m_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_phylogenetic_info_28S.csv")
# data_root_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_phylogenetic_info_28S.csv")
# 
# data_leaf_16s <- data.frame(lapply(data_leaf_16s, function(x){gsub(pattern = "p__", replacement = "", x)}), stringsAsFactors = FALSE)
# data_o_16s <- data.frame(lapply(data_o_16s, function(x){gsub(pattern = "p__", replacement = "", x)}), stringsAsFactors = FALSE)
# data_m_16s <- data.frame(lapply(data_m_16s, function(x){gsub(pattern = "p__", replacement = "", x)}), stringsAsFactors = FALSE)
# data_root_16s <- data.frame(lapply(data_root_16s, function(x){gsub(pattern = "p__", replacement = "", x)}), stringsAsFactors = FALSE)
# 
# # to change color by phylum
# groupInfo_leaf_16s <- split(data_leaf_16s$X, data_leaf_16s$Phylum)
# groupInfo_o_16s <- split(data_o_16s$X, data_o_16s$Phylum)
# groupInfo_m_16s <- split(data_m_16s$X, data_m_16s$Phylum)
# groupInfo_root_16s <- split(data_root_16s$X, data_root_16s$Phylum)
# 
# tree_leaf_16s <- groupOTU(tree_16s, groupInfo_leaf_16s, group_name = "Phylum")
# tree_o_16s <- groupOTU(tree_16s, groupInfo_o_16s, group_name = "Phylum")
# tree_m_16s <- groupOTU(tree_16s, groupInfo_m_16s, group_name = "Phylum")
# tree_root_16s <- groupOTU(tree_16s, groupInfo_root_16s, group_name = "Phylum")
# 
# # save
# ggtree(tree_leaf_16s, aes(color = Phylum), layout = 'circular')
# ggtree(tree_o_16s, aes(color = Phylum), layout = 'circular')
# ggtree(tree_m_16s, aes(color = Phylum), layout = 'circular')
# ggtree(tree_root_16s, aes(color = Phylum), layout = 'circular')