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
all_clr_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/CLR/16S/all_16s_asv_clr_concat.csv", row.names = 1)
# all_rare <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/all_16s_asv_rarefied.csv", row.names = 1)
leaf_rare_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Rarefaction/16S/leaf_16s_asv_rare_wtax_guild_200.csv", row.names = 1)
m_rare_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Rarefaction/16S/m_16s_asv_rare_wtax_guild.csv", row.names = 1)
o_rare_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Rarefaction/16S/o_16s_asv_rare_wtax_guild.csv", row.names = 1)
root_rare_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Rarefaction/16S/root_16s_asv_rare_wtax_guild.csv", row.names = 1)

leaf_rare_16s_tax <- leaf_rare_16s[,2:8]
m_rare_16s_tax <- m_rare_16s[,2:8]
o_rare_16s_tax <- o_rare_16s[,2:8]
root_rare_16s_tax <- root_rare_16s[,2:8]

bg_rare_16s_tax <- rbind(m_rare_16s_tax, o_rare_16s_tax, root_rare_16s_tax)
bg_rare_16s_tax <- distinct(bg_rare_16s_tax)
leaf_rare_16s_tax <- distinct(leaf_rare_16s_tax)

leaf_rare_16s <- leaf_rare_16s[,8:(ncol(leaf_rare_16s)-1)]
root_rare_16s <- root_rare_16s[,8:(ncol(root_rare_16s)-1)]
m_rare_16s <- m_rare_16s[,8:(ncol(m_rare_16s)-1)]
o_rare_16s <- o_rare_16s[,8:(ncol(o_rare_16s)-1)]

leaf_rare_16s <- aggregate(leaf_rare_16s, . ~ Species, FUN = "sum")
root_rare_16s <- aggregate(root_rare_16s, . ~ Species, FUN = "sum")
m_rare_16s <- aggregate(m_rare_16s, . ~ Species, FUN = "sum")
o_rare_16s <- aggregate(o_rare_16s, . ~ Species, FUN = "sum")

bg_rare_16s <- merge(root_rare_16s, m_rare_16s, by = "Species", all = T)
bg_rare_16s <- merge(bg_rare_16s, o_rare_16s, by = "Species", all = T)

metadata_16s$sample_type <- gsub("Soil", "Belowground", metadata_16s$sample_type)
metadata_16s$sample_type <- gsub("Root", "Belowground", metadata_16s$sample_type)


metadata_16s$sample_type <- factor(metadata_16s$sample_type, levels = c("Leaf", "Belowground"))


metadata_16s$tree_pit_type <- gsub("Pit", "Street Tree", metadata_16s$tree_pit_type)
metadata_16s$tree_pit_type <- gsub("Grass", "Street Tree", metadata_16s$tree_pit_type)

metadata_16s$tree_pit_type <- factor(metadata_16s$tree_pit_type, 
                                     levels = c("Street Tree", "Urban Edge", 
                                                "Urban Interior", "Rural Edge", 
                                                "Rural Interior"))

# format data
metadata_16s_leaf <- metadata_16s[which(metadata_16s$sample_name %in% colnames(leaf_rare_16s)),]
metadata_16s_bg <- metadata_16s[which(metadata_16s$sample_name %in% colnames(bg_rare_16s)),]
metadata_16s_bg <- distinct(metadata_16s_bg)

asv_leaf_16s <- leaf_rare_16s[,2:(ncol(leaf_rare_16s))]
asv_bg_16s <- bg_rare_16s[,2:(ncol(bg_rare_16s))]

rownames(asv_leaf_16s) <- leaf_rare_16s$Species
rownames(asv_bg_16s) <- bg_rare_16s$Species

asv_leaf.t_16s <- t(asv_leaf_16s)
asv_bg.t_16s <- t(asv_bg_16s)

asv_leaf.t_16s <- asv_leaf.t_16s[order(rownames(asv_leaf.t_16s)),]
asv_bg.t_16s <- asv_bg.t_16s[order(rownames(asv_bg.t_16s)),]

percent_leaf_16s <- asv_leaf_16s/mean(colSums(asv_leaf_16s)) * 100
percent_bg_16s <- asv_bg_16s/mean(colSums(asv_bg_16s)) * 100

asv_leaf.t_16s[is.na(asv_leaf.t_16s)] <- 0
asv_bg.t_16s[is.na(asv_bg.t_16s)] <- 0
# run indicspecies
inv_leaf_16s <- multipatt(asv_leaf.t_16s, cluster=factor(metadata_16s_leaf$tree_pit_type), func = "r.g", duleg=TRUE)
inv_bg_16s <- multipatt(asv_bg.t_16s, cluster=factor(metadata_16s_bg$tree_pit_type), func = "r.g", duleg=TRUE)

# save information
options(max.print=100000)
sink("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/summary.inv_leaf_16s_sp.csv")
summary(inv_leaf_16s)
sink()

sink("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/summary.inv_bg_16s.csv")
summary(inv_bg_16s)
sink()

str_leaf_16s <- inv_leaf_16s$str
write.csv(str_leaf_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/correlation.inv_leaf_16s_sp.csv")

str_bg_16s <- inv_bg_16s$str
write.csv(str_bg_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/correlation.inv_bg_16s.csv")

# get indicator asvs of each location group
summary_leaf_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/summary.inv_leaf_16s_sp.csv")
rows_leaf_16s <- grep(" Group", summary_leaf_16s$Multilevel.pattern.analysis)

summary_bg_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/summary.inv_bg_16s.csv")
rows_bg_16s <- grep(" Group", summary_bg_16s$Multilevel.pattern.analysis)

leaf_street_tree_16s <- data.frame(summary_leaf_16s[(rows_leaf_16s[1]+1):(rows_leaf_16s[2]-1),])
leaf_urban_interior_16s <- data.frame(summary_leaf_16s[(rows_leaf_16s[2]+1):(rows_leaf_16s[3]-1),])
leaf_rural_edge_16s <- data.frame(summary_leaf_16s[(rows_leaf_16s[3]+1):(rows_leaf_16s[4]-1),])
leaf_rural_interior_16s <- data.frame(summary_leaf_16s[(rows_leaf_16s[4]+1):(nrow(summary_leaf_16s)-2),])

bg_street_tree_16s <- data.frame(summary_bg_16s[(rows_bg_16s[1]+1):(rows_bg_16s[2]-1),])
bg_urban_edge_16s <- data.frame(summary_bg_16s[(rows_bg_16s[2]+1):(rows_bg_16s[3]-1),])
bg_urban_interior_16s <- data.frame(summary_bg_16s[(rows_bg_16s[3]+1):(rows_bg_16s[4]-1),])
bg_rural_edge_16s <- data.frame(summary_bg_16s[(rows_bg_16s[4]+1):(rows_bg_16s[5]-1),])
bg_rural_interior_16s <- data.frame(summary_bg_16s[(rows_bg_16s[5]+1):(nrow(summary_bg_16s)-2),])

name_16s_bg <- c("Street.Tree", "Urban.Edge", "Urban.Interior", "Rural.Edge", "Rural.Interior")
name_16s_leaf <- c("Street.Tree", "Urban.Interior", "Rural.Edge", "Rural.Interior")

sample.list_leaf_16s <- list()
sample.list_leaf_16s[[1]] <- leaf_street_tree_16s
sample.list_leaf_16s[[2]] <- leaf_urban_interior_16s
sample.list_leaf_16s[[3]] <- leaf_rural_edge_16s
sample.list_leaf_16s[[4]] <- leaf_rural_interior_16s

sample.list_bg_16s <- list()
sample.list_bg_16s[[1]] <- bg_street_tree_16s
sample.list_bg_16s[[2]] <- bg_urban_edge_16s
sample.list_bg_16s[[3]] <- bg_urban_interior_16s
sample.list_bg_16s[[4]] <- bg_rural_edge_16s
sample.list_bg_16s[[5]] <- bg_rural_interior_16s

# save
for(i in 1:4){
  df <- sample.list_leaf_16s[[i]]
  df <- as.data.frame(str_split_fixed(df[2:nrow(df),], " 0.",3))
  # rownames(sample.list_leaf_16s[[i]]) <- sample.list_leaf_16s[[i]][,1]
  # rownames(sample.list_leaf_16s[[i]]) <- gsub(" ", "", rownames(sample.list_leaf_16s[[i]]))
  # sample.list_leaf_16s[[i]][,1] <- rownames(sample.list_leaf_16s[[i]])
  colnames(df) <- c("ASV_ID", "stat", "p.value")
  df$ASV_ID <- trimws(df$ASV_ID)
  df$stat <- trimws(df$stat)
  df$p.value <- trimws(df$p.value)
  write.csv(df, paste0("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/",name_16s_leaf[i],"_indic_leaf_16s_sp.csv"), row.names = F)
}

for(i in 1:5){
  df <- sample.list_bg_16s[[i]]
  df <- as.data.frame(str_split_fixed(df[2:nrow(df),], " 0.",3))
  # rownames(sample.list_bg_16s[[i]]) <- sample.list_bg_16s[[i]][,1]
  # rownames(sample.list_bg_16s[[i]]) <- gsub(" ", "", rownames(sample.list_bg_16s[[i]]))
  # sample.list_bg_16s[[i]][,1] <- rownames(sample.list_bg_16s[[i]])
  colnames(df) <- c("ASV_ID", "stat", "p.value")
  df$ASV_ID <- trimws(df$ASV_ID)
  df$stat <- trimws(df$stat)
  df$p.value <- trimws(df$p.value)
  write.csv(df, paste0("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/",name_16s_bg[i],"_indic_bg_16s.csv"), row.names = F)
}

leaf_street_tree_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/Street.Tree_indic_leaf_16s_sp.csv", header = T)
leaf_urban_interior_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/Urban.Interior_indic_leaf_16s_sp.csv", header = T)
leaf_rural_edge_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/Rural.Edge_indic_leaf_16s_sp.csv", header = T)
leaf_rural_interior_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/Rural.Interior_indic_leaf_16s_sp.csv", header = T)

bg_street_tree_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/Street.Tree_indic_bg_16s.csv", header = T)
bg_urban_edge_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/Urban.Edge_indic_bg_16s.csv", header = T)
bg_urban_interior_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/Urban.Interior_indic_bg_16s.csv", header = T)
bg_rural_edge_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/Rural.Edge_indic_bg_16s.csv", header = T)
bg_rural_interior_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/Rural.Interior_indic_bg_16s.csv", header = T)

leaf_street_tree_16s_merged <- merge(leaf_street_tree_16s, leaf_rare_16s_tax, by.x = "ASV_ID", by.y = "Species", all.x = T)
leaf_urban_interior_16s_merged <- merge(leaf_urban_interior_16s, leaf_rare_16s_tax, by.x = "ASV_ID", by.y = "Species")
leaf_rural_edge_16s_merged <- merge(leaf_rural_edge_16s, leaf_rare_16s_tax, by.x = "ASV_ID", by.y = "Species")
leaf_rural_interior_16s_merged <- merge(leaf_rural_interior_16s, leaf_rare_16s_tax, by.x = "ASV_ID", by.y = "Species")

bg_street_tree_16s_merged <- merge(bg_street_tree_16s, bg_rare_16s_tax, by.x = "ASV_ID", by.y = "Species")
bg_urban_edge_16s_merged <- merge(bg_urban_edge_16s, bg_rare_16s_tax, by.x = "ASV_ID", by.y = "Species")
bg_urban_interior_16s_merged <- merge(bg_urban_interior_16s, bg_rare_16s_tax, by.x = "ASV_ID", by.y = "Species")
bg_rural_edge_16s_merged <- merge(bg_rural_edge_16s, bg_rare_16s_tax, by.x = "ASV_ID", by.y = "Species")
bg_rural_interior_16s_merged <- merge(bg_rural_interior_16s, bg_rare_16s_tax, by.x = "ASV_ID", by.y = "Species")

rownames(leaf_street_tree_16s_merged) <- leaf_street_tree_16s_merged[,1]
rownames(leaf_urban_interior_16s_merged) <- leaf_urban_interior_16s_merged[,1]
rownames(leaf_rural_edge_16s_merged) <- leaf_rural_edge_16s_merged[,1]
rownames(leaf_rural_interior_16s_merged) <- leaf_rural_interior_16s_merged[,1]

rownames(bg_street_tree_16s_merged) <- bg_street_tree_16s_merged[,1]
rownames(bg_urban_edge_16s_merged) <- bg_urban_edge_16s_merged[,1]
rownames(bg_urban_interior_16s_merged) <- bg_urban_interior_16s_merged[,1]
rownames(bg_rural_edge_16s_merged) <- bg_rural_edge_16s_merged[,1]
rownames(bg_rural_interior_16s_merged) <- bg_rural_interior_16s_merged[,1]

leaf_street_tree_16s_merged$Group <- "Street.tree"
leaf_urban_interior_16s_merged$Group <- "Urban.interior"
leaf_rural_edge_16s_merged$Group <- "Rural.edge"
leaf_rural_interior_16s_merged$Group <- "Rural.interior"

bg_street_tree_16s_merged$Group <- "Street.tree"
bg_urban_edge_16s_merged$Group <- "Urban.edge"
bg_urban_interior_16s_merged$Group <- "Urban.interior"
bg_rural_edge_16s_merged$Group <- "Rural.edge"
bg_rural_interior_16s_merged$Group <- "Rural.interior"

leaf_street_tree_16s_merged$Sample.Type <- "Leaf"
leaf_urban_interior_16s_merged$Sample.Type <- "Leaf"
leaf_rural_edge_16s_merged$Sample.Type <- "Leaf"
leaf_rural_interior_16s_merged$Sample.Type <- "Leaf"

bg_street_tree_16s_merged$Sample.Type <- "Belowground"
bg_urban_edge_16s_merged$Sample.Type <- "Belowground"
bg_urban_interior_16s_merged$Sample.Type <- "Belowground"
bg_rural_edge_16s_merged$Sample.Type <- "Belowground"
bg_rural_interior_16s_merged$Sample.Type <- "Belowground"

leaf_street_tree_16s <- leaf_street_tree_16s_merged[,-c(1:2)]
leaf_urban_interior_16s <- leaf_urban_interior_16s_merged[,-c(1:2)]
leaf_rural_edge_16s <- leaf_rural_edge_16s_merged[,-c(1:2)]
leaf_rural_interior_16s <- leaf_rural_interior_16s_merged[,-c(1:2)]

bg_street_tree_16s <- bg_street_tree_16s_merged[,-c(2)]
bg_urban_edge_16s <- bg_urban_edge_16s_merged[,-c(2)]
bg_urban_interior_16s <- bg_urban_interior_16s_merged[,-c(2)]
bg_rural_edge_16s <- bg_rural_edge_16s_merged[,-c(2)]
bg_rural_interior_16s <- bg_rural_interior_16s_merged[,-c(2)]

leaf_16s_indic <- rbind(leaf_street_tree_16s, leaf_urban_interior_16s, leaf_rural_edge_16s, leaf_rural_interior_16s)
bg_16s_indic <- rbind(bg_street_tree_16s, bg_urban_edge_16s, bg_urban_interior_16s, bg_rural_edge_16s, bg_rural_interior_16s)

write.csv(leaf_16s_indic,"/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/leaf_indic_phylogenetic_info_sp.csv")
write.csv(bg_16s_indic,"/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/16S/bg_indic_phylogenetic_info_sp.csv")
