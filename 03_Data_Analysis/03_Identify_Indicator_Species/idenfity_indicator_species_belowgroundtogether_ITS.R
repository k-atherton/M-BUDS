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

leaf_rare_its_tax <- leaf_rare_its[,1:10]
m_rare_its_tax <- m_rare_its[,1:10]
o_rare_its_tax <- o_rare_its[,1:10]
root_rare_its_tax <- root_rare_its[,1:10]

bg_rare_its_tax <- rbind(m_rare_its_tax, o_rare_its_tax, root_rare_its_tax)
bg_rare_its_tax <- distinct(bg_rare_its_tax)
leaf_rare_its_tax <- distinct(leaf_rare_its_tax)

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

# format data
metadata_its_leaf <- metadata_its[which(metadata_its$sample_type == "Leaf"),]
metadata_its_bg <- metadata_its[which(metadata_its$sample_type == "Belowground"),]

asv_leaf_its <- leaf_rare_its[,2:(ncol(leaf_rare_its))]
asv_bg_its <- bg_rare_its[,2:(ncol(bg_rare_its))]

rownames(asv_leaf_its) <- leaf_rare_its$Species
rownames(asv_bg_its) <- bg_rare_its$Species

asv_leaf.t_its <- t(asv_leaf_its)
asv_bg.t_its <- t(asv_bg_its)

asv_leaf.t_its <- asv_leaf.t_its[order(rownames(asv_leaf.t_its)),]
asv_bg.t_its <- asv_bg.t_its[order(rownames(asv_bg.t_its)),]

percent_leaf_its <- asv_leaf_its/mean(colSums(asv_leaf_its)) * 100
percent_bg_its <- asv_bg_its/mean(colSums(asv_bg_its)) * 100

asv_leaf.t_its[is.na(asv_leaf.t_its)] <- 0
asv_bg.t_its[is.na(asv_bg.t_its)] <- 0
# run indicspecies
inv_leaf_its <- multipatt(asv_leaf.t_its, cluster=factor(metadata_its_leaf$tree_pit_type), func = "r.g", duleg=TRUE)
inv_bg_its <- multipatt(asv_bg.t_its, cluster=factor(metadata_its_bg$tree_pit_type), func = "r.g", duleg=TRUE)

# save information
options(max.print=100000)
sink("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/summary.inv_leaf_its_sp.csv")
summary(inv_leaf_its)
sink()

sink("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/summary.inv_bg_its.csv")
summary(inv_bg_its)
sink()

str_leaf_its <- inv_leaf_its$str
write.csv(str_leaf_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/correlation.inv_leaf_its_sp.csv")

str_bg_its <- inv_bg_its$str
write.csv(str_bg_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/correlation.inv_bg_its.csv")

# get indicator asvs of each location group
summary_leaf_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/summary.inv_leaf_its_sp.csv")
rows_leaf_its <- grep(" Group", summary_leaf_its$Multilevel.pattern.analysis)

summary_bg_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/summary.inv_bg_its.csv")
rows_bg_its <- grep(" Group", summary_bg_its$Multilevel.pattern.analysis)

leaf_street_tree_its <- data.frame(summary_leaf_its[(rows_leaf_its[1]+1):(rows_leaf_its[2]-1),])
leaf_urban_edge_its <- data.frame(summary_leaf_its[(rows_leaf_its[2]+1):(rows_leaf_its[3]-1),])
leaf_urban_interior_its <- data.frame(summary_leaf_its[(rows_leaf_its[3]+1):(rows_leaf_its[4]-1),])
leaf_rural_edge_its <- data.frame(summary_leaf_its[(rows_leaf_its[4]+1):(rows_leaf_its[5]-1),])
leaf_rural_interior_its <- data.frame(summary_leaf_its[(rows_leaf_its[5]+1):(nrow(summary_leaf_its)-2),])

bg_street_tree_its <- data.frame(summary_bg_its[(rows_bg_its[1]+1):(rows_bg_its[2]-1),])
bg_urban_edge_its <- data.frame(summary_bg_its[(rows_bg_its[2]+1):(rows_bg_its[3]-1),])
bg_urban_interior_its <- data.frame(summary_bg_its[(rows_bg_its[3]+1):(rows_bg_its[4]-1),])
bg_rural_edge_its <- data.frame(summary_bg_its[(rows_bg_its[4]+1):(rows_bg_its[5]-1),])
bg_rural_interior_its <- data.frame(summary_bg_its[(rows_bg_its[5]+1):(nrow(summary_bg_its)-2),])

name_its <- c("Street.Tree", "Urban.Edge", "Urban.Interior", "Rural.Edge", "Rural.Interior")

sample.list_leaf_its <- list()
sample.list_leaf_its[[1]] <- leaf_street_tree_its
sample.list_leaf_its[[2]] <- leaf_urban_edge_its
sample.list_leaf_its[[3]] <- leaf_urban_interior_its
sample.list_leaf_its[[4]] <- leaf_rural_edge_its
sample.list_leaf_its[[5]] <- leaf_rural_interior_its

sample.list_bg_its <- list()
sample.list_bg_its[[1]] <- bg_street_tree_its
sample.list_bg_its[[2]] <- bg_urban_edge_its
sample.list_bg_its[[3]] <- bg_urban_interior_its
sample.list_bg_its[[4]] <- bg_rural_edge_its
sample.list_bg_its[[5]] <- bg_rural_interior_its

# save
for(i in 1:5){
  df <- sample.list_leaf_its[[i]]
  df <- as.data.frame(str_split_fixed(df[2:nrow(df),], " 0.",3))
  # rownames(sample.list_leaf_its[[i]]) <- sample.list_leaf_its[[i]][,1]
  # rownames(sample.list_leaf_its[[i]]) <- gsub(" ", "", rownames(sample.list_leaf_its[[i]]))
  # sample.list_leaf_its[[i]][,1] <- rownames(sample.list_leaf_its[[i]])
  colnames(df) <- c("ASV_ID", "stat", "p.value")
  df$ASV_ID <- trimws(df$ASV_ID)
  df$stat <- trimws(df$stat)
  df$p.value <- trimws(df$p.value)
  write.csv(df, paste0("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/ITS/",name_its[i],"_indic_leaf_its_sp.csv"), row.names = F)
}

for(i in 1:5){
  df <- sample.list_bg_its[[i]]
  df <- as.data.frame(str_split_fixed(df[2:nrow(df),], " 0.",3))
  # rownames(sample.list_bg_its[[i]]) <- sample.list_bg_its[[i]][,1]
  # rownames(sample.list_bg_its[[i]]) <- gsub(" ", "", rownames(sample.list_bg_its[[i]]))
  # sample.list_bg_its[[i]][,1] <- rownames(sample.list_bg_its[[i]])
  colnames(df) <- c("ASV_ID", "stat", "p.value")
  df$ASV_ID <- trimws(df$ASV_ID)
  df$stat <- trimws(df$stat)
  df$p.value <- trimws(df$p.value)
  write.csv(df, paste0("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/ITS/",name_its[i],"_indic_bg_its.csv"), row.names = F)
}

leaf_street_tree_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/ITS/Street.Tree_indic_leaf_its_sp.csv", header = T)
leaf_urban_edge_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/ITS/Urban.Edge_indic_leaf_its_sp.csv", header = T)
leaf_urban_interior_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/ITS/Urban.Interior_indic_leaf_its_sp.csv", header = T)
leaf_rural_edge_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/ITS/Rural.Edge_indic_leaf_its_sp.csv", header = T)
leaf_rural_interior_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/ITS/Rural.Interior_indic_leaf_its_sp.csv", header = T)

bg_street_tree_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/ITS/Street.Tree_indic_bg_its.csv", header = T)
bg_urban_edge_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/ITS/Urban.Edge_indic_bg_its.csv", header = T)
bg_urban_interior_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/ITS/Urban.Interior_indic_bg_its.csv", header = T)
bg_rural_edge_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/ITS/Rural.Edge_indic_bg_its.csv", header = T)
bg_rural_interior_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/ITS/Rural.Interior_indic_bg_its.csv", header = T)

leaf_street_tree_its_merged <- merge(leaf_street_tree_its, leaf_rare_its_tax, by.x = "ASV_ID", by.y = "Species", all.x = T)
leaf_urban_edge_its_merged <- merge(leaf_urban_edge_its, leaf_rare_its_tax, by.x = "ASV_ID", by.y = "Species")
leaf_urban_interior_its_merged <- merge(leaf_urban_interior_its, leaf_rare_its_tax, by.x = "ASV_ID", by.y = "Species")
leaf_rural_edge_its_merged <- merge(leaf_rural_edge_its, leaf_rare_its_tax, by.x = "ASV_ID", by.y = "Species")
leaf_rural_interior_its_merged <- merge(leaf_rural_interior_its, leaf_rare_its_tax, by.x = "ASV_ID", by.y = "Species")

bg_street_tree_its_merged <- merge(bg_street_tree_its, bg_rare_its_tax, by.x = "ASV_ID", by.y = "Species")
bg_urban_edge_its_merged <- merge(bg_urban_edge_its, bg_rare_its_tax, by.x = "ASV_ID", by.y = "Species")
bg_urban_interior_its_merged <- merge(bg_urban_interior_its, bg_rare_its_tax, by.x = "ASV_ID", by.y = "Species")
bg_rural_edge_its_merged <- merge(bg_rural_edge_its, bg_rare_its_tax, by.x = "ASV_ID", by.y = "Species")
bg_rural_interior_its_merged <- merge(bg_rural_interior_its, bg_rare_its_tax, by.x = "ASV_ID", by.y = "Species")

rownames(leaf_street_tree_its_merged) <- leaf_street_tree_its_merged[,1]
rownames(leaf_urban_edge_its_merged) <- leaf_urban_edge_its_merged[,1]
rownames(leaf_urban_interior_its_merged) <- leaf_urban_interior_its_merged[,1]
rownames(leaf_rural_edge_its_merged) <- leaf_rural_edge_its_merged[,1]
rownames(leaf_rural_interior_its_merged) <- leaf_rural_interior_its_merged[,1]

rownames(bg_street_tree_its_merged) <- bg_street_tree_its_merged[,1]
rownames(bg_urban_edge_its_merged) <- bg_urban_edge_its_merged[,1]
rownames(bg_urban_interior_its_merged) <- bg_urban_interior_its_merged[,1]
rownames(bg_rural_edge_its_merged) <- bg_rural_edge_its_merged[,1]
rownames(bg_rural_interior_its_merged) <- bg_rural_interior_its_merged[,1]

leaf_street_tree_its_merged$Group <- "Street.tree"
leaf_urban_edge_its_merged$Group <- "Urban.edge"
leaf_urban_interior_its_merged$Group <- "Urban.interior"
leaf_rural_edge_its_merged$Group <- "Rural.edge"
leaf_rural_interior_its_merged$Group <- "Rural.interior"

bg_street_tree_its_merged$Group <- "Street.tree"
bg_urban_edge_its_merged$Group <- "Urban.edge"
bg_urban_interior_its_merged$Group <- "Urban.interior"
bg_rural_edge_its_merged$Group <- "Rural.edge"
bg_rural_interior_its_merged$Group <- "Rural.interior"

leaf_street_tree_its_merged$Sample.Type <- "Leaf"
leaf_urban_edge_its_merged$Sample.Type <- "Leaf"
leaf_urban_interior_its_merged$Sample.Type <- "Leaf"
leaf_rural_edge_its_merged$Sample.Type <- "Leaf"
leaf_rural_interior_its_merged$Sample.Type <- "Leaf"

bg_street_tree_its_merged$Sample.Type <- "Belowground"
bg_urban_edge_its_merged$Sample.Type <- "Belowground"
bg_urban_interior_its_merged$Sample.Type <- "Belowground"
bg_rural_edge_its_merged$Sample.Type <- "Belowground"
bg_rural_interior_its_merged$Sample.Type <- "Belowground"

leaf_street_tree_its <- leaf_street_tree_its_merged[,-c(1:2)]
leaf_urban_edge_its <- leaf_urban_edge_its_merged[,-c(1:2)]
leaf_urban_interior_its <- leaf_urban_interior_its_merged[,-c(1:2)]
leaf_rural_edge_its <- leaf_rural_edge_its_merged[,-c(1:2)]
leaf_rural_interior_its <- leaf_rural_interior_its_merged[,-c(1:2)]

bg_street_tree_its <- bg_street_tree_its_merged[,-c(1:2)]
bg_urban_edge_its <- bg_urban_edge_its_merged[,-c(1:2)]
bg_urban_interior_its <- bg_urban_interior_its_merged[,-c(1:2)]
bg_rural_edge_its <- bg_rural_edge_its_merged[,-c(1:2)]
bg_rural_interior_its <- bg_rural_interior_its_merged[,-c(1:2)]

leaf_its_indic <- rbind(leaf_street_tree_its, leaf_urban_edge_its, leaf_urban_interior_its, leaf_rural_edge_its, leaf_rural_interior_its)
bg_its_indic <- rbind(bg_street_tree_its, bg_urban_edge_its, bg_urban_interior_its, bg_rural_edge_its, bg_rural_interior_its)

write.csv(leaf_its_indic,"/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/ITS/leaf_indic_phylogenetic_info_sp.csv")
write.csv(bg_its_indic,"/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/ITS/bg_indic_phylogenetic_info_sp.csv")

# # save sequences
# rep_seqs_its <- readDNAStringSet("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/its_representative_seqs_unique.fasta")
# seq_name_its <- names(rep_seqs_its)
# sequence_its <- paste(rep_seqs_its)
# rep_seqs_its <- data.frame(seq_name_its, sequence_its)
# all_leaf_its <- rbind(leaf_street_tree_its, leaf_urban_edge_its, leaf_urban_interior_its, leaf_rural_edge_its, leaf_rural_interior_its)
# all_o_its <- rbind(o_street_tree_its, o_urban_edge_its, o_urban_interior_its, o_rural_edge_its, o_rural_interior_its)
# all_m_its <- rbind(m_street_tree_its, m_urban_edge_its, m_urban_interior_its, m_rural_edge_its, m_rural_interior_its)
# all_root_its <- rbind(root_street_tree_its, root_urban_edge_its, root_urban_interior_its, root_rural_edge_its, root_rural_interior_its)
# 
# write.csv(all_leaf_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_indic_phylogenetic_info.csv")
# write.csv(all_o_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_indic_phylogenetic_info.csv")
# write.csv(all_m_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_indic_phylogenetic_info.csv")
# write.csv(all_root_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_indic_phylogenetic_info.csv")
# 
# all_leaf_fasta_its <- merge(all_leaf_its, rep_seqs_its, by.x = 0, by.y = "seq_name_its")
# all_o_fasta_its <- merge(all_o_its, rep_seqs_its, by.x = 0, by.y = "seq_name_its")
# all_m_fasta_its <- merge(all_m_its, rep_seqs_its, by.x = 0, by.y = "seq_name_its")
# all_root_fasta_its <- merge(all_root_its, rep_seqs_its, by.x = 0, by.y = "seq_name_its")
# 
# write.fasta(sequences = as.list(all_leaf_fasta_its$sequence_its), names = all_leaf_fasta_its$Row.names, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_indic_seqs.fasta")
# write.fasta(sequences = as.list(all_o_fasta_its$sequence_its), names = all_o_fasta_its$Row.names, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_indic_seqs.fasta")
# write.fasta(sequences = as.list(all_m_fasta_its$sequence_its), names = all_m_fasta_its$Row.names, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_indic_seqs.fasta")
# write.fasta(sequences = as.list(all_root_fasta_its$sequence_its), names = all_root_fasta_its$Row.names, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_indic_seqs.fasta")
# 
# # get 28S sequence
# sample_leaf_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_indic_phylogenetic_info.csv", header = T)
# sample_o_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_indic_phylogenetic_info.csv", header = T)
# sample_m_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_indic_phylogenetic_info.csv", header = T)
# sample_root_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_indic_phylogenetic_info.csv", header = T)
# 
# sample_leaf_its <- data.frame(lapply(sample_leaf_its, function(x){gsub(pattern="o__", replacement = "", x)}), stringsAsFactors = FALSE)
# sample_o_its <- data.frame(lapply(sample_o_its, function(x){gsub(pattern="o__", replacement = "", x)}), stringsAsFactors = FALSE)
# sample_m_its <- data.frame(lapply(sample_m_its, function(x){gsub(pattern="o__", replacement = "", x)}), stringsAsFactors = FALSE)
# sample_root_its <- data.frame(lapply(sample_root_its, function(x){gsub(pattern="o__", replacement = "", x)}), stringsAsFactors = FALSE)
# 
# colnames(sample_leaf_its)[1] <- "ASV_ID"
# colnames(sample_o_its)[1] <- "ASV_ID"
# colnames(sample_m_its)[1] <- "ASV_ID"
# colnames(sample_root_its)[1] <- "ASV_ID"
# 
# database_its <- read.fasta("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/SILVA_138.1_LSURef_NR99_tax_silva.fasta", as.string=TRUE)
# at_its <- data.frame(attributes(database_its))
# 
# # format database
# tax.list_its <- NULL
# id.list_its <- NULL
# seq.list_its <- NULL
# for(i in 1:nrow(at_its)){
#   description <- attributes(database_its[[at_its[i,]]])[["Annot"]]
#   annot <- strsplit(description, " +")[[1]][2]
#   if(grepl("Fungi", annot, fixed = TRUE)){
#     id <- at_its[i,]
#     id.list_its <- c(id.list_its, id)
#     seq <- database_its[[at_its[i,]]]
#     seq.list_its <- rbind(seq.list_its, seq)
#     tax.name <- strsplit(annot, ";")[[1]]
#     tax.name <- c(tax.name, rep(NA, (14-length(tax.name))))
#     tax.list_its <- rbind(tax.list_its, tax.name)}
# }
# tax.data_its <- data.frame(tax.list_its)
# seq.data_its <- data.frame(seq.list_its)
# 
# database_its <- tax.data_its[,6:14]
# database_its$ID <- id.list_its
# database_its$seq <- seq.data_its$seq.list_its
# database_its <- data.frame(database_its)
# 
# colnames(database_its) <- c("Kingdom","Subkingdom","Phylum","Subphylum","Class","Order","Family","Genus","Species","ID","Sequence")
# 
# # assign 28S sequences to ASVs referring to Genus
# database_its.genus <- na.omit(cbind(database_its$Genus, database_its$Sequence))
# colnames(database_its.genus) <- c("Genus", "Sequence")
# database_its.genus <- data.frame(database_its.genus)
# db_its.genus <- database_its.genus[!duplicated(database_its.genus$Genus),]
# 
# merge_leaf_its.genus <- merge(sample_leaf_its, db_its.genus, by.x = "GENUS", by.y = "Genus", all.x = TRUE)
# merge_leaf_its.genus <- merge_leaf_its.genus[order(merge_leaf_its.genus$ASV_ID, decreasing = FALSE),]
# merge_o_its.genus <- merge(sample_o_its, db_its.genus, by.x = "GENUS", by.y = "Genus", all.x = TRUE)
# merge_o_its.genus <- merge_o_its.genus[order(merge_o_its.genus$ASV_ID, decreasing = FALSE),]
# merge_m_its.genus <- merge(sample_m_its, db_its.genus, by.x = "GENUS", by.y = "Genus", all.x = TRUE)
# merge_m_its.genus <- merge_m_its.genus[order(merge_m_its.genus$ASV_ID, decreasing = FALSE),]
# merge_root_its.genus <- merge(sample_root_its, db_its.genus, by.x = "GENUS", by.y = "Genus", all.x = TRUE)
# merge_root_its.genus <- merge_root_its.genus[order(merge_root_its.genus$ASV_ID, decreasing = FALSE),]
# 
# # assign 28S sequences to ASVs referring to Family
# database_its.family <- na.omit(cbind(database_its$Family, database_its$Sequence))
# colnames(database_its.family) <- c("Family", "Sequence")
# database_its.family <- data.frame(database_its.family)
# db_its.family <- database_its.family[!duplicated(database_its.family$Family),]
# 
# merge_leaf_its.family <- merge(sample_leaf_its, db_its.family, by.x = "Family", by.y = "Family", all.x = TRUE)
# merge_leaf_its.family <- merge_leaf_its.family[order(merge_leaf_its.family$ASV_ID, decreasing = FALSE),]
# merge_o_its.family <- merge(sample_o_its, db_its.family, by.x = "Family", by.y = "Family", all.x = TRUE)
# merge_o_its.family <- merge_o_its.family[order(merge_o_its.family$ASV_ID, decreasing = FALSE),]
# merge_m_its.family <- merge(sample_m_its, db_its.family, by.x = "Family", by.y = "Family", all.x = TRUE)
# merge_m_its.family <- merge_m_its.family[order(merge_m_its.family$ASV_ID, decreasing = FALSE),]
# merge_root_its.family <- merge(sample_root_its, db_its.family, by.x = "Family", by.y = "Family", all.x = TRUE)
# merge_root_its.family <- merge_root_its.family[order(merge_root_its.family$ASV_ID, decreasing = FALSE),]
# 
# # assign 28S sequences to ASVs referring to Order
# database_its.order <- na.omit(cbind(database_its$Order, database_its$Sequence))
# colnames(database_its.order) <- c("Order", "Sequence")
# database_its.order <- data.frame(database_its.order)
# db_its.order <- database_its.order[!duplicated(database_its.order$Order),]
# 
# merge_leaf_its.order <- merge(sample_leaf_its, db_its.order, by.x = "Order", by.y = "Order", all.x = TRUE)
# merge_leaf_its.order <- merge_leaf_its.order[order(merge_leaf_its.order$ASV_ID, decreasing = FALSE),]
# merge_o_its.order <- merge(sample_o_its, db_its.order, by.x = "Order", by.y = "Order", all.x = TRUE)
# merge_o_its.order <- merge_o_its.order[order(merge_o_its.order$ASV_ID, decreasing = FALSE),]
# merge_m_its.order <- merge(sample_m_its, db_its.order, by.x = "Order", by.y = "Order", all.x = TRUE)
# merge_m_its.order <- merge_m_its.order[order(merge_m_its.order$ASV_ID, decreasing = FALSE),]
# merge_root_its.order <- merge(sample_root_its, db_its.order, by.x = "Order", by.y = "Order", all.x = TRUE)
# merge_root_its.order <- merge_root_its.order[order(merge_root_its.order$ASV_ID, decreasing = FALSE),]
# 
# # combine taxonomy information
# data_leaf_its <- sample_leaf_its[order(sample_leaf_its$ASV_ID, decreasing = FALSE),]
# data_leaf_its$seqs <- merge_leaf_its.genus$Sequence
# data_o_its <- sample_o_its[order(sample_o_its$ASV_ID, decreasing = FALSE),]
# data_o_its$seqs <- merge_o_its.genus$Sequence
# data_m_its <- sample_m_its[order(sample_m_its$ASV_ID, decreasing = FALSE),]
# data_m_its$seqs <- merge_m_its.genus$Sequence
# data_root_its <- sample_root_its[order(sample_root_its$ASV_ID, decreasing = FALSE),]
# data_root_its$seqs <- merge_root_its.genus$Sequence
# 
# for(i in 1:nrow(data_leaf_its)){
#   if(is.na(data_leaf_its$seqs[i])){data_leaf_its$seqs[i] <- merge_leaf_its.family$Sequence[i]}
#   if(is.na(data_leaf_its$seqs[i])){data_leaf_its$seqs[i] <- merge_leaf_its.order$Sequence[i]}}
# for(i in 1:nrow(data_o_its)){
#   if(is.na(data_o_its$seqs[i])){data_o_its$seqs[i] <- merge_o_its.family$Sequence[i]}
#   if(is.na(data_o_its$seqs[i])){data_o_its$seqs[i] <- merge_o_its.order$Sequence[i]}}
# for(i in 1:nrow(data_m_its)){
#   if(is.na(data_m_its$seqs[i])){data_m_its$seqs[i] <- merge_m_its.family$Sequence[i]}
#   if(is.na(data_m_its$seqs[i])){data_m_its$seqs[i] <- merge_m_its.order$Sequence[i]}}
# for(i in 1:nrow(data_root_its)){
#   if(is.na(data_root_its$seqs[i])){data_root_its$seqs[i] <- merge_root_its.family$Sequence[i]}
#   if(is.na(data_root_its$seqs[i])){data_root_its$seqs[i] <- merge_root_its.order$Sequence[i]}}
# 
# data_leaf_its <- data_leaf_its[!is.na(data_leaf_its$seqs),]
# data_o_its <- data_o_its[!is.na(data_o_its$seqs),]
# data_m_its <- data_m_its[!is.na(data_m_its$seqs),]
# data_root_its <- data_root_its[!is.na(data_root_its$seqs),]
# 
# rownames(data_leaf_its) <- data_leaf_its$ASV_ID
# rownames(data_o_its) <- data_o_its$ASV_ID
# rownames(data_m_its) <- data_m_its$ASV_ID
# rownames(data_root_its) <- data_root_its$ASV_ID
# 
# write.csv(data_leaf_its[,-1], "/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_phylogenetic_info_28S.csv")
# write.csv(data_o_its[,-1], "/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_phylogenetic_info_28S.csv")
# write.csv(data_m_its[,-1], "/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_phylogenetic_info_28S.csv")
# write.csv(data_root_its[,-1], "/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_phylogenetic_info_28S.csv")
# 
# write.fasta(sequences = as.list(data_leaf_its$seqs), names = data_leaf_its$ASV_ID, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/leaf_seqs_28s.fasta")
# write.fasta(sequences = as.list(data_o_its$seqs), names = data_o_its$ASV_ID, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/o_seqs_28s.fasta")
# write.fasta(sequences = as.list(data_m_its$seqs), names = data_m_its$ASV_ID, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/m_seqs_28s.fasta")
# write.fasta(sequences = as.list(data_root_its$seqs), names = data_root_its$ASV_ID, file.out = "/projectnb/talbot-lab-data/Katies_data/Street_Trees/root_seqs_28s.fasta")
# 
# # to visualize in iTOL
# # make a binary file to show their presence or absence
# data_leaf_its$Street.tree <- 0
# data_leaf_its$Urban.edge <- 0
# data_leaf_its$Urban.interior <- 0
# data_leaf_its$Rural.edge <- 0
# data_leaf_its$Rural.interior <- 0
# 
# data_o_its$Street.tree <- 0
# data_o_its$Urban.edge <- 0
# data_o_its$Urban.interior <- 0
# data_o_its$Rural.edge <- 0
# data_o_its$Rural.interior <- 0
# 
# data_m_its$Street.tree <- 0
# data_m_its$Urban.edge <- 0
# data_m_its$Urban.interior <- 0
# data_m_its$Rural.edge <- 0
# data_m_its$Rural.interior <- 0
# 
# data_root_its$Street.tree <- 0
# data_root_its$Urban.edge <- 0
# data_root_its$Urban.interior <- 0
# data_root_its$Rural.edge <- 0
# data_root_its$Rural.interior <- 0
# 
# for(i in 1:nrow(data_leaf_its)){
#   if(data_leaf_its$Group[i] == "Street.tree"){data_leaf_its$Street.tree[i] <- 1}else{data_leaf_its$Street.tree[i] <- 0}
#   if(data_leaf_its$Group[i] == "Urban.edge"){data_leaf_its$Urban.edge[i] <- 1}else{data_leaf_its$Urban.edge[i] <- 0}
#   if(data_leaf_its$Group[i] == "Urban.interior"){data_leaf_its$Urban.interior[i] <- 1}else{data_leaf_its$Urban.interior[i] <- 0}
#   if(data_leaf_its$Group[i] == "Rural.edge"){data_leaf_its$Rural.edge[i] <- 1}else{data_leaf_its$Rural.edge[i] <- 0}
#   if(data_leaf_its$Group[i] == "Rural.interior"){data_leaf_its$Rural.interior[i] <- 1}else{data_leaf_its$Rural.interior[i] <- 0}
# }
# 
# for(i in 1:nrow(data_o_its)){
#   if(data_o_its$Group[i] == "Street.tree"){data_o_its$Street.tree[i] <- 1}else{data_o_its$Street.tree[i] <- 0}
#   if(data_o_its$Group[i] == "Urban.edge"){data_o_its$Urban.edge[i] <- 1}else{data_o_its$Urban.edge[i] <- 0}
#   if(data_o_its$Group[i] == "Urban.interior"){data_o_its$Urban.interior[i] <- 1}else{data_o_its$Urban.interior[i] <- 0}
#   if(data_o_its$Group[i] == "Rural.edge"){data_o_its$Rural.edge[i] <- 1}else{data_o_its$Rural.edge[i] <- 0}
#   if(data_o_its$Group[i] == "Rural.interior"){data_o_its$Rural.interior[i] <- 1}else{data_o_its$Rural.interior[i] <- 0}
# }
# 
# for(i in 1:nrow(data_m_its)){
#   if(data_m_its$Group[i] == "Street.tree"){data_m_its$Street.tree[i] <- 1}else{data_m_its$Street.tree[i] <- 0}
#   if(data_m_its$Group[i] == "Urban.edge"){data_m_its$Urban.edge[i] <- 1}else{data_m_its$Urban.edge[i] <- 0}
#   if(data_m_its$Group[i] == "Urban.interior"){data_m_its$Urban.interior[i] <- 1}else{data_m_its$Urban.interior[i] <- 0}
#   if(data_m_its$Group[i] == "Rural.edge"){data_m_its$Rural.edge[i] <- 1}else{data_m_its$Rural.edge[i] <- 0}
#   if(data_m_its$Group[i] == "Rural.interior"){data_m_its$Rural.interior[i] <- 1}else{data_m_its$Rural.interior[i] <- 0}
# }
# 
# for(i in 1:nrow(data_root_its)){
#   if(data_root_its$Group[i] == "Street.tree"){data_root_its$Street.tree[i] <- 1}else{data_root_its$Street.tree[i] <- 0}
#   if(data_root_its$Group[i] == "Urban.edge"){data_root_its$Urban.edge[i] <- 1}else{data_root_its$Urban.edge[i] <- 0}
#   if(data_root_its$Group[i] == "Urban.interior"){data_root_its$Urban.interior[i] <- 1}else{data_root_its$Urban.interior[i] <- 0}
#   if(data_root_its$Group[i] == "Rural.edge"){data_root_its$Rural.edge[i] <- 1}else{data_root_its$Rural.edge[i] <- 0}
#   if(data_root_its$Group[i] == "Rural.interior"){data_root_its$Rural.interior[i] <- 1}else{data_root_its$Rural.interior[i] <- 0}
# }
# 
# heatmap_leaf_its <- cbind(data_leaf_its$ASV_ID, data_leaf_its$Street.tree, data_leaf_its$Urban.edge, data_leaf_its$Urban.interior, data_leaf_its$Rural.edge, data_leaf_its$Rural.interior)
# heatmap_o_its <- cbind(data_o_its$ASV_ID, data_o_its$Street.tree, data_o_its$Urban.edge, data_o_its$Urban.interior, data_o_its$Rural.edge, data_o_its$Rural.interior)
# heatmap_m_its <- cbind(data_m_its$ASV_ID, data_m_its$Street.tree, data_m_its$Urban.edge, data_m_its$Urban.interior, data_m_its$Rural.edge, data_m_its$Rural.interior)
# heatmap_root_its <- cbind(data_root_its$ASV_ID, data_root_its$Street.tree, data_root_its$Urban.edge, data_root_its$Urban.interior, data_root_its$Rural.edge, data_root_its$Rural.interior)
# 
# write.csv(heatmap_leaf_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/heatmap_leaf_its_iTOL_28S.csv")
# write.csv(heatmap_o_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/heatmap_o_its_iTOL_28S.csv")
# write.csv(heatmap_m_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/heatmap_m_its_iTOL_28S.csv")
# write.csv(heatmap_root_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Indicator_Species/heatmap_root_its_iTOL_28S.csv")

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