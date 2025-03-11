# load data
leaf_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Decontam/ITS/leaf_its_asv_decontam.csv", row.names = 1)
root_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Decontam/ITS/root_its_asv_decontam.csv", row.names = 1)
o_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Decontam/ITS/o_its_asv_decontam.csv", row.names = 1)
m_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Decontam/ITS/m_its_asv_decontam.csv", row.names = 1)

leaf_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Decontam/16S/leaf_16s_asv_decontam.csv", row.names = 1)
root_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Decontam/16S/root_16s_asv_decontam.csv", row.names = 1)
o_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Decontam/16S/o_16s_asv_decontam.csv", row.names = 1)
m_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Decontam/16S/m_16s_asv_decontam.csv", row.names = 1)

metadata_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Data/metadata_its.csv")
metadata_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Data/metadata_16s_200.csv")

metadata_its <- metadata_its[which(metadata_its$sequences_dropped == "No"),]
metadata_16s <- metadata_16s[which(metadata_16s$sequences_dropped == "No"),]

leaf_its <- leaf_its[,which(colnames(leaf_its) %in% metadata_its$sample_name)]
root_its <- root_its[,which(colnames(root_its) %in% metadata_its$sample_name)]
m_its <- m_its[,which(colnames(m_its) %in% metadata_its$sample_name)]
o_its <- o_its[,which(colnames(o_its) %in% metadata_its$sample_name)]

leaf_16s <- leaf_16s[,which(colnames(leaf_16s) %in% metadata_16s$sample_name)]
root_16s <- root_16s[,which(colnames(root_16s) %in% metadata_16s$sample_name)]
m_16s <- m_16s[,which(colnames(m_16s) %in% metadata_16s$sample_name)]
o_16s <- o_16s[,which(colnames(o_16s) %in% metadata_16s$sample_name)]

tax_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Raw_Data/16S/tax_16s_allruns_wpathogen.csv", row.names = 1)
row.names(tax_16s) <- tax_16s$asv_id
tax_16s <- tax_16s[,-1]
tax_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Raw_Data/ITS/tax_its_allruns.csv", row.names = 1)

# merge asvs and taxonomic data
leaf_its_tax <- merge(tax_its, leaf_its, by.x = 0, by.y = 0, all.y = TRUE)
root_its_tax <- merge(tax_its, root_its, by.x = 0, by.y = 0, all.y = TRUE)
m_its_tax <- merge(tax_its, m_its, by.x = 0, by.y = 0, all.y = TRUE)
o_its_tax <- merge(tax_its, o_its, by.x = 0, by.y = 0, all.y = TRUE)

leaf_16s_tax <- merge(tax_16s, leaf_16s, by.x = 0, by.y = 0, all.y = TRUE)
root_16s_tax <- merge(tax_16s, root_16s, by.x = 0, by.y = 0, all.y = TRUE)
m_16s_tax <- merge(tax_16s, m_16s, by.x = 0, by.y = 0, all.y = TRUE)
o_16s_tax <- merge(tax_16s, o_16s, by.x = 0, by.y = 0, all.y = TRUE)

# keep genus labels
leaf_its_tax <- leaf_its_tax[,c(7,9:ncol(leaf_its_tax))]
root_its_tax <- root_its_tax[,c(7,9:ncol(root_its_tax))]
m_its_tax <- m_its_tax[,c(7,9:ncol(m_its_tax))]
o_its_tax <- o_its_tax[,c(7,9:ncol(o_its_tax))]

leaf_16s_tax <- leaf_16s_tax[,c(7,10:ncol(leaf_16s_tax))]
root_16s_tax <- root_16s_tax[,c(7,10:ncol(root_16s_tax))]
m_16s_tax <- m_16s_tax[,c(7,10:ncol(m_16s_tax))]
o_16s_tax <- o_16s_tax[,c(7,10:ncol(o_16s_tax))]

# collapse to genus level
leaf_its_genus <- aggregate(leaf_its_tax, . ~ Genus, FUN = "sum")
root_its_genus <- aggregate(root_its_tax, . ~ Genus, FUN = "sum")
o_its_genus <- aggregate(o_its_tax, . ~ Genus, FUN = "sum")
m_its_genus <- aggregate(m_its_tax, . ~ Genus, FUN = "sum")

leaf_16s_genus <- aggregate(leaf_16s_tax, . ~ Genus, FUN = "sum")
root_16s_genus <- aggregate(root_16s_tax, . ~ Genus, FUN = "sum")
o_16s_genus <- aggregate(o_16s_tax, . ~ Genus, FUN = "sum")
m_16s_genus <- aggregate(m_16s_tax, . ~ Genus, FUN = "sum")

# move genus labels to row names
row.names(leaf_its_genus) <- leaf_its_genus[,1]
row.names(root_its_genus) <- root_its_genus[,1]
row.names(m_its_genus) <- m_its_genus[,1]
row.names(o_its_genus) <- o_its_genus[,1]

row.names(leaf_16s_genus) <- leaf_16s_genus[,1]
row.names(root_16s_genus) <- root_16s_genus[,1]
row.names(m_16s_genus) <- m_16s_genus[,1]
row.names(o_16s_genus) <- o_16s_genus[,1]

leaf_its_genus <- leaf_its_genus[,-1]
root_its_genus <- root_its_genus[,-1]
m_its_genus <- m_its_genus[,-1]
o_its_genus <- o_its_genus[,-1]

leaf_16s_genus <- leaf_16s_genus[,-1]
root_16s_genus <- root_16s_genus[,-1]
m_16s_genus <- m_16s_genus[,-1]
o_16s_genus <- o_16s_genus[,-1]

# save asv table
write.csv(leaf_its_genus, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/leaf_its_genus.csv")
write.csv(root_its_genus, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/root_its_genus.csv")
write.csv(m_its_genus, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/m_its_genus.csv")
write.csv(o_its_genus, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/o_its_genus.csv")

write.csv(leaf_16s_genus, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/leaf_16s_genus.csv")
write.csv(root_16s_genus, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/root_16s_genus.csv")
write.csv(m_16s_genus, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/m_16s_genus.csv")
write.csv(o_16s_genus, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/o_16s_genus.csv")
