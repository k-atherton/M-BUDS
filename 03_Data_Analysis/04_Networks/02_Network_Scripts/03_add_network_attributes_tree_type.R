library(SpiecEasi)
library(Matrix) 
library(igraph)
library(dplyr)

# load data
bacterial_groups <- read.csv("/projectnb/talbot-lab-data/zrwerbin/soil_bacteria_functional_groups/reference_data/bacteria_func_groups.csv")
genus_groups <- bacterial_groups[which(bacterial_groups$Taxonomic.level == "Genus"),]
genus_groups$Taxon <- paste0("g__", genus_groups$Taxon)
genus_groups <- genus_groups[,c(2,4)]

genus_groups <- distinct(genus_groups)
classification <- c()
for(i in 1:length(unique(genus_groups$Taxon))){
  taxa <- unique(genus_groups$Taxon)[i]
  classes <- genus_groups$Classification[which(genus_groups$Taxon == taxa)]
  classification <- c(classification, paste(classes, collapse = ", "))
}
genus_groups <- data.frame(unique(genus_groups$Taxon), classification)
colnames(genus_groups) <- c("Taxon", "Classification")

tax_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/CLR/ITS/all_its_asv_clr_wtax_guild.csv")
tax_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/CLR/16S/leaf_16s_asv_clr_wtax_guild.csv")
path_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Data_Cleaning/Raw_Data/16S/tax_16s_allruns_wpathogen.csv")
path_16s <- path_16s[,c(8,10)]
path_16s <- distinct(path_16s)

path_genera <- as.data.frame(unique(path_16s$Genus))
colnames(path_genera) <- "Genus"
path_genera$pathogen_type <- NA
for(i in 1:nrow(path_genera)){
  genus <- path_genera$Genus[i]
  pathogen_types <- paste(path_16s$pathogen_type[which(path_16s$Genus == genus)], collapse = ",")
  path_genera$pathogen_type[i] <- pathogen_types
}
path_genera$pathogen_type <- gsub("NA,","", path_genera$pathogen_type)
path_genera$pathogen_type <- gsub(",NA","", path_genera$pathogen_type)
path_16s <- path_genera
path_16s$Genus <- gsub(" g__", "g__", path_16s$Genus)
path_16s$pathogen_type <- gsub("NA", NA, path_16s$pathogen_type)

leaf_networks <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/leaf_networks.RDS")
root_networks <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/root_networks.RDS")
m_networks <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/m_networks.RDS")
o_networks <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/o_networks.RDS") 

leaf_its <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/leaf_its_treetype.RDS")
root_its <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/root_its_treetype.RDS")
m_its <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/m_its_treetype.RDS")
o_its <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/o_its_treetype.RDS")

leaf_16s <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/leaf_16s_treetype.RDS")
root_16s <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/root_16s_treetype.RDS")
m_16s <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/m_16s_treetype.RDS")
o_16s <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/o_16s_treetype.RDS")

genera_its <- c(rownames(leaf_its[[1]]), rownames(leaf_its[[2]]), rownames(leaf_its[[3]]), rownames(leaf_its[[4]]), rownames(leaf_its[[5]]), rownames(leaf_its[[6]]))
genera_its <- c(genera_its, rownames(root_its[[1]]), rownames(root_its[[2]]), rownames(root_its[[3]]), rownames(root_its[[4]]), rownames(root_its[[5]]), rownames(root_its[[6]]))
genera_its <- c(genera_its, rownames(m_its[[1]]), rownames(m_its[[2]]), rownames(m_its[[3]]), rownames(m_its[[4]]), rownames(m_its[[5]]), rownames(m_its[[6]]))
genera_its <- c(genera_its, rownames(o_its[[1]]), rownames(o_its[[2]]), rownames(o_its[[3]]), rownames(o_its[[4]]), rownames(o_its[[5]]), rownames(o_its[[6]]))
genera_its <- unique(genera_its)

genera_16s <- c(rownames(leaf_16s[[1]]), rownames(leaf_16s[[2]]), rownames(leaf_16s[[3]]), rownames(leaf_16s[[4]]), rownames(leaf_16s[[5]]), rownames(leaf_16s[[6]]))
genera_16s <- c(genera_16s, rownames(root_16s[[1]]), rownames(root_16s[[2]]), rownames(root_16s[[3]]), rownames(root_16s[[4]]), rownames(root_16s[[5]]), rownames(root_16s[[6]]))
genera_16s <- c(genera_16s, rownames(m_16s[[1]]), rownames(m_16s[[2]]), rownames(m_16s[[3]]), rownames(m_16s[[4]]), rownames(m_16s[[5]]), rownames(m_16s[[6]]))
genera_16s <- c(genera_16s, rownames(o_16s[[1]]), rownames(o_16s[[2]]), rownames(o_16s[[3]]), rownames(o_16s[[4]]), rownames(o_16s[[5]]), rownames(o_16s[[6]]))
genera_16s <- unique(genera_16s)

tax_its <- tax_its[,c(2:3)]
tax_16s <- tax_16s[,c(2,7)]
tax_16s <- distinct(tax_16s)
tax_16s$Genus <- gsub(" g__", "g__", tax_16s$Genus)
tax_16s <- merge(tax_16s, path_16s, by = "Genus", all.x = TRUE)
tax_16s$guild <- tax_16s$pathogen_type

 for(i in 1:nrow(tax_16s)){
  if(is.na(tax_16s$guild[i])){
    g <- tax_16s$Genus[i]
    k <- tax_16s$Kingdom[i]
    if(g %in% genus_groups$Taxon){
      tax_16s$guild[i] <- paste(genus_groups$Classification[which(genus_groups$Taxon == g)], collapse = ",")
    }
    if(k == "d__Archaea"){
      tax_16s$guild[i] <- "Archaea"
    }
  }
}

tax_its <- distinct(tax_its)
tax_16s <- distinct(tax_16s)

tax_its$GENUS <- gsub("g__", "f__", tax_its$GENUS)
tax_16s$Genus <- gsub("g__", "b__", tax_16s$Genus)

genera_16s <- gsub(" b__", "b__", genera_16s)

tax_its <- tax_its[which(tax_its$GENUS %in% genera_its),]
tax_16s <- tax_16s[which(tax_16s$Genus %in% genera_16s),]

tax_16s$guild[which(tax_16s$Genus == "b__uncultured")] <- NA
tax_16s <- distinct(tax_16s)
tax_16s$pathogen_type[which(tax_16s$Kingdom == "d__Archaea")] <- NA

tax_its$primary_lifestyle[is.na(tax_its$primary_lifestyle)] <- "Other Fungi"
tax_its$primary_lifestyle <- gsub("soil_saprotroph", "Soil Saprotroph", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("epiphyte", "Epiphyte", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("plant_pathogen", "Fungal Plant Pathogen", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("ectomycorrhizal", "Ectomycorrhizal", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("unspecified_saprotroph", "Other Saprotroph", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("litter_saprotroph", "Litter Saprotroph", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("sooty_mold", "Sooty Mold", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("foliar_endophyte", "Foliar Endophyte", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("wood_saprotroph", "Wood Saprotroph", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("animal_parasite", "Fungal Animal/Human Pathogen", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("animal_endosymbiont", "Fungal Animal/Human Endosymbiont", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("mycoparasite", "Mycoparasite", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("arbuscular_mycorrhizal", "Arbuscular Mycorrhizal", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("root_endophyte", "Root Endophyte", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("dung_saprotroph", "Dung Saprotroph", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("pollen_saprotroph", "Pollen Saprotroph", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("nectar/tap_saprotroph", "Nectar/Tap Saprotroph", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("lichenized", "Lichenized Fungi", tax_its$primary_lifestyle)
tax_its$primary_lifestyle <- gsub("algal_parasite", "Fungal Algal Parasite", tax_its$primary_lifestyle)

tax_16s$guild[is.na(tax_16s$guild)] <- "Other Bacteria"
tax_16s$guild <- gsub("t__Zoonotic,t__Animal t__Zoonotic t__Plant,t__Zoonotic t__Animal t__Plant,t__Plant t__Animal t__Zoonotic,t__Zoonotic t__Animal,t__Animal,t__Animal t__Zoonotic,t__Plant t__Animal,t__Plant t__Zoonotic t__Animal,t__Plant,t__Animal t__Plant t__Zoonotic", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Animal t__Plant,t__Plant,t__Plant t__Animal,t__Animal t__Zoonotic,t__Animal t__Plant t__Zoonotic", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal t__Zoonotic t__Plant,t__Plant,t__Zoonotic,t__Animal,t__Zoonotic t__Animal,t__Zoonotic t__Plant", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Zoonotic t__Animal,t__Animal,t__Zoonotic,t__Animal t__Zoonotic,t__Zoonotic t__Animal t__Plant", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal t__Zoonotic t__Plant,t__Animal,t__Zoonotic t__Animal,t__Animal t__Zoonotic,t__Zoonotic", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal t__Zoonotic t__Plant,t__Zoonotic t__Plant t__Animal,t__Zoonotic", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Animal t__Zoonotic,t__Zoonotic t__Animal,t__Zoonotic", "Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Animal t__Zoonotic,t__Animal t__Plant t__Zoonotic,t__Plant", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Animal t__Zoonotic,t__Zoonotic,t__Zoonotic t__Animal", "Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal t__Zoonotic,t__Zoonotic,t__Zoonotic t__Animal,t__Animal", "Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Plant,t__Animal t__Plant,t__Plant t__Animal", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Zoonotic,t__Plant,t__Animal t__Zoonotic", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Plant,t__Plant t__Zoonotic,t__Zoonotic,t__Animal", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Plant,t__Zoonotic,t__Plant t__Animal", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal t__Zoonotic,t__Animal,t__Plant t__Animal", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal t__Zoonotic,t__Animal,t__Zoonotic", "Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Zoonotic t__Animal,t__Zoonotic", "Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Zoonotic,t__Animal,t__Animal t__Zoonotic", "Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Animal t__Zoonotic t__Plant", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal t__Zoonotic,t__Zoonotic", "Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Plant,t__Animal t__Zoonotic", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Plant t__Animal,t__Animal,t__Plant", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal t__Zoonotic,t__Animal", "Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Plant,t__Plant t__Animal", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Plant t__Animal,t__Plant,t__Animal", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Zoonotic t__Animal,t__Animal", "Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Zoonotic t__Animal", "Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Plant,t__Plant t__Animal", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Plant t__Animal", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Zoonotic", "Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Zoonotic,t__Animal", "Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Plant,t__Animal", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Animal,t__Plant", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Plant t__Animal", "General Pathogen", tax_16s$guild)

tax_16s$guild <- gsub("t__Animal", "Animal Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Plant", "Plant Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("t__Zoonotic", "Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("N_fixation", "Nitrogen Fixer", tax_16s$guild)
tax_16s$guild <- gsub("Partial_N2", "Nitrogen Fixer", tax_16s$guild)
tax_16s$guild <- gsub("Nitrite oxidation", "Nitrite Oxidizer", tax_16s$guild)
tax_16s$guild <- gsub("Dissim_nitrate_reduction", "Nitrate Reducer", tax_16s$guild)
tax_16s$guild <- gsub("Cellulolytic, Other N-cycling", "Cellulolytic Bacteria", tax_16s$guild)
tax_16s$guild <- gsub("Cellulolytic$", "Cellulolytic Bacteria", tax_16s$guild)
tax_16s$guild <- gsub("Denitrification", "Denitrifier", tax_16s$guild)
tax_16s$guild <- gsub("Hydrocarbon degradation", "Hydrocarbon Degrader", tax_16s$guild)
tax_16s$guild <- gsub("Plant Pathogen,General Pathogen", "General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("Other N-cycling", "Nitrogen Cycler", tax_16s$guild)
tax_16s$guild <- gsub("Chitinolytic", "Chitinolytic Bacteria", tax_16s$guild)
tax_16s$guild <- gsub("Methanotroph, Partial_Nitrification", "Methanotroph", tax_16s$guild)


tax_16s$guild <- gsub("General Pathogen", "Bacterial General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("Animal/Human Pathogen", "Bacterial Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("Plant Pathogen", "Bacterial Plant Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("Animal Pathogen", "Bacterial Animal Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("^Human Pathogen", "Bacterial Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("Bacterial Animal Pathogen,Bacterial Animal Pathogen Bacterial Plant Pathogen,Human Pathogen", "Bacterial General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("Bacterial Animal Pathogen Bacterial Plant Pathogen,Human Pathogen", "Bacterial General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("Bacterial General Pathogen Bacterial Plant Pathogen", "Bacterial General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("Bacterial Animal Pathogen Human Pathogen Bacterial Plant Pathogen", "Bacterial General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("Bacterial General Pathogen,Bacterial Animal Pathogen", "Bacterial General Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("Bacterial Animal/Human Pathogen,Bacterial Animal Pathogen Human Pathogen", "Bacterial Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("Bacterial Human Pathogen Bacterial Animal/Human Pathogen", "Bacterial Animal/Human Pathogen", tax_16s$guild)
tax_16s$guild <- gsub("Bacterial Animal/Human Pathogen", "Bacterial Animal/Human Pathogen", tax_16s$guild)


attributes_cross_genus <- function(se, site_ITS, site_16S, taxonomy_ITS, taxonomy_16S){
  network <- adj2igraph(getRefit(se))
  network <- set_vertex_attr(network, "shape", index = V(network), 
                             c(rep("square",ncol(site_ITS)), 
                               rep('circle',ncol(site_16S))))
  network <- set_vertex_attr(network, "Genus", index = V(network), 
                             c(colnames(site_ITS), 
                               colnames(site_16S)))
  f_genera <- colnames(site_ITS)
  b_genera <- colnames(site_16S)
  for(i in 1:length(f_genera)){
    f_genera[i] <- unique(taxonomy_ITS$primary_lifestyle[which(taxonomy_ITS$GENUS == f_genera[[i]])])
  }
  for(i in 1:length(b_genera)){
    b_genera[i] <- unique(taxonomy_16S$guild[which(taxonomy_16S$Genus == b_genera[[i]])])
  }
  network <- set_vertex_attr(network, "guild", index = V(network),
                             c(f_genera, b_genera))
  V(network)[V(network)$guild == "Fungal Algal Parasite"]$color <- "firebrick4"
  V(network)[V(network)$guild == "Fungal Animal/Human Pathogen"]$color <- "red3"
  V(network)[V(network)$guild == "Mycoparasite"]$color <- "tomato3"
  V(network)[V(network)$guild == "Fungal Plant Pathogen"]$color <- "red"
  V(network)[V(network)$guild == "Sooty Mold"]$color <- "indianred1"
  V(network)[V(network)$guild == "Dung Saprotroph"]$color <- "darkgoldenrod"
  V(network)[V(network)$guild == "Litter Saprotroph"]$color <- "darkgoldenrod1"
  V(network)[V(network)$guild == "Nectar/Tap Saprotroph"]$color <- "gold"
  V(network)[V(network)$guild == "Pollen Saprotroph"]$color <- "yellow2"
  V(network)[V(network)$guild == "Soil Saprotroph"]$color <- "yellow1"
  V(network)[V(network)$guild == "Wood Saprotroph"]$color <- "lightgoldenrod2"
  V(network)[V(network)$guild == "Other Saprotroph"]$color <- "lemonchiffon"
  V(network)[V(network)$guild == "Arbuscular Mycorrhizal"]$color <- "blue4"
  V(network)[V(network)$guild == "Ectomycorrhizal"]$color <- "#3f51b5"
  V(network)[V(network)$guild == "Epiphyte"]$color <- "blue"
  V(network)[V(network)$guild == "Foliar Endophyte"]$color <- "#03a9f4"
  V(network)[V(network)$guild == "Fungal Animal/Human Endosymbiont"]$color <- "#00bcd4"
  V(network)[V(network)$guild == "Lichenized Fungi"]$color <- "lightskyblue"
  V(network)[V(network)$guild == "Root Endophyte"]$color <- "cadetblue1"
  V(network)[V(network)$guild == "Other Fungi"]$color <- "white"
  V(network)[V(network)$guild == "Bacterial Animal Pathogen"]$color <- "deeppink4"
  V(network)[V(network)$guild == "Bacterial Human Pathogen"]$color <- "deeppink2"
  V(network)[V(network)$guild == "Bacterial Plant Pathogen"]$color <- "hotpink"
  V(network)[V(network)$guild == "Bacterial Animal/Human Pathogen"]$color <- "lightpink"
  V(network)[V(network)$guild == "Bacterial General Pathogen"]$color <- "mistyrose"
  V(network)[V(network)$guild == "Cellulolytic Bacteria"]$color <- "orangered"
  V(network)[V(network)$guild == "Chitinolytic Bacteria"]$color <- "darkorange2"
  V(network)[V(network)$guild == "Methanotroph"]$color <- "orange"
  V(network)[V(network)$guild == "Denitrifier"]$color <- "darkgreen"
  V(network)[V(network)$guild == "Nitrate Reducer"]$color <- "chartreuse4"
  V(network)[V(network)$guild == "Nitrite Oxidizer"]$color <- "chartreuse2"
  V(network)[V(network)$guild == "Nitrogen Cycler"]$color <- "darkolivegreen1"
  V(network)[V(network)$guild == "Nitrogen Fixer"]$color <- "palegreen"
  V(network)[V(network)$guild == "Archaea"]$color <- "purple"
  V(network)[V(network)$guild == "Hydrocarbon Degrader"]$color <- "#795548"
  V(network)[V(network)$guild == "Other Bacteria"]$color <- "black"
  secor <- cov2cor(getOptCov(se))
  elist.gl <- summary(triu(secor*getRefit(se), k=1))
  E(network)$wt <- elist.gl$x
  E(network)$color <- ifelse(E(network)$wt < 0,'lightcoral','steelblue')
  return(network)
}

leaf_network_atts <- vector("list", length(leaf_networks))
leaf_tax_index <- c(1:6)
for(i in 1:length(leaf_network_atts)){
  rownames(leaf_16s[[i]]) <- gsub(" b__", "b__", rownames(leaf_16s[[i]]))
  leaf_network_atts[[i]] <- attributes_cross_genus(leaf_networks[[i]],
                                                       t(leaf_its[[i]]),
                                                       t(leaf_16s[[i]]),
                                                       tax_its[which(tax_its$GENUS %in% colnames(t(leaf_its[[i]]))),], 
                                                       tax_16s[which(tax_16s$Genus %in% colnames(t(leaf_16s[[i]]))),])
}
saveRDS(leaf_network_atts, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/leaf_network_atts_treetype.RDS")

root_network_atts <- vector("list", length(root_networks))
for(i in 1:length(root_network_atts)){
  rownames(root_16s[[i]]) <- gsub(" b__", "b__", rownames(root_16s[[i]]))
  root_network_atts[[i]] <- attributes_cross_genus(root_networks[[i]],
                                                t(root_its[[i]]),
                                                t(root_16s[[i]]),
                                                tax_its[which(tax_its$GENUS %in% colnames(t(root_its[[i]]))),], 
                                                tax_16s[which(tax_16s$Genus %in% colnames(t(root_16s[[i]]))),])
}
saveRDS(root_network_atts, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/root_network_atts.RDS")

m_network_atts <- vector("list", length(m_networks))
for(i in 1:length(m_network_atts)){
  rownames(m_16s[[i]]) <- gsub(" b__", "b__", rownames(m_16s[[i]]))
  m_network_atts[[i]] <- attributes_cross_genus(m_networks[[i]],
                                                   t(m_its[[i]]),
                                                   t(m_16s[[i]]),
                                                   tax_its[which(tax_its$GENUS %in% colnames(t(m_its[[i]]))),], 
                                                   tax_16s[which(tax_16s$Genus %in% colnames(t(m_16s[[i]]))),])
}
saveRDS(m_network_atts, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/m_network_atts.RDS")

o_network_atts <- vector("list", length(o_networks))
for(i in 1:length(o_network_atts)){
  rownames(o_16s[[i]]) <- gsub(" b__", "b__", rownames(o_16s[[i]]))
  o_network_atts[[i]] <- attributes_cross_genus(o_networks[[i]],
                                                   t(o_its[[i]]),
                                                   t(o_16s[[i]]),
                                                   tax_its[which(tax_its$GENUS %in% colnames(t(o_its[[i]]))),], 
                                                   tax_16s[which(tax_16s$Genus %in% colnames(t(o_16s[[i]]))),])
}
saveRDS(o_network_atts, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/o_network_atts.RDS")


