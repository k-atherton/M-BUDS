library(igraph)
source('/projectnb2/talbot-lab-data/Katies_data/UNENetworks/UNE-Networks/Functions/08_Calculate_Network_Metrics.R')

leaf_network_atts <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/leaf_network_atts.RDS")
root_network_atts <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/root_network_atts.RDS")
m_network_atts <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/m_network_atts.RDS")
o_network_atts <- readRDS("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/o_network_atts.RDS")

df <- data.frame(matrix(ncol = 13, nrow = 9))
colnames(df) <- c("Sample_Type", "Plot", "Tree_Type", "Edges",
                  "Density","Complexity",
                  "Transitivity","Short_Path",
                  "Diameter","Betweenness","Degree","Betweenness_otu",
                  "Degree_otu")

df$Sample_Type <- c(rep("Leaf",9))
df$Plot <- c("AA01", "BM03", "HF04", "HF06", "HF069C", "HF069D", "HW07", "JP", 
             "SW08")
df$Tree_Type <- c("Urban Forest", "Urban Forest", "Rural Forest", 
                  "Rural Forest", "Rural Forest", "Rural Forest", 
                  "Urban Forest", "Street Tree", "Urban Forest")
metrics_leaf <- df
metrics_leaf <- network_attributes_cross(leaf_network_atts, metrics_leaf)

df <- data.frame(matrix(ncol = 13, nrow = 9))
colnames(df) <- c("Sample_Type", "Plot", "Tree_Type", "Edges",
                  "Density","Complexity",
                  "Transitivity","Short_Path",
                  "Diameter","Betweenness","Degree","Betweenness_otu",
                  "Degree_otu")

df$Sample_Type <- c(rep("Root",9))
df$Plot <- c("AA01", "BM03", "DO/SB", "HF069C", "HW07", "JP", "MC/SE", 
             "SW08", "WR")
df$Tree_Type <- c("Urban Forest", "Urban Forest", "Street Tree", "Rural Forest", 
                  "Urban Forest", "Street Tree", "Street Tree", "Urban Forest", 
                  "Street Tree")
metrics_root <- df
metrics_root <- network_attributes_cross(root_network_atts, metrics_root)

df <- data.frame(matrix(ncol = 13, nrow = 14))
colnames(df) <- c("Sample_Type", "Plot", "Tree_Type", "Edges",
                  "Density","Complexity",
                  "Transitivity","Short_Path",
                  "Diameter","Betweenness","Degree","Betweenness_otu",
                  "Degree_otu")

df$Sample_Type <- c(rep("Soil, 15-30 cm",14))
df$Plot <- c("AA01", "AB/BB", "BM03", "DO/SB", "HF04", "HF06", "HF069C", "HF069D", 
             "HF069E", "HW07", "JP", "MC/SE", "MP", "WR")
df$Tree_Type <- c("Urban Forest", "Street Tree", "Urban Forest", "Street Tree", 
                  "Rural Forest", "Rural Forest", "Rural Forest", "Rural Forest", 
                  "Rural Forest", "Urban Forest", rep("Street Tree",4))
metrics_m <- df
metrics_m <- network_attributes_cross(m_network_atts, metrics_m)

df <- data.frame(matrix(ncol = 13, nrow = 14))
colnames(df) <- c("Sample_Type", "Plot", "Tree_Type", "Edges",
                  "Density","Complexity",
                  "Transitivity","Short_Path",
                  "Diameter","Betweenness","Degree","Betweenness_otu",
                  "Degree_otu")

df$Sample_Type <- c(rep("Soil, 0-15 cm",14))
df$Plot <- c("AA01", "AB/BB", "BM03", "DO/SB", "HF04", "HF069C", "HF069D", 
             "HF069E", "HW07", "JP", "MC/SE", "MP", "SW08", "WR")
df$Tree_Type <- c("Urban Forest", "Street Tree", "Urban Forest", "Street Tree", 
                  "Rural Forest", "Rural Forest", "Rural Forest", 
                  "Rural Forest", "Urban Forest", rep("Street Tree",3), 
                  "Urban Forest", "Street Tree")
metrics_o <- df
metrics_o <- network_attributes_cross(o_network_atts, metrics_o)

metrics <- rbind(metrics_leaf, metrics_root, metrics_o, metrics_m)

all_networks <- c(leaf_network_atts, root_network_atts, o_network_atts, m_network_atts)

metrics$n_ECM <- NA
metrics$ECM_betweenness <- NA
metrics$ECM_degree <- NA
metrics$n_plantpath <- NA
metrics$plantpath_betweenness <- NA
metrics$plantpath_degree <- NA
metrics$n_animalpath <- NA
metrics$animalpath_betweenness <- NA
metrics$animalpath_degree <- NA
metrics$n_epiphyte <- NA
metrics$epiphyte_betweenness <- NA
metrics$epiphyte_degree <- NA
metrics$n_sap <- NA
metrics$sap_betweenness <- NA
metrics$sap_degree <- NA
metrics$n_endophyte <- NA
metrics$endophyte_betweenness <- NA
metrics$endophyte_degree <- NA

for(i in 1:length(all_networks)){
  net <- all_networks[[i]]
  metrics$n_ECM[i] <- length(V(net)[which(V(net)$guild == "Ectomycorrhizal")])
  metrics$ECM_betweenness[i] <- mean(betweenness(net, V(net)[which(V(net)$guild == "Ectomycorrhizal")]))
  metrics$ECM_degree[i] <- mean(degree(net, V(net)[which(V(net)$guild == "Ectomycorrhizal")]))
  metrics$n_plantpath[i] <- length(V(net)[which(V(net)$guild %in% c("Fungal Plant Pathogen", "Sooty Mold", "Bacterial Plant Pathogen", "Bacterial General Pathogen"))])
  metrics$plantpath_betweenness[i] <- mean(betweenness(net, V(net)[which(V(net)$guild %in% c("Fungal Plant Pathogen", "Sooty Mold", "Bacterial Plant Pathogen", "Bacterial General Pathogen"))]))
  metrics$plantpath_degree[i] <- mean(degree(net, V(net)[which(V(net)$guild %in% c("Fungal Plant Pathogen", "Sooty Mold", "Bacterial Plant Pathogen", "Bacterial General Pathogen"))]))
  metrics$n_animalpath[i] <- length(V(net)[which(V(net)$guild %in% c("Fungal Animal/Human Pathogen", "Bacterial Animal Pathogen", "Bacterial Human Pathogen", "Bacterial Animal/Human Pathogen", "Bacterial General Pathogen"))])
  metrics$animalpath_betweenness[i] <- mean(betweenness(net, V(net)[which(V(net)$guild %in% c("Fungal Animal/Human Pathogen", "Bacterial Animal Pathogen", "Bacterial Human Pathogen", "Bacterial Animal/Human Pathogen", "Bacterial General Pathogen"))]))
  metrics$animalpath_degree[i] <- mean(degree(net, V(net)[which(V(net)$guild %in% c("Fungal Animal/Human Pathogen", "Bacterial Animal Pathogen", "Bacterial Human Pathogen", "Bacterial Animal/Human Pathogen", "Bacterial General Pathogen"))]))
  metrics$n_epiphyte[i] <- length(V(net)[which(V(net)$guild == "Epiphyte")])
  metrics$epiphyte_betweenness[i] <- mean(betweenness(net, V(net)[which(V(net)$guild == "Epiphyte")]))
  metrics$epiphyte_degree[i] <- mean(degree(net, V(net)[which(V(net)$guild == "Epiphyte")]))
  metrics$n_sap[i] <- length(V(net)[which(V(net)$guild %in% c("Dung Saprotroph", "Litter Saprotroph", "Pollen Saprotroph", "Soil Saprotroph", "Wood Saprotroph", "Other Saprotroph"))])
  metrics$sap_betweenness[i] <- mean(betweenness(net, V(net)[which(V(net)$guild %in% c("Dung Saprotroph", "Litter Saprotroph", "Pollen Saprotroph", "Soil Saprotroph", "Wood Saprotroph", "Other Saprotroph"))]))
  metrics$sap_degree[i] <- mean(degree(net, V(net)[which(V(net)$guild %in% c("Dung Saprotroph", "Litter Saprotroph", "Pollen Saprotroph", "Soil Saprotroph", "Wood Saprotroph", "Other Saprotroph"))]))
  metrics$n_endophyte[i] <- length(V(net)[which(V(net)$guild == "Foliar Endophyte")])
  metrics$endophyte_betweenness[i] <- mean(betweenness(net, V(net)[which(V(net)$guild == "Foliar Endophyte")]))
  metrics$endophyte_degree[i] <- mean(degree(net, V(net)[which(V(net)$guild == "Foliar Endophyte")]))
}

metrics$Degree_otu <- gsub("c\\(", "", metrics$Degree_otu)
metrics$Degree_otu <- gsub("\\)", "", metrics$Degree_otu)
metrics$Degree_otu <- gsub(",", "&", metrics$Degree_otu)
metrics$Betweenness_otu <- gsub("c\\(", "", metrics$Betweenness_otu)
metrics$Betweenness_otu <- gsub("\\)", "", metrics$Betweenness_otu)
metrics$Betweenness_otu <- gsub(",", "&", metrics$Betweenness_otu)

write.csv(metrics, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Results/network_metrics_plot.csv", row.names = F)
