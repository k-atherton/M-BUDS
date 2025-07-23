### LOAD IN LIBRARIES #########################################################
library(dplyr)
library(ggplot2)
library(nlme)
library(emmeans)
library(ggtext)
library(MuMIn)
library(reshape2)

### SET UP SCRIPT #############################################################
shapes <- c(16, 17, 15, 8)
linear_pal <- c("#669BBC", "#003049")
categorical_pal <- c("#FF8080","#FF0000","#C00000","#800000","#000000")

mean_se <- function(x, mult = 1) {  
  x <- na.omit(x)
  se <- mult * sqrt(var(x) / length(x))
  mean <- mean(x)
  data.frame(y = mean, ymin = mean - se, ymax = mean + se)
}
standard_error <- function(x) sd(x) / sqrt(length(x)) 

# read in metadata
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/01_Collect_Data/01_Sample_Metadata/")
metadata_16s <- read.csv("atherton_sample_metadata_16S_rareto200_20250121.csv")
metadata_16s <- distinct(metadata_16s)
metadata_its <- read.csv("atherton_sample_metadata_ITS_20250125.csv")
metadata_its <- distinct(metadata_its)

# read in sample names of cleaned data
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/02_Clean_Data/05_Transform_Data/ITS/Aitchison_Distance_Tables/")
all_aitch_its <- read.csv("atherton_ITS_allsampletypes_aitchisondistance_20240919.csv", 
                          row.names = 1)
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/02_Clean_Data/05_Transform_Data/16S/Aitchison_Distance_Tables/")
bg_aitch_16s <- read.csv("atherton_16S_allsampletypes_aitchisondistance_20240925.csv", 
                          row.names = 1)
leaf_aitch_16s <- read.csv("atherton_16S_leaf_aitchisondistance_20241112.csv", 
                          row.names = 1)

### FORMAT METADATA ###########################################################
# Combine ITS soil horizon with sample type
metadata_its$sample_type<- paste0(metadata_its$soil_horizon, " ", 
                                  metadata_its$sample_type)
metadata_its$sample_type <- gsub("NA ", "", 
                                 metadata_its$sample_type)

# Replace ITS soil horizon category with depth measurement
metadata_its$sample_type <- gsub("O Soil", "Soil, 0-15 cm", 
                                 metadata_its$sample_type)

# Replace ITS street tree separate sample types with just Street Tree
metadata_its$tree_pit_type <- gsub("Pit", "Street Tree",
                                   metadata_its$tree_pit_type)
metadata_its$tree_pit_type <- gsub("Grass", "Street Tree",
                                   metadata_its$tree_pit_type)

# Filter the sample names for just the samples from the cleaned ITS dataset
metadata_its <- metadata_its[which(metadata_its$sample_name %in% 
                                     colnames(all_aitch_its)),]

# Make tree pit type a factor with the levels order
metadata_its$tree_pit_type <- factor(metadata_its$tree_pit_type,
                                     levels = c("Street Tree", "Urban Edge",
                                                "Urban Interior", "Rural Edge",
                                                "Rural Interior"))

# Format the ITS environmental data
metadata_its$soil_percent_moisture <- gsub("%", "", 
                                           metadata_its$soil_percent_moisture)
metadata_its$soil_percent_organic_matter <- gsub("%", "",
                                                 metadata_its$soil_percent_organic_matter)
metadata_its$soil_percent_moisture <- as.numeric(metadata_its$soil_percent_moisture)
metadata_its$soil_percent_organic_matter <- as.numeric(metadata_its$soil_percent_organic_matter)
metadata_its$tree_distance_from_edge <- as.numeric(metadata_its$tree_distance_from_edge)

# Combine 16S soil horizon with sample type
metadata_16s$sample_type<- paste0(metadata_16s$soil_horizon, " ",
                                      metadata_16s$sample_type)
metadata_16s$sample_type <- gsub("NA ", "", 
                                     metadata_16s$sample_type)

# Replace 16S soil horizon category with depth measurement
metadata_16s$sample_type <- gsub("O Soil", "Soil, 0-15 cm",
                                     metadata_16s$sample_type)

# Replace 16S street tree separate sample types with just Street Tree
metadata_16s$tree_pit_type <- gsub("Pit", "Street Tree",
                                       metadata_16s$tree_pit_type)
metadata_16s$tree_pit_type <- gsub("Grass", "Street Tree",
                                       metadata_16s$tree_pit_type)

# Filter the sample names for just the samples from the cleaned ITS dataset
metadata_16s <- metadata_16s[which(metadata_16s$sample_name %in% 
                                     c(colnames(bg_aitch_16s), 
                                       colnames(leaf_aitch_16s))),]
metadata_16s <- metadata_16s[which(metadata_16s$sequences_dropped == "No"),]

# Make tree pit type a factor with the levels order
metadata_16s$tree_pit_type <- factor(metadata_16s$tree_pit_type,
                                         levels = c("Street Tree", "Urban Edge",
                                                    "Urban Interior", 
                                                    "Rural Edge",
                                                    "Rural Interior"))

# Format the 16S environmental data
metadata_16s$soil_percent_moisture <- gsub("%", "",
                                           metadata_16s$soil_percent_moisture)
metadata_16s$soil_percent_organic_matter <- gsub("%", "",
                                                 metadata_16s$soil_percent_organic_matter)
metadata_16s$soil_percent_moisture <- as.numeric(metadata_16s$soil_percent_moisture)
metadata_16s$soil_percent_organic_matter <- as.numeric(metadata_16s$soil_percent_organic_matter)
metadata_16s$tree_distance_from_edge <- as.numeric(metadata_16s$tree_distance_from_edge)

# Normalize the data to be used in this figure
metadata_its$perc_epiphyte <- log(metadata_its$perc_epiphyte+1)
metadata_its$perc_sooty_mold <- log(metadata_its$perc_sooty_mold+1)
metadata_its$perc_wood_saprotroph <- log(metadata_its$perc_wood_saprotroph+1)

metadata_16s$perc_cellulolytic <- log(metadata_16s$perc_cellulolytic+1)
metadata_16s$perc_lignolytic <- log(metadata_16s$perc_lignolytic+1)
metadata_16s$perc_plant_pathogen <- log(metadata_16s$perc_plant_pathogen+1)

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE 1D #################################################################
# Select soil, top layer data from ITS
metadata_soil_its <- metadata_its[which(metadata_its$sample_type == 
                                          "Soil, 0-15 cm"),]

# Select soil, top layer data from 16S
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == 
                                          "Soil, 0-15 cm"),]

# Select leaf data from ITS
metadata_leaf_its <- metadata_its[which(metadata_its$sample_type == "Leaf"),]

# Select leaf data from 16S
metadata_leaf_16s <- metadata_16s[which(metadata_16s$sample_type == "Leaf"),]

# Take the functional guild abundances part of the metadata from ITS and 16S
func_soil_its <- metadata_soil_its[,colnames(metadata_soil_its) %in%
                                             c("perc_ecm", 
                                               "perc_wood_saprotroph")]
func_soil_16s <- metadata_soil_16s[,colnames(metadata_soil_16s) %in% 
                                     c("perc_cellulolytic", "perc_lignolytic")]

func_leaf_its <- metadata_leaf_its[,colnames(metadata_leaf_its) %in%
                                             c("perc_epiphyte", 
                                               "perc_sooty_mold")]
func_leaf_16s <- as.data.frame(metadata_leaf_16s[,colnames(metadata_leaf_16s) == 
                                     "perc_plant_pathogen"])

# Rename the column names for the figure
colnames(func_soil_its) <- c("ECM","Wood decomposers")
colnames(func_soil_its) <- paste0(colnames(func_soil_its), " in soil")
colnames(func_leaf_its) <- c("Fungal epiphytes", "Sooty molds")
colnames(func_leaf_its) <- paste0(colnames(func_leaf_its), " in leaves")

colnames(func_soil_16s) <- c("Cellulolytic bacteria", "Lignolytic bacteria")
colnames(func_soil_16s) <- paste0(colnames(func_soil_16s), " in soil")
colnames(func_leaf_16s) <- c("Bacterial plant pathogens")
colnames(func_leaf_16s) <- paste0(colnames(func_leaf_16s), " in leaves")

# Keep the environmental factors for the figure
metadata_soil_heatmap_its <- metadata_soil_its[,c(17,18,20,23,30:36,38,49,53:55,59:61,63:65,67:76)]
metadata_soil_heatmap_16s <- metadata_soil_16s[,c(17,18,20,23,30:36,38,49,53:55,59:61,63:65,67:76)]
metadata_leaf_heatmap_its <- metadata_leaf_its[,c(17,18,20,23,30:36,38,49,53:55,59:61,63:65,67:76)]
metadata_leaf_heatmap_16s <- metadata_leaf_16s[,c(17,18,20,23,30:36,38,49,53:55,59:61,63:65,67:76)]

# Rename the columns for the figure
colnames(metadata_soil_heatmap_its) <- c("Tree distance from edge", 
                                         "Tree distance from Boston", 
                                         "Tree DBH", "Tree growth rate",
                                         "Soil temperature", 
                                         "Soil bulk density", 
                                         "Soil pH", "Soil moisture", 
                                         "Soil organic matter", 
                                         "Soil available NH<sub>4</sub><sup>+</sup>",
                                         "Soil available NO<sub>3</sub><sup>-</sup>", 
                                         "Normalized root biomass", 
                                         "Soil soluble salt",
                                         "NO concentration",
                                         "NO<sub>2</sub> concentration",
                                         "NO<sub>x</sub> concentration",
                                         "NH<sub>4</sub><sup>+</sup> annual throughfall",
                                         "NO<sub>3</sub><sup>-</sup> annual throughfall",
                                         "Total inorganic N annual throughfall",
                                         "O<sub>3</sub> concentration",
                                         "Foliar N %",
                                         "Foliar C:N",
                                         "Ca<sup>2+</sup> annual wet deposition",
                                         "Cl<sup>-</sup> annual wet deposition",
                                         "H<sup>+</sup> annual wet deposition",
                                         "K<sup>+</sup> annual wet deposition",
                                         "Mg<sup>2+</sup> annual wet deposition",
                                         "Na<sup>+</sup> annual wet deposition",
                                         "NH<sub>4</sub><sup>+</sup> annual wet deposition",
                                         "NO<sub>3</sub><sup>-</sup> annual wet deposition",
                                         "SO<sub>4</sub><sup>2-</sup> annual wet deposition",
                                         "Total inorganic N annual wet deposition")

# Rename the columns
colnames(metadata_soil_heatmap_16s) <- c("Tree distance from edge", 
                                         "Tree distance from Boston", 
                                         "Tree DBH", "Tree growth rate",
                                         "Soil temperature", 
                                         "Soil bulk density", 
                                         "Soil pH", "Soil moisture", 
                                         "Soil organic matter", 
                                         "Soil available NH<sub>4</sub><sup>+</sup>",
                                         "Soil available NO<sub>3</sub><sup>-</sup>", 
                                         "Normalized root biomass", 
                                         "Soil soluble salt",
                                         "NO concentration",
                                         "NO<sub>2</sub> concentration",
                                         "NO<sub>x</sub> concentration",
                                         "NH<sub>4</sub><sup>+</sup> annual throughfall",
                                         "NO<sub>3</sub><sup>-</sup> annual throughfall",
                                         "Total inorganic N annual throughfall",
                                         "O<sub>3</sub> concentration",
                                         "Foliar N %",
                                         "Foliar C:N",
                                         "Ca<sup>2+</sup> annual wet deposition",
                                         "Cl<sup>-</sup> annual wet deposition",
                                         "H<sup>+</sup> annual wet deposition",
                                         "K<sup>+</sup> annual wet deposition",
                                         "Mg<sup>2+</sup> annual wet deposition",
                                         "Na<sup>+</sup> annual wet deposition",
                                         "NH<sub>4</sub><sup>+</sup> annual wet deposition",
                                         "NO<sub>3</sub><sup>-</sup> annual wet deposition",
                                         "SO<sub>4</sub><sup>2-</sup> annual wet deposition",
                                         "Total inorganic N annual wet deposition")

# Rename the columns
colnames(metadata_leaf_heatmap_its) <- c("Tree distance from edge", 
                                         "Tree distance from Boston", 
                                         "Tree DBH", "Tree growth rate",
                                         "Soil temperature", 
                                         "Soil bulk density", 
                                         "Soil pH", "Soil moisture", 
                                         "Soil organic matter", 
                                         "Soil available NH<sub>4</sub><sup>+</sup>",
                                         "Soil available NO<sub>3</sub><sup>-</sup>", 
                                         "Normalized root biomass", 
                                         "Soil soluble salt",
                                         "NO concentration",
                                         "NO<sub>2</sub> concentration",
                                         "NO<sub>x</sub> concentration",
                                         "NH<sub>4</sub><sup>+</sup> annual throughfall",
                                         "NO<sub>3</sub><sup>-</sup> annual throughfall",
                                         "Total inorganic N annual throughfall",
                                         "O<sub>3</sub> concentration",
                                         "Foliar N %",
                                         "Foliar C:N",
                                         "Ca<sup>2+</sup> annual wet deposition",
                                         "Cl<sup>-</sup> annual wet deposition",
                                         "H<sup>+</sup> annual wet deposition",
                                         "K<sup>+</sup> annual wet deposition",
                                         "Mg<sup>2+</sup> annual wet deposition",
                                         "Na<sup>+</sup> annual wet deposition",
                                         "NH<sub>4</sub><sup>+</sup> annual wet deposition",
                                         "NO<sub>3</sub><sup>-</sup> annual wet deposition",
                                         "SO<sub>4</sub><sup>2-</sup> annual wet deposition",
                                         "Total inorganic N annual wet deposition")

# Rename the columns
colnames(metadata_leaf_heatmap_16s) <- c("Tree distance from edge", 
                                         "Tree distance from Boston", 
                                         "Tree DBH", "Tree growth rate",
                                         "Soil temperature", 
                                         "Soil bulk density", 
                                         "Soil pH", "Soil moisture", 
                                         "Soil organic matter", 
                                         "Soil available NH<sub>4</sub><sup>+</sup>",
                                         "Soil available NO<sub>3</sub><sup>-</sup>", 
                                         "Normalized root biomass", 
                                         "Soil soluble salt",
                                         "NO concentration",
                                         "NO<sub>2</sub> concentration",
                                         "NO<sub>x</sub> concentration",
                                         "NH<sub>4</sub><sup>+</sup> annual throughfall",
                                         "NO<sub>3</sub><sup>-</sup> annual throughfall",
                                         "Total inorganic N annual throughfall",
                                         "O<sub>3</sub> concentration",
                                         "Foliar N %",
                                         "Foliar C:N",
                                         "Ca<sup>2+</sup> annual wet deposition",
                                         "Cl<sup>-</sup> annual wet deposition",
                                         "H<sup>+</sup> annual wet deposition",
                                         "K<sup>+</sup> annual wet deposition",
                                         "Mg<sup>2+</sup> annual wet deposition",
                                         "Na<sup>+</sup> annual wet deposition",
                                         "NH<sub>4</sub><sup>+</sup> annual wet deposition",
                                         "NO<sub>3</sub><sup>-</sup> annual wet deposition",
                                         "SO<sub>4</sub><sup>2-</sup> annual wet deposition",
                                         "Total inorganic N annual wet deposition")

# Make coefficient matrices
# Soil ITS
Result_soil_its <- NULL
for (i in 1:ncol(func_soil_its)){
  estimate_soil_its <- NULL
  for (j in 1:ncol(metadata_soil_heatmap_its)){
    cor_soil_its <- cor.test(func_soil_its[,i],metadata_soil_heatmap_its[,j], method="p")
    estimate_soil_its <- c(estimate_soil_its, cor_soil_its$estimate)
  }
  Result_soil_its <- rbind(Result_soil_its, estimate_soil_its)
}
colnames(Result_soil_its) <- colnames(metadata_soil_heatmap_its)
rownames(Result_soil_its) <- colnames(func_soil_its) 

# Leaf ITS
Result_leaf_its <- NULL
for (i in 1:ncol(func_leaf_its)){
  estimate_leaf_its <- NULL
  for (j in 1:ncol(metadata_leaf_heatmap_its)){
    cor_leaf_its <- cor.test(func_leaf_its[,i],metadata_leaf_heatmap_its[,j], method="p")
    estimate_leaf_its <- c(estimate_leaf_its, cor_leaf_its$estimate)
  }
  Result_leaf_its <- rbind(Result_leaf_its, estimate_leaf_its)
}
colnames(Result_leaf_its) <- colnames(metadata_leaf_heatmap_its)
rownames(Result_leaf_its) <- colnames(func_leaf_its) 

# Soil 16S
Result_soil_16s <- NULL
for (i in 1:ncol(func_soil_16s)){
  estimate_soil_16s <- NULL
  for (j in 1:ncol(metadata_soil_heatmap_16s)){
    cor_soil_16s <- cor.test(func_soil_16s[,i],metadata_soil_heatmap_16s[,j], method="p")
    estimate_soil_16s <- c(estimate_soil_16s, cor_soil_16s$estimate)
  }
  Result_soil_16s <- rbind(Result_soil_16s, estimate_soil_16s)
}
colnames(Result_soil_16s) <- colnames(metadata_soil_heatmap_16s)
rownames(Result_soil_16s) <- colnames(func_soil_16s) 

# Leaf 16S
Result_leaf_16s <- NULL
estimate_leaf_16s <- NULL
for (j in 1:ncol(metadata_leaf_heatmap_16s)){
  cor_leaf_16s <- cor.test(func_leaf_16s[,1],metadata_leaf_heatmap_16s[,j], method="p")
  estimate_leaf_16s <- c(estimate_leaf_16s, cor_leaf_16s$estimate)
}
Result_leaf_16s <- rbind(Result_leaf_16s, estimate_leaf_16s)
colnames(Result_leaf_16s) <- colnames(metadata_leaf_heatmap_16s)
rownames(Result_leaf_16s) <- "Bacterial plant pathogens in leaves" 

# Make Pval matrices
# Soil ITS
Result.Pval_soil_its <- NULL
for (i in 1:ncol(func_soil_its)){
  Pval_soil_its <- NULL
  for (j in 1:ncol(metadata_soil_heatmap_its)){
    cor_soil_its <- cor.test(func_soil_its[,i],metadata_soil_heatmap_its[,j], method="p")
    Pval_soil_its <- c(Pval_soil_its, cor_soil_its$p.value)
  }
  Result.Pval_soil_its <- rbind (Result.Pval_soil_its, Pval_soil_its)
}
colnames(Result.Pval_soil_its) <- colnames(metadata_soil_heatmap_its)
rownames(Result.Pval_soil_its) <- colnames(func_soil_its) 

# Leaf ITS
Result.Pval_leaf_its <- NULL
for (i in 1:ncol(func_leaf_its)){
  Pval_leaf_its <- NULL
  for (j in 1:ncol(metadata_leaf_heatmap_its)){
    cor_leaf_its <- cor.test(func_leaf_its[,i],metadata_leaf_heatmap_its[,j], method="p")
    Pval_leaf_its <- c(Pval_leaf_its, cor_leaf_its$p.value)
  }
  Result.Pval_leaf_its <- rbind (Result.Pval_leaf_its, Pval_leaf_its)
}
colnames(Result.Pval_leaf_its) <- colnames(metadata_leaf_heatmap_its)
rownames(Result.Pval_leaf_its) <- colnames(func_leaf_its) 

# Soil 16S
Result.Pval_soil_16s <- NULL
for (i in 1:ncol(func_soil_16s)){
  Pval_soil_16s <- NULL
  for (j in 1:ncol(metadata_soil_heatmap_16s)){
    cor_soil_16s <- cor.test(func_soil_16s[,i],metadata_soil_heatmap_16s[,j], method="p")
    Pval_soil_16s <- c(Pval_soil_16s, cor_soil_16s$p.value)
  }
  Result.Pval_soil_16s <- rbind (Result.Pval_soil_16s, Pval_soil_16s)
}
colnames(Result.Pval_soil_16s) <- colnames(metadata_soil_heatmap_16s)
rownames(Result.Pval_soil_16s) <- colnames(func_soil_16s) 

# Leaf 16S
Result.Pval_leaf_16s <- NULL
Pval_leaf_16s <- NULL
for (j in 1:ncol(metadata_leaf_heatmap_16s)){
  cor_leaf_16s <- cor.test(func_leaf_16s[,1],metadata_leaf_heatmap_16s[,j], method="p")
  Pval_leaf_16s <- c(Pval_leaf_16s, cor_leaf_16s$p.value)
}
Result.Pval_leaf_16s <- rbind (Result.Pval_leaf_16s, Pval_leaf_16s)
colnames(Result.Pval_leaf_16s) <- colnames(metadata_leaf_heatmap_16s)
rownames(Result.Pval_leaf_16s) <- "Bacterial plant pathogens in leaves"

# Make Asterisk matrices
# Soil ITS
Result.Pval.combined_soil_its <- Result.Pval_soil_its
Result.asterisk_soil_its <- Result.Pval.combined_soil_its
for (k in 1:ncol(Result.Pval.combined_soil_its)){
  for (i in 1:nrow(Result.Pval.combined_soil_its)){
    if(!is.na(Result.Pval.combined_soil_its[i,k])){
      if (Result.Pval.combined_soil_its[i,k] < 0.001/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_soil_its[i,k] <- "***"
      } else if (Result.Pval.combined_soil_its[i,k] < 0.01/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_soil_its[i,k] <- "**"
      } else if (Result.Pval.combined_soil_its[i,k] < 0.05/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_soil_its[i,k] <- "*"
      } else {
        Result.asterisk_soil_its[i,k]<- " "
      }
    } 
  }
}

# Leaf ITS
Result.Pval.combined_leaf_its <- Result.Pval_leaf_its
Result.asterisk_leaf_its <- Result.Pval.combined_leaf_its
for (k in 1:ncol(Result.Pval.combined_leaf_its)){
  for (i in 1:nrow(Result.Pval.combined_leaf_its)){
    if(!is.na(Result.Pval.combined_leaf_its[i,k])){
      if (Result.Pval.combined_leaf_its[i,k] < 0.001/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_leaf_its[i,k] <- "***"
      } else if (Result.Pval.combined_leaf_its[i,k] < 0.01/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_leaf_its[i,k] <- "**"
      } else if (Result.Pval.combined_leaf_its[i,k] < 0.05/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_leaf_its[i,k] <- "*"
      } else {
        Result.asterisk_leaf_its[i,k]<- " "
      }
    } 
  }
}

# Soil 16S
Result.Pval.combined_soil_16s <- Result.Pval_soil_16s
Result.asterisk_soil_16s <- Result.Pval.combined_soil_16s
for (k in 1:ncol(Result.Pval.combined_soil_16s)){
  for (i in 1:nrow(Result.Pval.combined_soil_16s)){
    if(!is.na(Result.Pval.combined_soil_16s[i,k])){
      if (Result.Pval.combined_soil_16s[i,k] < 0.001/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_soil_16s[i,k] <- "***"
      } else if (Result.Pval.combined_soil_16s[i,k] < 0.01/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_soil_16s[i,k] <- "**"
      } else if (Result.Pval.combined_soil_16s[i,k] < 0.05/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_soil_16s[i,k] <- "*"
      } else {
        Result.asterisk_soil_16s[i,k]<- " "
      }
    } 
  }
}

# Leaf 16S
Result.Pval.combined_leaf_16s <- Result.Pval_leaf_16s
Result.asterisk_leaf_16s <- Result.Pval.combined_leaf_16s
for (k in 1:ncol(Result.Pval.combined_leaf_16s)){
  for (i in 1:nrow(Result.Pval.combined_leaf_16s)){
    if(!is.na(Result.Pval.combined_leaf_16s[i,k])){
      if (Result.Pval.combined_leaf_16s[i,k] < 0.001/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_leaf_16s[i,k] <- "***"
      } else if (Result.Pval.combined_leaf_16s[i,k] < 0.01/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_leaf_16s[i,k] <- "**"
      } else if (Result.Pval.combined_leaf_16s[i,k] < 0.05/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_leaf_16s[i,k] <- "*"
      } else {
        Result.asterisk_leaf_16s[i,k]<- " "
      }
    } 
  }
}

# Combine the p values with the correlation values
Result.combined_soil_its <- Result_soil_its
melt_soil_its <- melt(Result.combined_soil_its)
melt.asterisk_soil_its <- melt(Result.asterisk_soil_its)

Result.combined_leaf_its <- Result_leaf_its
melt_leaf_its <- melt(Result.combined_leaf_its)
melt.asterisk_leaf_its <- melt(Result.asterisk_leaf_its)

Result.combined_soil_16s <- Result_soil_16s
melt_soil_16s <- melt(Result.combined_soil_16s)
melt.asterisk_soil_16s <- melt(Result.asterisk_soil_16s)

Result.combined_leaf_16s <- Result_leaf_16s
melt_leaf_16s <- melt(Result.combined_leaf_16s)
melt.asterisk_leaf_16s <- melt(Result.asterisk_leaf_16s)

melt_all <- rbind(melt_soil_its, melt_soil_16s, melt_leaf_its, melt_leaf_16s)
melt_asterisk_all <- rbind(melt.asterisk_soil_its, melt.asterisk_soil_16s, 
                           melt.asterisk_leaf_its, melt.asterisk_leaf_16s)

# Combine the correlation values with the p-values
vdata_its <- cbind(melt_all, melt_asterisk_all$value)

# Rename the columns
colnames(vdata_its) <- c("microbe", "environment", "coef", "Pval")

# Order the microbial guilds
vdata_its$microbe <- factor(vdata_its$microbe, 
                            levels = c("ECM in soil",
                                       "Fungal epiphytes in leaves",
                                       "Wood decomposers in soil", 
                                       "Cellulolytic bacteria in soil", 
                                       "Lignolytic bacteria in soil",
                                       "Sooty molds in leaves",
                                       "Bacterial plant pathogens in leaves"))

# Make an order of the environmental factors based on their correlation with the first microbial group
order_of_env_factors <- vdata_its[which(vdata_its$microbe == "ECM in soil"),]
vdata_its$environment <- factor(vdata_its$environment, 
                                levels = order_of_env_factors$environment[order(order_of_env_factors$coef)])

# plot figure 1d
figure1d <- ggplot(vdata_its) +
  # this makes the heatmap
  geom_tile(aes(x=microbe, y=environment, fill=coef)) +
  # this adds the asterisks
  geom_text(aes(x=microbe, y=environment, label = Pval), size = 3) +
  labs(x="",y="") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.3)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.3)) + 
  scale_fill_gradient2(low = "#003049", mid = "#FDF0D5", high = "#C00000", 
                       name = "Correlation\ncoefficient") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_markdown(size = 8), 
        axis.text.y = element_markdown(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8), 
        plot.title = element_text(size = 8)) + 
  ggtitle("") + coord_fixed()

figure1d
ggsave(paste0("atherton_figure1d_functionalguilds_vs_environment_heatmap_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), plot = figure1d, 
       device = cairo_pdf, width = 5, height = 10, units = "in", dpi = 300)
