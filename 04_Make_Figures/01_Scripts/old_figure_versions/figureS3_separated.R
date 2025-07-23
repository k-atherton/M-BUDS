### LOAD IN LIBRARIES #########################################################
library(dplyr)
library(ggplot2)
library(nlme)
library(emmeans)
library(ggtext)
library(MuMIn)
library(reshape2)
library(fields)

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
metadata_its$sample_type <- gsub("M Soil", "Soil, 15-30 cm", 
                                 metadata_its$sample_type)

metadata_its$sample_type <- gsub("O Root", "Roots, 0-15 cm", 
                                 metadata_its$sample_type)
metadata_its$sample_type <- gsub("M Root", "Roots, 15-30 cm", 
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
metadata_16s$sample_type <- gsub("M Soil", "Soil, 15-30 cm",
                                 metadata_16s$sample_type)

metadata_16s$sample_type <- gsub("O Root", "Roots, 0-15 cm",
                                 metadata_16s$sample_type)
metadata_16s$sample_type <- gsub("M Root", "Roots, 15-30 cm",
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
metadata_its$perc_dung_saprotroph <- log(metadata_its$perc_dung_saprotroph+1)
metadata_its$perc_wood_saprotroph <- log(metadata_its$perc_wood_saprotroph+1)
metadata_its$perc_litter_saprotroph <- log(metadata_its$perc_litter_saprotroph+1)
metadata_its$perc_soil_saprotroph <- log(metadata_its$perc_soil_saprotroph+1)

metadata_16s$perc_chitinolytic <- log(metadata_16s$perc_chitinolytic+1)
metadata_16s$perc_partial_nitrification <- log(metadata_16s$perc_partial_nitrification+1)
metadata_16s$perc_dissim_nitrate_reduction <- log(metadata_16s$perc_dissim_nitrate_reduction+1)
metadata_16s$perc_denitrification <- log(metadata_16s$perc_denitrification+1)
metadata_16s$perc_methanotroph <- log(metadata_16s$perc_methanotroph+1)

# Make a list of fungal guild names
fungal_trait_names <- c("perc_algal_parasite", "perc_animal_endosymbiont", 
                        "perc_animal_parasite", 
                        "perc_amf", "perc_arthropod_associated", 
                        "perc_dung_saprotroph", "perc_ecm", "perc_epiphyte",
                        "perc_foliar_endophyte", "perc_lichen_parasite",
                        "perc_lichenized", "perc_litter_saprotroph",
                        "perc_moss_symbiont", "perc_mycoparasite",
                        "perc_nectar_tap_saprotroph", 
                        "perc_plant_pathogen", 
                        "perc_pollen_saprotroph", "perc_root_endophyte",
                        "perc_soil_saprotroph", "perc_sooty_mold",
                        "perc_wood_saprotroph", 
                        "perc_animal_pathogen", 
                        "perc_human_pathogen", "perc_opportunistic_pathogen",
                        "perc_root_associated_pathogen", 
                        "perc_leaf_fruit_seed_pathogen", "perc_wood_pathogen",
                        "shannon", "shannon_nopath", "perc_all_saprotroph", 
                        "perc_all_pathotroph", "perc_all_symbiont")

# Make a list of bacterial guild names
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
                           "shannon", "shannon.x", "shannon_nopath")

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE S3A #################################################################
# Select the leaf ITS metadata
metadata_leaf_its <- metadata_its[which(metadata_its$sample_type == 
                                               "Leaf"),]

# Select the leaf samples from the aitchison data
leaf_aitch_its <- all_aitch_its[which(row.names(all_aitch_its) %in%
                                         metadata_leaf_its$sample_name),
                                 which(colnames(all_aitch_its) %in%
                                         metadata_leaf_its$sample_name)]

# Run PCA
pca_result_leaf_its <- prcomp(leaf_aitch_its, scale. = TRUE)

# Make a dataframe of the PCs
pca_leaf_df_its <- data.frame(Sample = rownames(leaf_aitch_its), 
                              PC1 = pca_result_leaf_its$x[,1], 
                              PC2 = pca_result_leaf_its$x[,2])

# Merge the PCs with the leaf metadata
pca_leaf_df_its <- merge(pca_leaf_df_its, metadata_leaf_its, 
                         by.x = "Sample", by.y = "sample_name", all.x = T)

# Summarize the leaf PC data
pca_summary_leaf_its <- summary(pca_result_leaf_its)
pca_variance_leaf_its <- pca_summary_leaf_its$importance[2, ]  # Proportion of variance
pca_variance_cumulative_leaf_its <- pca_summary_leaf_its$importance[3, ]  # Cumulative proportion of variance

### STATISTICS LEAF PCA ###
# define the spatial matrix
dec.degrees.its<-as.matrix(cbind(pca_leaf_df_its$tree_longitude, 
                                 pca_leaf_df_its$tree_latitude))
rownames(dec.degrees.its) <-pca_leaf_df_its$Sample 

# find the spatial distance between samples
space.matrix_its<-rdist.earth(dec.degrees.its[,1:2],miles=FALSE)

# define the design matrix
design_its <- as.data.frame(matrix(nrow = nrow(pca_leaf_df_its), ncol = 5, data = 0))
rownames(design_its) <- pca_leaf_df_its$Sample
colnames(design_its) <- c("street_tree", "urban_edge", "urban_interior", 
                          "rural_edge", "rural_interior")
for(i in 1:nrow(design_its)){
  tree_type <- pca_leaf_df_its$tree_pit_type[which(pca_leaf_df_its$Sample == rownames(design_its)[i])]
  if(tree_type == "Street Tree"){
    design_its$street_tree[i] <- 1
  } else if(tree_type == "Urban Edge"){
    design_its$urban_edge[i] <- 1
  } else if(tree_type == "Urban Interior"){
    design_its$urban_interior[i] <- 1
  } else if(tree_type == "Rural Edge"){
    design_its$rural_edge[i] <- 1
  } else if(tree_type == "Rural Interior"){
    design_its$rural_interior[i] <- 1
  }
}

# Find the design distance between samples
design.dist_its <- dist(design_its)

# Test the effect of space and design on the PCA groupings
MRM <- lm(as.numeric(as.dist(leaf_aitch_its)) ~  as.numeric(design.dist_its) + as.numeric(as.dist(space.matrix_its)))
anova(MRM)

# Plot figure s3a
figs3a <- ggplot(pca_leaf_df_its, 
                         aes(x = PC1, y = PC2, color = tree_pit_type)) +
  geom_point(shape = 18, alpha = 0.5, size = 2.5) +
  xlab(paste("Principal Component 1\n(", 
             round(pca_variance_leaf_its[1] * 100, 2), "% variance)", 
             sep = "")) +
  ylab(paste("Principal Component 2\n(", 
             round(pca_variance_leaf_its[2] * 100, 2), "% variance)", 
             sep = "")) +
  theme_classic() + stat_ellipse() + 
  labs(shape = "", color = "") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), legend.position = "bottom") +
  ggtitle("Fungal community in leaves") + 
  scale_color_manual(values = categorical_pal)

figs3a
ggsave(paste0("atherton_figures3a_ITSpcaleaves_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs3a, dpi = 300, width = 6, height = 3.5, units = "in")

rm(list = c("dec.degrees.its", "design_its", "metadata_leaf_its", "MRM", 
            "pca_leaf_df_its", "pca_result_leaf_its", "pca_summary_leaf_its", 
            "space.matrix_its"))

### FIGURE S3B #################################################################
# Select the leaf ITS metadata
metadata_leaf_16s <- metadata_16s[which(metadata_16s$sample_type == 
                                          "Leaf"),]

# Select the leaf samples from the aitchison data
leaf_aitch_16s <- leaf_aitch_16s[which(row.names(leaf_aitch_16s) %in%
                                         metadata_leaf_16s$sample_name),
                                 which(colnames(leaf_aitch_16s) %in%
                                         metadata_leaf_16s$sample_name)]

# Run PCA
pca_result_leaf_16s <- prcomp(leaf_aitch_16s, scale. = TRUE)

# Make a dataframe of the PCs
pca_leaf_df_16s <- data.frame(Sample = rownames(leaf_aitch_16s), 
                              PC1 = pca_result_leaf_16s$x[,1], 
                              PC2 = pca_result_leaf_16s$x[,2])

# Merge the PCs with the leaf metadata
pca_leaf_df_16s <- merge(pca_leaf_df_16s, metadata_leaf_16s, 
                         by.x = "Sample", by.y = "sample_name", all = T)

# Summarize the leaf PC data
pca_summary_leaf_16s <- summary(pca_result_leaf_16s)
pca_variance_leaf_16s <- pca_summary_leaf_16s$importance[2, ]  # Proportion of variance
pca_variance_cumulative_leaf_16s <- pca_summary_leaf_16s$importance[3, ]  # Cumulative proportion of variance

### STATISTICS LEAF PCA ###
# define the spatial matrix
dec.degrees.16s<-as.matrix(cbind(pca_leaf_df_16s$tree_longitude, 
                                 pca_leaf_df_16s$tree_latitude))
rownames(dec.degrees.16s) <-pca_leaf_df_16s$Sample 

# find the spatial distance between samples
space.matrix_16s<-rdist.earth(dec.degrees.16s[,1:2],miles=FALSE)

# define the design matrix
design_16s <- as.data.frame(matrix(nrow = nrow(pca_leaf_df_16s), ncol = 5, data = 0))
rownames(design_16s) <- pca_leaf_df_16s$Sample
colnames(design_16s) <- c("street_tree", "urban_edge", "urban_interior", 
                          "rural_edge", "rural_interior")
for(i in 1:nrow(design_16s)){
  tree_type <- pca_leaf_df_16s$tree_pit_type[which(pca_leaf_df_16s$Sample == rownames(design_16s)[i])]
  if(tree_type == "Street Tree"){
    design_16s$street_tree[i] <- 1
  } else if(tree_type == "Urban Edge"){
    design_16s$urban_edge[i] <- 1
  } else if(tree_type == "Urban Interior"){
    design_16s$urban_interior[i] <- 1
  } else if(tree_type == "Rural Edge"){
    design_16s$rural_edge[i] <- 1
  } else if(tree_type == "Rural Interior"){
    design_16s$rural_interior[i] <- 1
  }
}

# Find the design distance between samples
design.dist_16s <- dist(design_16s)

# Test the effect of space and design on the PCA groupings
MRM <- lm(as.numeric(as.dist(leaf_aitch_16s)) ~  as.numeric(design.dist_16s) + as.numeric(as.dist(space.matrix_16s)))
anova(MRM)

# Plot figure s3b
figs3b <- ggplot(pca_leaf_df_16s, 
                       aes(x = PC1, y = PC2, color = tree_pit_type)) +
  geom_point(shape = 18, alpha = 0.5, size = 2.5) +
  xlab(paste("Principal Component 1\n(", 
             round(pca_variance_leaf_16s[1] * 100, 2), "% variance)", sep = "")) +
  ylab(paste("Principal Component 2\n(", 
             round(pca_variance_leaf_16s[2] * 100, 2), "% variance)", sep = "")) +
  theme_classic() + stat_ellipse() + 
  labs(shape = "Sample Type", color = "Tree Type") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), legend.position = "none") +
  ggtitle("Bacterial community in leaves") + 
  scale_color_manual(values = categorical_pal)

figs3b
ggsave(paste0("atherton_figures3b_16Spcaleaves_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs3a, dpi = 300, width = 6, height = 3.5, units = "in")

rm(list = c("dec.degrees.16s", "design_16s", "metadata_leaf_16s", "MRM", 
            "pca_leaf_df_16s", "pca_result_leaf_16s", "pca_summary_leaf_16s", 
            "space.matrix_16s"))
### FIGURE S3C #################################################################
# Select the belowground ITS metadata
metadata_osoil_its <- metadata_its[which(metadata_its$sample_type == "Soil, 0-15 cm"),]
metadata_msoil_its <- metadata_its[which(metadata_its$sample_type == "Soil, 15-30 cm"),]
metadata_root_its <- metadata_its[which(metadata_its$sample_type %in% 
                                          c("Roots, 0-15 cm", 
                                            "Roots, 15-30 cm")),]

# Select the belowground samples from the aitchison data
osoil_aitch_its <- all_aitch_its[which(row.names(all_aitch_its) %in%
                                        metadata_osoil_its$sample_name),
                                which(colnames(all_aitch_its) %in%
                                        metadata_osoil_its$sample_name)]
msoil_aitch_its <- all_aitch_its[which(row.names(all_aitch_its) %in%
                                         metadata_msoil_its$sample_name),
                                 which(colnames(all_aitch_its) %in%
                                         metadata_msoil_its$sample_name)]
root_aitch_its <- all_aitch_its[which(row.names(all_aitch_its) %in%
                                         metadata_root_its$sample_name),
                                 which(colnames(all_aitch_its) %in%
                                         metadata_root_its$sample_name)]

# Run PCA
pca_result_osoil_its <- prcomp(osoil_aitch_its, scale. = TRUE)
pca_result_msoil_its <- prcomp(msoil_aitch_its, scale. = TRUE)
pca_result_root_its <- prcomp(root_aitch_its, scale. = TRUE)

# Make a dataframe of the PCs
pca_osoil_df_its <- data.frame(Sample = rownames(osoil_aitch_its), 
                              PC1 = pca_result_osoil_its$x[,1], 
                              PC2 = pca_result_osoil_its$x[,2])
pca_msoil_df_its <- data.frame(Sample = rownames(msoil_aitch_its), 
                               PC1 = pca_result_msoil_its$x[,1], 
                               PC2 = pca_result_msoil_its$x[,2])
pca_root_df_its <- data.frame(Sample = rownames(root_aitch_its), 
                               PC1 = pca_result_root_its$x[,1], 
                               PC2 = pca_result_root_its$x[,2])

# Merge the PCs with the belowground metadata
pca_osoil_df_its <- merge(pca_osoil_df_its, metadata_osoil_its, 
                         by.x = "Sample", by.y = "sample_name", all.x = T)
pca_msoil_df_its <- merge(pca_msoil_df_its, metadata_msoil_its, 
                          by.x = "Sample", by.y = "sample_name", all.x = T)
pca_root_df_its <- merge(pca_root_df_its, metadata_root_its, 
                          by.x = "Sample", by.y = "sample_name", all.x = T)

# Summarize the belowground PC data
pca_summary_osoil_its <- summary(pca_result_osoil_its)
pca_variance_osoil_its <- pca_summary_osoil_its$importance[2, ]  # Proportion of variance
pca_variance_cumulative_osoil_its <- pca_summary_osoil_its$importance[3, ]  # Cumulative proportion of variance

pca_summary_msoil_its <- summary(pca_result_msoil_its)
pca_variance_msoil_its <- pca_summary_msoil_its$importance[2, ]  # Proportion of variance
pca_variance_cumulative_msoil_its <- pca_summary_msoil_its$importance[3, ]  # Cumulative proportion of variance

pca_summary_root_its <- summary(pca_result_root_its)
pca_variance_root_its <- pca_summary_root_its$importance[2, ]  # Proportion of variance
pca_variance_cumulative_root_its <- pca_summary_root_its$importance[3, ]  # Cumulative proportion of variance

### STATISTICS BELOWGROUND PCA ###
# define the spatial matrix
dec.degrees.its_osoil <- as.matrix(cbind(pca_osoil_df_its$tree_longitude, 
                                 pca_osoil_df_its$tree_latitude))
rownames(dec.degrees.its_osoil) <- pca_osoil_df_its$Sample 

dec.degrees.its_msoil <- as.matrix(cbind(pca_msoil_df_its$tree_longitude, 
                                 pca_msoil_df_its$tree_latitude))
rownames(dec.degrees.its_msoil) <- pca_msoil_df_its$Sample 

dec.degrees.its_root <- as.matrix(cbind(pca_root_df_its$tree_longitude, 
                                         pca_root_df_its$tree_latitude))
rownames(dec.degrees.its_root) <- pca_root_df_its$Sample 

# find the spatial distance between samples
space.matrix_its_osoil <-rdist.earth(dec.degrees.its_osoil[,1:2], miles=FALSE)
space.matrix_its_msoil <-rdist.earth(dec.degrees.its_msoil[,1:2], miles=FALSE)
space.matrix_its_root <-rdist.earth(dec.degrees.its_root[,1:2], miles=FALSE)

# define the design matrix
design_its_osoil <- as.data.frame(matrix(nrow = nrow(pca_osoil_df_its), 
                                         ncol = 5, data = 0))
rownames(design_its_osoil) <- pca_osoil_df_its$Sample
colnames(design_its_osoil) <- c("street_tree", "urban_edge", "urban_interior", 
                          "rural_edge", "rural_interior")
for(i in 1:nrow(design_its_osoil)){
  tree_type <- pca_osoil_df_its$tree_pit_type[which(pca_osoil_df_its$Sample == rownames(design_its_osoil)[i])]
  if(tree_type == "Street Tree"){
    design_its_osoil$street_tree[i] <- 1
  } else if(tree_type == "Urban Edge"){
    design_its_osoil$urban_edge[i] <- 1
  } else if(tree_type == "Urban Interior"){
    design_its_osoil$urban_interior[i] <- 1
  } else if(tree_type == "Rural Edge"){
    design_its_osoil$rural_edge[i] <- 1
  } else if(tree_type == "Rural Interior"){
    design_its_osoil$rural_interior[i] <- 1
  }
}

design_its_msoil <- as.data.frame(matrix(nrow = nrow(pca_msoil_df_its), 
                                         ncol = 5, data = 0))
rownames(design_its_msoil) <- pca_msoil_df_its$Sample
colnames(design_its_msoil) <- c("street_tree", "urban_edge", "urban_interior", 
                                "rural_edge", "rural_interior")
for(i in 1:nrow(design_its_msoil)){
  tree_type <- pca_msoil_df_its$tree_pit_type[which(pca_msoil_df_its$Sample == rownames(design_its_msoil)[i])]
  if(tree_type == "Street Tree"){
    design_its_msoil$street_tree[i] <- 1
  } else if(tree_type == "Urban Edge"){
    design_its_msoil$urban_edge[i] <- 1
  } else if(tree_type == "Urban Interior"){
    design_its_msoil$urban_interior[i] <- 1
  } else if(tree_type == "Rural Edge"){
    design_its_msoil$rural_edge[i] <- 1
  } else if(tree_type == "Rural Interior"){
    design_its_msoil$rural_interior[i] <- 1
  }
}

design_its_root <- as.data.frame(matrix(nrow = nrow(pca_root_df_its), 
                                         ncol = 5, data = 0))
rownames(design_its_root) <- pca_root_df_its$Sample
colnames(design_its_root) <- c("street_tree", "urban_edge", "urban_interior", 
                                "rural_edge", "rural_interior")
for(i in 1:nrow(design_its_root)){
  tree_type <- pca_root_df_its$tree_pit_type[which(pca_root_df_its$Sample == rownames(design_its_root)[i])]
  if(tree_type == "Street Tree"){
    design_its_root$street_tree[i] <- 1
  } else if(tree_type == "Urban Edge"){
    design_its_root$urban_edge[i] <- 1
  } else if(tree_type == "Urban Interior"){
    design_its_root$urban_interior[i] <- 1
  } else if(tree_type == "Rural Edge"){
    design_its_root$rural_edge[i] <- 1
  } else if(tree_type == "Rural Interior"){
    design_its_root$rural_interior[i] <- 1
  }
}

# Find the design distance between samples
design.dist_its_osoil <- dist(design_its_osoil)
design.dist_its_msoil <- dist(design_its_msoil)
design.dist_its_root <- dist(design_its_root)

# Test the effect of space and design on the PCA groupings
MRM_osoil <- lm(as.numeric(as.dist(osoil_aitch_its)) ~  as.numeric(design.dist_its_osoil) + as.numeric(as.dist(space.matrix_its_osoil)))
anova(MRM_osoil)

MRM_msoil <- lm(as.numeric(as.dist(msoil_aitch_its)) ~  as.numeric(design.dist_its_msoil) + as.numeric(as.dist(space.matrix_its_msoil)))
anova(MRM_msoil)

MRM_root <- lm(as.numeric(as.dist(root_aitch_its)) ~  as.numeric(design.dist_its_root) + as.numeric(as.dist(space.matrix_its_root)))
anova(MRM_root)

# Plot figure s3c
figs3c_osoil <- ggplot(pca_osoil_df_its, 
                 aes(x = PC1, y = PC2, color = tree_pit_type)) +
  geom_point(aes(shape = sample_type), alpha = 0.5, size = 2) +
  xlab(paste("Principal Component 1 \n(", 
             round(pca_variance_osoil_its[1] * 100, 2), "% variance)", sep = "")) +
  ylab(paste("Principal Component 2 \n(", 
             round(pca_variance_osoil_its[2] * 100, 2), "% variance)", sep = "")) +
  theme_classic() + stat_ellipse() + 
  labs(shape = "Sample Type", color = "Tree Type") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), legend.position = "none") +
  ggtitle("Fungal community in soil, 0-15 cm") + 
  scale_color_manual(values = categorical_pal, guide = "none") + 
  scale_shape_manual(values = shapes)

figs3c_osoil
ggsave(paste0("atherton_figures3c_ITSpcaosoil_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs3c_osoil, dpi = 300, width = 6, height = 3.5, units = "in")

figs3c_msoil <- ggplot(pca_msoil_df_its, 
                       aes(x = PC1, y = PC2, color = tree_pit_type)) +
  geom_point(aes(shape = sample_type), alpha = 0.5, size = 2) +
  xlab(paste("Principal Component 1 \n(", 
             round(pca_variance_msoil_its[1] * 100, 2), "% variance)", sep = "")) +
  ylab(paste("Principal Component 2 \n(", 
             round(pca_variance_msoil_its[2] * 100, 2), "% variance)", sep = "")) +
  theme_classic() + stat_ellipse() + 
  labs(shape = "Sample Type", color = "Tree Type") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), legend.position = "none") +
  ggtitle("Fungal community in soil, 15-30 cm") + 
  scale_color_manual(values = categorical_pal, guide = "none") + 
  scale_shape_manual(values = shapes)

figs3c_msoil
ggsave(paste0("atherton_figures3c_ITSpcamsoil_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs3c_msoil, dpi = 300, width = 6, height = 3.5, units = "in")

figs3c_root <- ggplot(pca_root_df_its, 
                       aes(x = PC1, y = PC2, color = tree_pit_type)) +
  geom_point(aes(shape = sample_type), alpha = 0.5, size = 2) +
  xlab(paste("Principal Component 1 \n(", 
             round(pca_variance_root_its[1] * 100, 2), "% variance)", sep = "")) +
  ylab(paste("Principal Component 2 \n(", 
             round(pca_variance_root_its[2] * 100, 2), "% variance)", sep = "")) +
  theme_classic() + stat_ellipse() + 
  labs(shape = "Sample Type", color = "Tree Type") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), legend.position = "none") +
  ggtitle("Fungal community in roots") + 
  scale_color_manual(values = categorical_pal, guide = "none") + 
  scale_shape_manual(values = shapes)

figs3c_root
ggsave(paste0("atherton_figures3c_ITSpcaroot_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs3c_root, dpi = 300, width = 6, height = 3.5, units = "in")

rm(list = c("dec.degrees.its", "design_its", "metadata_soil_its", "MRM", 
            "pca_soil_df_its", "pca_result_soil_its", "pca_summary_soil_its", 
            "space.matrix_its", "soil_aitch_its"))

### FIGURE S3D ################################################################
# Select the belowground 16s metadata
metadata_osoil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]
metadata_msoil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 15-30 cm"),]
metadata_root_16s <- metadata_16s[which(metadata_16s$sample_type %in% 
                                          c("Roots, 0-15 cm", 
                                            "Roots, 15-30 cm")),]

# Select the belowground samples from the aitchison data
osoil_aitch_16s <- bg_aitch_16s[which(row.names(bg_aitch_16s) %in%
                                         metadata_osoil_16s$sample_name),
                                 which(colnames(bg_aitch_16s) %in%
                                         metadata_osoil_16s$sample_name)]
msoil_aitch_16s <- bg_aitch_16s[which(row.names(bg_aitch_16s) %in%
                                         metadata_msoil_16s$sample_name),
                                 which(colnames(bg_aitch_16s) %in%
                                         metadata_msoil_16s$sample_name)]
root_aitch_16s <- bg_aitch_16s[which(row.names(bg_aitch_16s) %in%
                                        metadata_root_16s$sample_name),
                                which(colnames(bg_aitch_16s) %in%
                                        metadata_root_16s$sample_name)]

# Run PCA
pca_result_osoil_16s <- prcomp(osoil_aitch_16s, scale. = TRUE)
pca_result_msoil_16s <- prcomp(msoil_aitch_16s, scale. = TRUE)
pca_result_root_16s <- prcomp(root_aitch_16s, scale. = TRUE)

# Make a dataframe of the PCs
pca_osoil_df_16s <- data.frame(Sample = rownames(osoil_aitch_16s), 
                               PC1 = pca_result_osoil_16s$x[,1], 
                               PC2 = pca_result_osoil_16s$x[,2])
pca_msoil_df_16s <- data.frame(Sample = rownames(msoil_aitch_16s), 
                               PC1 = pca_result_msoil_16s$x[,1], 
                               PC2 = pca_result_msoil_16s$x[,2])
pca_root_df_16s <- data.frame(Sample = rownames(root_aitch_16s), 
                              PC1 = pca_result_root_16s$x[,1], 
                              PC2 = pca_result_root_16s$x[,2])

# Merge the PCs with the belowground metadata
pca_osoil_df_16s <- merge(pca_osoil_df_16s, metadata_osoil_16s, 
                          by.x = "Sample", by.y = "sample_name", all.x = T)
pca_msoil_df_16s <- merge(pca_msoil_df_16s, metadata_msoil_16s, 
                          by.x = "Sample", by.y = "sample_name", all.x = T)
pca_root_df_16s <- merge(pca_root_df_16s, metadata_root_16s, 
                         by.x = "Sample", by.y = "sample_name", all.x = T)

# Summarize the belowground PC data
pca_summary_osoil_16s <- summary(pca_result_osoil_16s)
pca_variance_osoil_16s <- pca_summary_osoil_16s$importance[2, ]  # Proportion of variance
pca_variance_cumulative_osoil_16s <- pca_summary_osoil_16s$importance[3, ]  # Cumulative proportion of variance

pca_summary_msoil_16s <- summary(pca_result_msoil_16s)
pca_variance_msoil_16s <- pca_summary_msoil_16s$importance[2, ]  # Proportion of variance
pca_variance_cumulative_msoil_16s <- pca_summary_msoil_16s$importance[3, ]  # Cumulative proportion of variance

pca_summary_root_16s <- summary(pca_result_root_16s)
pca_variance_root_16s <- pca_summary_root_16s$importance[2, ]  # Proportion of variance
pca_variance_cumulative_root_16s <- pca_summary_root_16s$importance[3, ]  # Cumulative proportion of variance

### STATISTICS BELOWGROUND PCA ###
# define the spatial matrix
dec.degrees.16s_osoil <- as.matrix(cbind(pca_osoil_df_16s$tree_longitude, 
                                         pca_osoil_df_16s$tree_latitude))
rownames(dec.degrees.16s_osoil) <- pca_osoil_df_16s$Sample 

dec.degrees.16s_msoil <- as.matrix(cbind(pca_msoil_df_16s$tree_longitude, 
                                         pca_msoil_df_16s$tree_latitude))
rownames(dec.degrees.16s_msoil) <- pca_msoil_df_16s$Sample 

dec.degrees.16s_root <- as.matrix(cbind(pca_root_df_16s$tree_longitude, 
                                        pca_root_df_16s$tree_latitude))
rownames(dec.degrees.16s_root) <- pca_root_df_16s$Sample 

# find the spatial distance between samples
space.matrix_16s_osoil <-rdist.earth(dec.degrees.16s_osoil[,1:2], miles=FALSE)
space.matrix_16s_msoil <-rdist.earth(dec.degrees.16s_msoil[,1:2], miles=FALSE)
space.matrix_16s_root <-rdist.earth(dec.degrees.16s_root[,1:2], miles=FALSE)

# define the design matrix
design_16s_osoil <- as.data.frame(matrix(nrow = nrow(pca_osoil_df_16s), 
                                         ncol = 5, data = 0))
rownames(design_16s_osoil) <- pca_osoil_df_16s$Sample
colnames(design_16s_osoil) <- c("street_tree", "urban_edge", "urban_interior", 
                                "rural_edge", "rural_interior")
for(i in 1:nrow(design_16s_osoil)){
  tree_type <- pca_osoil_df_16s$tree_pit_type[which(pca_osoil_df_16s$Sample == rownames(design_16s_osoil)[i])]
  if(tree_type == "Street Tree"){
    design_16s_osoil$street_tree[i] <- 1
  } else if(tree_type == "Urban Edge"){
    design_16s_osoil$urban_edge[i] <- 1
  } else if(tree_type == "Urban Interior"){
    design_16s_osoil$urban_interior[i] <- 1
  } else if(tree_type == "Rural Edge"){
    design_16s_osoil$rural_edge[i] <- 1
  } else if(tree_type == "Rural Interior"){
    design_16s_osoil$rural_interior[i] <- 1
  }
}

design_16s_msoil <- as.data.frame(matrix(nrow = nrow(pca_msoil_df_16s), 
                                         ncol = 5, data = 0))
rownames(design_16s_msoil) <- pca_msoil_df_16s$Sample
colnames(design_16s_msoil) <- c("street_tree", "urban_edge", "urban_interior", 
                                "rural_edge", "rural_interior")
for(i in 1:nrow(design_16s_msoil)){
  tree_type <- pca_msoil_df_16s$tree_pit_type[which(pca_msoil_df_16s$Sample == rownames(design_16s_msoil)[i])]
  if(tree_type == "Street Tree"){
    design_16s_msoil$street_tree[i] <- 1
  } else if(tree_type == "Urban Edge"){
    design_16s_msoil$urban_edge[i] <- 1
  } else if(tree_type == "Urban Interior"){
    design_16s_msoil$urban_interior[i] <- 1
  } else if(tree_type == "Rural Edge"){
    design_16s_msoil$rural_edge[i] <- 1
  } else if(tree_type == "Rural Interior"){
    design_16s_msoil$rural_interior[i] <- 1
  }
}

design_16s_root <- as.data.frame(matrix(nrow = nrow(pca_root_df_16s), 
                                        ncol = 5, data = 0))
rownames(design_16s_root) <- pca_root_df_16s$Sample
colnames(design_16s_root) <- c("street_tree", "urban_edge", "urban_interior", 
                               "rural_edge", "rural_interior")
for(i in 1:nrow(design_16s_root)){
  tree_type <- pca_root_df_16s$tree_pit_type[which(pca_root_df_16s$Sample == rownames(design_16s_root)[i])]
  if(tree_type == "Street Tree"){
    design_16s_root$street_tree[i] <- 1
  } else if(tree_type == "Urban Edge"){
    design_16s_root$urban_edge[i] <- 1
  } else if(tree_type == "Urban Interior"){
    design_16s_root$urban_interior[i] <- 1
  } else if(tree_type == "Rural Edge"){
    design_16s_root$rural_edge[i] <- 1
  } else if(tree_type == "Rural Interior"){
    design_16s_root$rural_interior[i] <- 1
  }
}

# Find the design distance between samples
design.dist_16s_osoil <- dist(design_16s_osoil)
design.dist_16s_msoil <- dist(design_16s_msoil)
design.dist_16s_root <- dist(design_16s_root)

# Test the effect of space and design on the PCA groupings
MRM_osoil <- lm(as.numeric(as.dist(osoil_aitch_16s)) ~  as.numeric(design.dist_16s_osoil) + as.numeric(as.dist(space.matrix_16s_osoil)))
anova(MRM_osoil)

MRM_msoil <- lm(as.numeric(as.dist(msoil_aitch_16s)) ~  as.numeric(design.dist_16s_msoil) + as.numeric(as.dist(space.matrix_16s_msoil)))
anova(MRM_msoil)

MRM_root <- lm(as.numeric(as.dist(root_aitch_16s)) ~  as.numeric(design.dist_16s_root) + as.numeric(as.dist(space.matrix_16s_root)))
anova(MRM_root)

# Plot figure s3d
figs3d_osoil <- ggplot(pca_osoil_df_16s, 
                       aes(x = PC1, y = PC2, color = tree_pit_type)) +
  geom_point(aes(shape = sample_type), alpha = 0.5, size = 2) +
  xlab(paste("Principal Component 1 \n(", 
             round(pca_variance_osoil_16s[1] * 100, 2), "% variance)", sep = "")) +
  ylab(paste("Principal Component 2 \n(", 
             round(pca_variance_osoil_16s[2] * 100, 2), "% variance)", sep = "")) +
  theme_classic() + stat_ellipse() + 
  labs(shape = "Sample Type", color = "Tree Type") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), legend.position = "none") +
  ggtitle("Bacterial community in soil, 0-15 cm") + 
  scale_color_manual(values = categorical_pal, guide = "none") + 
  scale_shape_manual(values = shapes)

figs3d_osoil
ggsave(paste0("atherton_figures3d_16spcaosoil_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs3d_osoil, dpi = 300, width = 6, height = 3.5, units = "in")

figs3d_msoil <- ggplot(pca_msoil_df_16s, 
                       aes(x = PC1, y = PC2, color = tree_pit_type)) +
  geom_point(aes(shape = sample_type), alpha = 0.5, size = 2) +
  xlab(paste("Principal Component 1 \n(", 
             round(pca_variance_msoil_16s[1] * 100, 2), "% variance)", sep = "")) +
  ylab(paste("Principal Component 2 \n(", 
             round(pca_variance_msoil_16s[2] * 100, 2), "% variance)", sep = "")) +
  theme_classic() + stat_ellipse() + 
  labs(shape = "Sample Type", color = "Tree Type") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), legend.position = "none") +
  ggtitle("Bacterial community in soil, 15-30 cm") + 
  scale_color_manual(values = categorical_pal, guide = "none") + 
  scale_shape_manual(values = shapes)

figs3d_msoil
ggsave(paste0("atherton_figures3d_16spcamsoil_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs3d_msoil, dpi = 300, width = 6, height = 3.5, units = "in")

figs3d_root <- ggplot(pca_root_df_16s, 
                      aes(x = PC1, y = PC2, color = tree_pit_type)) +
  geom_point(aes(shape = sample_type), alpha = 0.5, size = 2) +
  xlab(paste("Principal Component 1 \n(", 
             round(pca_variance_root_16s[1] * 100, 2), "% variance)", sep = "")) +
  ylab(paste("Principal Component 2 \n(", 
             round(pca_variance_root_16s[2] * 100, 2), "% variance)", sep = "")) +
  theme_classic() + stat_ellipse() + 
  labs(shape = "Sample Type", color = "Tree Type") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), legend.position = "none") +
  ggtitle("Bacterial community in roots") + 
  scale_color_manual(values = categorical_pal, guide = "none") + 
  scale_shape_manual(values = shapes)

figs3d_root
ggsave(paste0("atherton_figures3d_16spcaroot_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs3d_root, dpi = 300, width = 6, height = 3.5, units = "in")
