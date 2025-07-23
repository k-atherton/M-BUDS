### LOAD IN LIBRARIES #########################################################
library(dplyr)
library(ggplot2)
library(nlme)
library(emmeans)
library(ggtext)
library(MuMIn)
library(reshape2)
library(fields)
library(vegan)
library(multcompView)

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

### FIGURE S4A #################################################################
# Select the leaf ITS metadata
metadata_leaf_its <- metadata_its[which(metadata_its$sample_type == 
                                          "Leaf"),]

# Select the leaf samples from the aitchison data
leaf_aitch_its <- all_aitch_its[which(row.names(all_aitch_its) %in%
                                        metadata_leaf_its$sample_name),
                                which(colnames(all_aitch_its) %in%
                                        metadata_leaf_its$sample_name)]

# Order the aitchison data by sample ID
leaf_aitch_its <- leaf_aitch_its[order(row.names(leaf_aitch_its),
                                       colnames(leaf_aitch_its))]

# Run PCA
pca_result_leaf_its <- prcomp(leaf_aitch_its, scale. = TRUE)

# Make a dataframe of the PCs
pca_leaf_df_its <- data.frame(Sample = rownames(leaf_aitch_its), 
                              PC1 = pca_result_leaf_its$x[,1], 
                              PC2 = pca_result_leaf_its$x[,2])

# Merge the PCs with the leaf metadata
pca_leaf_df_its <- merge(pca_leaf_df_its, metadata_leaf_its, 
                         by.x = "Sample", by.y = "sample_name", all.x = T)

# Order the PCA by sample ID
pca_leaf_df_its <- pca_leaf_df_its[order(pca_leaf_df_its$Sample),]

# Calculate the beta dispersion of the ITS leaf PCA
beta_leaf_its <- betadisper(as.dist(leaf_aitch_its), 
                            group = pca_leaf_df_its$tree_pit_type, 
                            type = "centroid")

### STATISTICS BETA DISPERSION ITS LEAF PCA ###
anova(beta_leaf_its)

# Calculate and define the Tukey groups
tukey_its <- TukeyHSD(beta_leaf_its)
tukey_df_its <- as.data.frame(tukey_its$group)
x_its <- tukey_df_its$`p adj`
names(x_its) <- rownames(tukey_df_its)
tukey_groups_its <- multcompLetters(x_its)$Letters

# Make a dataframe of the beta distances
beta_leaf_dist_its <- as.data.frame(beta_leaf_its$distances)

# Merge the distances with the metadata
beta_leaf_df_its <- merge(beta_leaf_dist_its, metadata_its, 
                          by.x = 0, by.y = "sample_name", leaf = T)

# Order the dataframe by the tree pit type
beta_leaf_df_its <- beta_leaf_df_its[order(beta_leaf_df_its$tree_pit_type),]

# Format the tree pit type for plotting
beta_leaf_df_its$tree_pit_type <- gsub("Interior", "Forest Interior", 
                                       beta_leaf_df_its$tree_pit_type)
beta_leaf_df_its$tree_pit_type <- gsub("Edge", "Forest Edge", 
                                       beta_leaf_df_its$tree_pit_type)
beta_leaf_df_its$tree_pit_type <- gsub("Street", "Urban Street", 
                                       beta_leaf_df_its$tree_pit_type)
beta_leaf_df_its$tree_pit_type <- gsub(" ", "\n", 
                                       beta_leaf_df_its$tree_pit_type)
beta_leaf_df_its$tree_pit_type <- factor(beta_leaf_df_its$tree_pit_type, 
                                         levels = c("Urban\nStreet\nTree", 
                                                    "Urban\nForest\nEdge",
                                                    "Urban\nForest\nInterior",
                                                    "Rural\nForest\nEdge",
                                                    "Rural\nForest\nInterior"))

# Define where the Tukey group labels should go
text_positions <- beta_leaf_df_its %>%
  group_by(tree_pit_type) %>%
  summarise(max_value = max(`beta_leaf_its$distances`), .groups = 'drop') %>%
  mutate(text_position = max_value + 0.5) # Add buffer for text
text_positions$label <- c("a", "a", "a", "ab", "b")

figs4a <- ggplot(beta_leaf_df_its, 
                 aes(x = tree_pit_type, y = `beta_leaf_its$distances`, 
                     fill = tree_pit_type)) + 
  stat_summary(aes(group = 1),
               fun = median,
               geom = "line",
               linetype = "dashed",
               linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) + 
  geom_jitter(aes(color = tree_pit_type), alpha = 0.5, 
                              stroke = NA) +
  theme_classic() + xlab("") + 
  ylab("Beta dispersion") + 
  geom_text(data = text_positions, aes(x = tree_pit_type, y = text_position, 
                                       label = label), 
            inherit.aes = FALSE, vjust = -0.5, family = "Arial") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10),
        legend.position = "none") + 
  scale_color_manual(values = categorical_pal) +
  scale_fill_manual(values = categorical_pal) + ylim(c(20, 55)) +
  ggtitle("Fungal community in leaves") 
figs4a

ggsave(paste0("atherton_figures4a_ITSbetadisperleaves_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs4a, dpi = 300, width = 6, height = 3.5, units = "in")

rm(list = c("beta_leaf_df_its", "beta_leaf_dist_its", "beta_leaf_its", 
            "metadata_leaf_its", "pca_leaf_df_its", "pca_result_leaf_its",
            "text_positions", "tukey_df_its", "tukey_its"))

### FIGURE S4B #################################################################
# Select the leaf ITS metadata
metadata_leaf_16s <- metadata_16s[which(metadata_16s$sample_type == 
                                          "Leaf"),]

# Select the leaf samples from the aitchison data
leaf_aitch_16s <- leaf_aitch_16s[which(row.names(leaf_aitch_16s) %in%
                                         metadata_leaf_16s$sample_name),
                                 which(colnames(leaf_aitch_16s) %in%
                                         metadata_leaf_16s$sample_name)]

# Order the aitchison data by the sample ID
leaf_aitch_16s <- leaf_aitch_16s[order(row.names(leaf_aitch_16s),
                                       colnames(leaf_aitch_16s))]

# Run PCA
pca_result_leaf_16s <- prcomp(leaf_aitch_16s, scale. = TRUE)

# Make a dataframe of the PCs
pca_leaf_df_16s <- data.frame(Sample = rownames(leaf_aitch_16s), 
                              PC1 = pca_result_leaf_16s$x[,1], 
                              PC2 = pca_result_leaf_16s$x[,2])

# Merge the PCs with the leaf metadata
pca_leaf_df_16s <- merge(pca_leaf_df_16s, metadata_leaf_16s, 
                         by.x = "Sample", by.y = "sample_name", all = T)

# Order the PC data by the sample ID
pca_leaf_df_16s <- pca_leaf_df_16s[order(pca_leaf_df_16s$Sample),]

# Calculate beta disperson of the PCA
beta_leaf_16s <- betadisper(as.dist(leaf_aitch_16s), 
                            group = pca_leaf_df_16s$tree_pit_type, 
                            type = "centroid")

### STATISTICS BETA DISPERSION 16S LEAF PCA ###
anova(beta_leaf_16s)

# Calculate and define the Tukey groups
tukey_16s <- TukeyHSD(beta_leaf_16s)
tukey_df_16s <- as.data.frame(tukey_16s$group)
x_16s <- tukey_df_16s$`p adj`
names(x_16s) <- rownames(tukey_df_16s)
tukey_groups_16s <- multcompLetters(x_16s)$Letters

# Make a dataframe of the beta distances
beta_leaf_dist_16s <- as.data.frame(beta_leaf_16s$distances)

# Merge the distances with the metadata
beta_leaf_df_16s <- merge(beta_leaf_dist_16s,
                          metadata_16s[which(metadata_16s$sample_name %in% 
                                               colnames(leaf_aitch_16s)),], 
                          by.x = 0, by.y = "sample_name", leaf = T)

# Order the dataframe by the tree pit type
beta_leaf_df_16s <- beta_leaf_df_16s[order(beta_leaf_df_16s$tree_pit_type),]

# Format the tree pit type for plotting
beta_leaf_df_16s$tree_pit_type <- gsub("Interior", "Forest Interior", 
                                       beta_leaf_df_16s$tree_pit_type)
beta_leaf_df_16s$tree_pit_type <- gsub("Edge", "Forest Edge", 
                                       beta_leaf_df_16s$tree_pit_type)
beta_leaf_df_16s$tree_pit_type <- gsub("Street", "Urban Street",
                                       beta_leaf_df_16s$tree_pit_type)
beta_leaf_df_16s$tree_pit_type <- gsub(" ", "\n", 
                                       beta_leaf_df_16s$tree_pit_type)
beta_leaf_df_16s$tree_pit_type <- factor(beta_leaf_df_16s$tree_pit_type, 
                                         levels = c("Urban\nStreet\nTree", 
                                                    "Urban\nForest\nEdge",
                                                    "Urban\nForest\nInterior",
                                                    "Rural\nForest\nEdge",
                                                    "Rural\nForest\nInterior"))

# Define where the Tukey group labels should go
text_positions <- beta_leaf_df_16s %>%
  group_by(tree_pit_type) %>%
  summarise(max_value = max(`beta_leaf_16s$distances`), .groups = 'drop') %>%
  mutate(text_position = max_value + 0.5) # Add buffer for text
text_positions$label <- c("", "", "", "", "")

figs4b <- ggplot(beta_leaf_df_16s,
                 aes(x = tree_pit_type, y = `beta_leaf_16s$distances`, 
                     fill = tree_pit_type)) + 
  stat_summary(aes(group = 1),
               fun = median,
               geom = "line",
               linetype = "dashed",
               linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) + 
  geom_jitter(aes(color = tree_pit_type), alpha = 0.5, 
              stroke = NA) +
  theme_classic() + xlab("") + 
  ylab("Beta dispersion") + 
  geom_text(data = text_positions, aes(x = tree_pit_type, y = text_position, 
                                       label = label), 
            inherit.aes = FALSE, vjust = -0.5) +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10),
        legend.position = "none") + 
  scale_color_manual(values = categorical_pal) +
  scale_fill_manual(values = categorical_pal) + ylim(c(5, 35)) +
  ggtitle("Bacterial community in leaves") 

figs4b
ggsave(paste0("atherton_figures4b_16Sbetadisperleaves_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs4b, dpi = 300, width = 6, height = 3.5, units = "in")

rm(list = c("beta_leaf_df_16s", "beta_leaf_dist_16s", "beta_leaf_16s", 
            "metadata_leaf_16s", "pca_leaf_df_16s", "pca_result_leaf_16s",
            "text_positions", "tukey_df_16s", "tukey_16s"))

### FIGURE S4C #################################################################
# Select the belowground ITS metadata
metadata_osoil_its <- metadata_its[which(metadata_its$sample_type == 
                                          "Soil, 0-15 cm"),]
metadata_msoil_its <- metadata_its[which(metadata_its$sample_type == 
                                           "Soil, 15-30 cm"),]
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

# Order the aitchison data by the sample ID
osoil_aitch_its <- osoil_aitch_its[order(row.names(osoil_aitch_its),
                                       colnames(osoil_aitch_its))]
msoil_aitch_its <- msoil_aitch_its[order(row.names(msoil_aitch_its),
                                         colnames(msoil_aitch_its))]
root_aitch_its <- root_aitch_its[order(row.names(root_aitch_its),
                                         colnames(root_aitch_its))]

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

# Order the PC data by the sample ID
pca_osoil_df_its <- pca_osoil_df_its[order(pca_osoil_df_its$Sample),]
pca_msoil_df_its <- pca_msoil_df_its[order(pca_msoil_df_its$Sample),]
pca_root_df_its <- pca_root_df_its[order(pca_root_df_its$Sample),]

# Calculate beta disperson of the PCA
beta_osoil_its <- betadisper(as.dist(osoil_aitch_its), 
                            group = pca_osoil_df_its$tree_pit_type, 
                            type = "centroid")
beta_msoil_its <- betadisper(as.dist(msoil_aitch_its), 
                            group = pca_msoil_df_its$tree_pit_type, 
                            type = "centroid")
beta_root_its <- betadisper(as.dist(root_aitch_its), 
                            group = pca_root_df_its$tree_pit_type, 
                            type = "centroid")

### STATISTICS BETA DISPERSION ITS BELOWGROUND PCA ###
anova(beta_osoil_its)
anova(beta_msoil_its)
anova(beta_root_its)

# Calculate and define the Tukey groups
tukey_its_osoil <- TukeyHSD(beta_osoil_its)
tukey_df_its_osoil <- as.data.frame(tukey_its_osoil$group)
x_its_osoil <- tukey_df_its_osoil$`p adj`
names(x_its_osoil) <- rownames(tukey_df_its_osoil)
tukey_groups_its_osoil <- multcompLetters(x_its_osoil)$Letters

tukey_its_msoil <- TukeyHSD(beta_msoil_its)
tukey_df_its_msoil <- as.data.frame(tukey_its_msoil$group)
x_its_msoil <- tukey_df_its_msoil$`p adj`
names(x_its_msoil) <- rownames(tukey_df_its_msoil)
tukey_groups_its_msoil <- multcompLetters(x_its_msoil)$Letters

tukey_its_root <- TukeyHSD(beta_root_its)
tukey_df_its_root <- as.data.frame(tukey_its_root$group)
x_its_root <- tukey_df_its_root$`p adj`
names(x_its_root) <- rownames(tukey_df_its_root)
tukey_groups_its_root <- multcompLetters(x_its_root)$Letters

# Make a dataframe of the beta distances
beta_osoil_dist_its <- as.data.frame(beta_osoil_its$distances)
beta_msoil_dist_its <- as.data.frame(beta_msoil_its$distances)
beta_root_dist_its <- as.data.frame(beta_root_its$distances)

# Merge the distances with the metadata
beta_osoil_df_its <- merge(beta_osoil_dist_its, metadata_its, 
                          by.x = 0, by.y = "sample_name", soil = T)
beta_msoil_df_its <- merge(beta_msoil_dist_its, metadata_its, 
                          by.x = 0, by.y = "sample_name", soil = T)
beta_root_df_its <- merge(beta_root_dist_its, metadata_its, 
                          by.x = 0, by.y = "sample_name", soil = T)

# Order the dataframe by the tree pit type
beta_osoil_df_its <- beta_osoil_df_its[order(beta_osoil_df_its$tree_pit_type),]
beta_msoil_df_its <- beta_msoil_df_its[order(beta_msoil_df_its$tree_pit_type),]
beta_root_df_its <- beta_root_df_its[order(beta_root_df_its$tree_pit_type),]

# Format the tree pit type for plotting
beta_osoil_df_its$tree_pit_type <- gsub("Interior", "Forest Interior", 
                                       beta_osoil_df_its$tree_pit_type)
beta_osoil_df_its$tree_pit_type <- gsub("Edge", "Forest Edge", 
                                       beta_osoil_df_its$tree_pit_type)
beta_osoil_df_its$tree_pit_type <- gsub("Street", "Urban Street", 
                                       beta_osoil_df_its$tree_pit_type)
beta_osoil_df_its$tree_pit_type <- gsub(" ", "\n", 
                                       beta_osoil_df_its$tree_pit_type)
beta_osoil_df_its$tree_pit_type <- factor(beta_osoil_df_its$tree_pit_type, 
                                         levels = c("Urban\nStreet\nTree", 
                                                    "Urban\nForest\nEdge",
                                                    "Urban\nForest\nInterior",
                                                    "Rural\nForest\nEdge",
                                                    "Rural\nForest\nInterior"))

beta_msoil_df_its$tree_pit_type <- gsub("Interior", "Forest Interior", 
                                       beta_msoil_df_its$tree_pit_type)
beta_msoil_df_its$tree_pit_type <- gsub("Edge", "Forest Edge", 
                                       beta_msoil_df_its$tree_pit_type)
beta_msoil_df_its$tree_pit_type <- gsub("Street", "Urban Street", 
                                       beta_msoil_df_its$tree_pit_type)
beta_msoil_df_its$tree_pit_type <- gsub(" ", "\n", 
                                       beta_msoil_df_its$tree_pit_type)
beta_msoil_df_its$tree_pit_type <- factor(beta_msoil_df_its$tree_pit_type, 
                                         levels = c("Urban\nStreet\nTree", 
                                                    "Urban\nForest\nEdge",
                                                    "Urban\nForest\nInterior",
                                                    "Rural\nForest\nEdge",
                                                    "Rural\nForest\nInterior"))

beta_root_df_its$tree_pit_type <- gsub("Interior", "Forest Interior", 
                                       beta_root_df_its$tree_pit_type)
beta_root_df_its$tree_pit_type <- gsub("Edge", "Forest Edge", 
                                       beta_root_df_its$tree_pit_type)
beta_root_df_its$tree_pit_type <- gsub("Street", "Urban Street", 
                                       beta_root_df_its$tree_pit_type)
beta_root_df_its$tree_pit_type <- gsub(" ", "\n", 
                                       beta_root_df_its$tree_pit_type)
beta_root_df_its$tree_pit_type <- factor(beta_root_df_its$tree_pit_type, 
                                         levels = c("Urban\nStreet\nTree", 
                                                    "Urban\nForest\nEdge",
                                                    "Urban\nForest\nInterior",
                                                    "Rural\nForest\nEdge",
                                                    "Rural\nForest\nInterior"))

# Define where the Tukey group labels should go
text_positions_osoil <- beta_osoil_df_its %>%
  group_by(tree_pit_type) %>%
  summarise(max_value = max(`beta_osoil_its$distances`), .groups = 'drop') %>%
  mutate(text_position = max_value + 0.5) # Add buffer for text
text_positions_osoil$label <- c("c", "a", "ab", "b", "c")

text_positions_msoil <- beta_msoil_df_its %>%
  group_by(tree_pit_type) %>%
  summarise(max_value = max(`beta_msoil_its$distances`), .groups = 'drop') %>%
  mutate(text_position = max_value + 0.5) # Add buffer for text
text_positions_msoil$label <- c("c", "a", "bc", "ab", "c")

text_positions_root <- beta_root_df_its %>%
  group_by(tree_pit_type) %>%
  summarise(max_value = max(`beta_root_its$distances`), .groups = 'drop') %>%
  mutate(text_position = max_value + 0.5) # Add buffer for text
text_positions_root$label <- c("b", "a", "ab", "ab", "ab")

figs4c_osoil <- ggplot(beta_osoil_df_its,
                 aes(x = tree_pit_type, y = `beta_osoil_its$distances`, 
                     fill = tree_pit_type)) + 
  stat_summary(aes(group = 1),
               fun = median,
               geom = "line",
               linetype = "dashed",
               linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) + 
  geom_jitter(aes(color = tree_pit_type), alpha = 0.5, 
              stroke = NA) +
  theme_classic() + xlab("") + 
  ylab("Beta dispersion") + 
  geom_text(data = text_positions_osoil, aes(x = tree_pit_type, y = text_position, 
                                       label = label,
                                       family = "Arial"), 
            inherit.aes = FALSE, vjust = -0.5) +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10),
        legend.position = "none") + 
  scale_color_manual(values = categorical_pal) +
  scale_fill_manual(values = categorical_pal) +
  ggtitle("Fungal community belowground") 

figs4c_osoil
ggsave(paste0("atherton_figures4c_ITSbetadisperosoil_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs4c_osoil, dpi = 300, width = 6, height = 3.5, units = "in")


figs4c_msoil <- ggplot(beta_msoil_df_its,
                       aes(x = tree_pit_type, y = `beta_msoil_its$distances`, 
                           fill = tree_pit_type)) + 
  stat_summary(aes(group = 1),
               fun = median,
               geom = "line",
               linetype = "dashed",
               linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) + 
  geom_jitter(aes(color = tree_pit_type), alpha = 0.5, 
              stroke = NA) +
  theme_classic() + xlab("") + 
  ylab("Beta dispersion") + 
  geom_text(data = text_positions_msoil, aes(x = tree_pit_type, y = text_position, 
                                             label = label,
                                             family = "Arial"), 
            inherit.aes = FALSE, vjust = -0.5) +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10),
        legend.position = "none") + 
  scale_fill_manual(values = categorical_pal)+
  ggtitle("Fungal community belowground") 

figs4c_msoil
ggsave(paste0("atherton_figures4c_ITSbetadispermsoil_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs4c_msoil, dpi = 300, width = 6, height = 3.5, units = "in")

figs4c_root <- ggplot(beta_root_df_its,
                       aes(x = tree_pit_type, y = `beta_root_its$distances`, 
                           fill = tree_pit_type)) + 
  stat_summary(aes(group = 1),
               fun = median,
               geom = "line",
               linetype = "dashed",
               linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) + 
  geom_jitter(aes(color = tree_pit_type), alpha = 0.5, 
              stroke = NA) +
  theme_classic() + xlab("") + 
  ylab("Beta dispersion") + 
  geom_text(data = text_positions_root, aes(x = tree_pit_type, y = text_position, 
                                             label = label,
                                             family = "Arial"), 
            inherit.aes = FALSE, vjust = -0.5) +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10),
        legend.position = "none") + 
  scale_color_manual(values = categorical_pal) +
  scale_fill_manual(values = categorical_pal) +
  ggtitle("Fungal community belowground") 

figs4c_root
ggsave(paste0("atherton_figures4c_ITSbetadisperroot_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs4c_root, dpi = 300, width = 6, height = 3.5, units = "in")

rm(list = c("beta_soil_df_its", "beta_soil_dist_its", "beta_soil_its", 
            "metadata_soil_its", "pca_soil_df_its", "pca_result_soil_its",
            "text_positions", "tukey_df_its", "tukey_its"))

### FIGURE S4D ################################################################
# Select the belowground 16s metadata
metadata_osoil_16s <- metadata_16s[which(metadata_16s$sample_type == 
                                           "Soil, 0-15 cm"),]
metadata_msoil_16s <- metadata_16s[which(metadata_16s$sample_type == 
                                           "Soil, 15-30 cm"),]
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

# Order the aitchison data by the sample ID
osoil_aitch_16s <- osoil_aitch_16s[order(row.names(osoil_aitch_16s),
                                         colnames(osoil_aitch_16s))]
msoil_aitch_16s <- msoil_aitch_16s[order(row.names(msoil_aitch_16s),
                                         colnames(msoil_aitch_16s))]
root_aitch_16s <- root_aitch_16s[order(row.names(root_aitch_16s),
                                       colnames(root_aitch_16s))]

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

# Order the PC data by the sample ID
pca_osoil_df_16s <- pca_osoil_df_16s[order(pca_osoil_df_16s$Sample),]
pca_msoil_df_16s <- pca_msoil_df_16s[order(pca_msoil_df_16s$Sample),]
pca_root_df_16s <- pca_root_df_16s[order(pca_root_df_16s$Sample),]

# Calculate beta disperson of the PCA
beta_osoil_16s <- betadisper(as.dist(osoil_aitch_16s), 
                             group = pca_osoil_df_16s$tree_pit_type, 
                             type = "centroid")
beta_msoil_16s <- betadisper(as.dist(msoil_aitch_16s), 
                             group = pca_msoil_df_16s$tree_pit_type, 
                             type = "centroid")
beta_root_16s <- betadisper(as.dist(root_aitch_16s), 
                            group = pca_root_df_16s$tree_pit_type, 
                            type = "centroid")

### STATISTICS BETA DISPERSION 16s BELOWGROUND PCA ###
anova(beta_osoil_16s)
anova(beta_msoil_16s)
anova(beta_root_16s)

# Calculate and define the Tukey groups
tukey_16s_osoil <- TukeyHSD(beta_osoil_16s)
tukey_df_16s_osoil <- as.data.frame(tukey_16s_osoil$group)
x_16s_osoil <- tukey_df_16s_osoil$`p adj`
names(x_16s_osoil) <- rownames(tukey_df_16s_osoil)
tukey_groups_16s_osoil <- multcompLetters(x_16s_osoil)$Letters

tukey_16s_msoil <- TukeyHSD(beta_msoil_16s)
tukey_df_16s_msoil <- as.data.frame(tukey_16s_msoil$group)
x_16s_msoil <- tukey_df_16s_msoil$`p adj`
names(x_16s_msoil) <- rownames(tukey_df_16s_msoil)
tukey_groups_16s_msoil <- multcompLetters(x_16s_msoil)$Letters

tukey_16s_root <- TukeyHSD(beta_root_16s)
tukey_df_16s_root <- as.data.frame(tukey_16s_root$group)
x_16s_root <- tukey_df_16s_root$`p adj`
names(x_16s_root) <- rownames(tukey_df_16s_root)
tukey_groups_16s_root <- multcompLetters(x_16s_root)$Letters

# Make a dataframe of the beta distances
beta_osoil_dist_16s <- as.data.frame(beta_osoil_16s$distances)
beta_msoil_dist_16s <- as.data.frame(beta_msoil_16s$distances)
beta_root_dist_16s <- as.data.frame(beta_root_16s$distances)

# Merge the distances with the metadata
beta_osoil_df_16s <- merge(beta_osoil_dist_16s, metadata_16s, 
                           by.x = 0, by.y = "sample_name", soil = T)
beta_msoil_df_16s <- merge(beta_msoil_dist_16s, metadata_16s, 
                           by.x = 0, by.y = "sample_name", soil = T)
beta_root_df_16s <- merge(beta_root_dist_16s, metadata_16s, 
                          by.x = 0, by.y = "sample_name", soil = T)

# Order the dataframe by the tree pit type
beta_osoil_df_16s <- beta_osoil_df_16s[order(beta_osoil_df_16s$tree_pit_type),]
beta_msoil_df_16s <- beta_msoil_df_16s[order(beta_msoil_df_16s$tree_pit_type),]
beta_root_df_16s <- beta_root_df_16s[order(beta_root_df_16s$tree_pit_type),]

# Format the tree pit type for plotting
beta_osoil_df_16s$tree_pit_type <- gsub("Interior", "Forest Interior", 
                                        beta_osoil_df_16s$tree_pit_type)
beta_osoil_df_16s$tree_pit_type <- gsub("Edge", "Forest Edge", 
                                        beta_osoil_df_16s$tree_pit_type)
beta_osoil_df_16s$tree_pit_type <- gsub("Street", "Urban Street", 
                                        beta_osoil_df_16s$tree_pit_type)
beta_osoil_df_16s$tree_pit_type <- gsub(" ", "\n", 
                                        beta_osoil_df_16s$tree_pit_type)
beta_osoil_df_16s$tree_pit_type <- factor(beta_osoil_df_16s$tree_pit_type, 
                                          levels = c("Urban\nStreet\nTree", 
                                                     "Urban\nForest\nEdge",
                                                     "Urban\nForest\nInterior",
                                                     "Rural\nForest\nEdge",
                                                     "Rural\nForest\nInterior"))

beta_msoil_df_16s$tree_pit_type <- gsub("Interior", "Forest Interior", 
                                        beta_msoil_df_16s$tree_pit_type)
beta_msoil_df_16s$tree_pit_type <- gsub("Edge", "Forest Edge", 
                                        beta_msoil_df_16s$tree_pit_type)
beta_msoil_df_16s$tree_pit_type <- gsub("Street", "Urban Street", 
                                        beta_msoil_df_16s$tree_pit_type)
beta_msoil_df_16s$tree_pit_type <- gsub(" ", "\n", 
                                        beta_msoil_df_16s$tree_pit_type)
beta_msoil_df_16s$tree_pit_type <- factor(beta_msoil_df_16s$tree_pit_type, 
                                          levels = c("Urban\nStreet\nTree", 
                                                     "Urban\nForest\nEdge",
                                                     "Urban\nForest\nInterior",
                                                     "Rural\nForest\nEdge",
                                                     "Rural\nForest\nInterior"))

beta_root_df_16s$tree_pit_type <- gsub("Interior", "Forest Interior", 
                                       beta_root_df_16s$tree_pit_type)
beta_root_df_16s$tree_pit_type <- gsub("Edge", "Forest Edge", 
                                       beta_root_df_16s$tree_pit_type)
beta_root_df_16s$tree_pit_type <- gsub("Street", "Urban Street", 
                                       beta_root_df_16s$tree_pit_type)
beta_root_df_16s$tree_pit_type <- gsub(" ", "\n", 
                                       beta_root_df_16s$tree_pit_type)
beta_root_df_16s$tree_pit_type <- factor(beta_root_df_16s$tree_pit_type, 
                                         levels = c("Urban\nStreet\nTree", 
                                                    "Urban\nForest\nEdge",
                                                    "Urban\nForest\nInterior",
                                                    "Rural\nForest\nEdge",
                                                    "Rural\nForest\nInterior"))

# Define where the Tukey group labels should go
text_positions_osoil <- beta_osoil_df_16s %>%
  group_by(tree_pit_type) %>%
  summarise(max_value = max(`beta_osoil_16s$distances`), .groups = 'drop') %>%
  mutate(text_position = max_value + 0.5) # Add buffer for text
text_positions_osoil$label <- c("b", "a", "ab", "ab", "b")

text_positions_msoil <- beta_msoil_df_16s %>%
  group_by(tree_pit_type) %>%
  summarise(max_value = max(`beta_msoil_16s$distances`), .groups = 'drop') %>%
  mutate(text_position = max_value + 0.5) # Add buffer for text
text_positions_msoil$label <- c("bc", "a", "abc", "ab", "c")

text_positions_root <- beta_root_df_16s %>%
  group_by(tree_pit_type) %>%
  summarise(max_value = max(`beta_root_16s$distances`), .groups = 'drop') %>%
  mutate(text_position = max_value + 0.5) # Add buffer for text
text_positions_root$label <- c("b", "a", "a", "a", "a")

figs4c_osoil <- ggplot(beta_osoil_df_16s,
                       aes(x = tree_pit_type, y = `beta_osoil_16s$distances`, 
                           fill = tree_pit_type)) + 
  stat_summary(aes(group = 1),
               fun = median,
               geom = "line",
               linetype = "dashed",
               linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) + 
  geom_jitter(aes(color = tree_pit_type), alpha = 0.5, 
              stroke = NA) +
  theme_classic() + xlab("") + 
  ylab("Beta dispersion") + 
  geom_text(data = text_positions_osoil, aes(x = tree_pit_type, y = text_position, 
                                             label = label,
                                             family = "Arial"), 
            inherit.aes = FALSE, vjust = -0.5) +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10),
        legend.position = "none") + 
  scale_color_manual(values = categorical_pal) +
  scale_fill_manual(values = categorical_pal) +
  ggtitle("Fungal community belowground") 

figs4c_osoil
ggsave(paste0("atherton_figures4c_16sbetadisperosoil_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs4c_osoil, dpi = 300, width = 6, height = 3.5, units = "in")


figs4c_msoil <- ggplot(beta_msoil_df_16s,
                       aes(x = tree_pit_type, y = `beta_msoil_16s$distances`, 
                           fill = tree_pit_type)) + 
  stat_summary(aes(group = 1),
               fun = median,
               geom = "line",
               linetype = "dashed",
               linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) + 
  geom_jitter(aes(color = tree_pit_type), alpha = 0.5, 
              stroke = NA) +
  theme_classic() + xlab("") + 
  ylab("Beta dispersion") + 
  geom_text(data = text_positions_msoil, aes(x = tree_pit_type, y = text_position, 
                                             label = label,
                                             family = "Arial"), 
            inherit.aes = FALSE, vjust = -0.5) +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10),
        legend.position = "none") + 
  scale_color_manual(values = categorical_pal) +
  scale_fill_manual(values = categorical_pal) +
  ggtitle("Fungal community belowground") 

figs4c_msoil
ggsave(paste0("atherton_figures4c_16sbetadispermsoil_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs4c_msoil, dpi = 300, width = 6, height = 3.5, units = "in")

figs4c_root <- ggplot(beta_root_df_16s,
                      aes(x = tree_pit_type, y = `beta_root_16s$distances`, 
                          fill = tree_pit_type)) + 
  stat_summary(aes(group = 1),
               fun = median,
               geom = "line",
               linetype = "dashed",
               linewidth = 0.8) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) + 
  geom_jitter(aes(color = tree_pit_type), alpha = 0.5, 
              stroke = NA) +
  theme_classic() + xlab("") + 
  ylab("Beta dispersion") + 
  geom_text(data = text_positions_root, aes(x = tree_pit_type, y = text_position, 
                                            label = label,
                                            family = "Arial"), 
            inherit.aes = FALSE, vjust = -0.5) +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10),
        legend.position = "none") + 
  scale_color_manual(values = categorical_pal) +
  scale_fill_manual(values = categorical_pal) +
  ggtitle("Fungal community belowground") 

figs4c_root
ggsave(paste0("atherton_figures4c_16sbetadisperroot_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs4c_root, dpi = 300, width = 6, height = 3.5, units = "in")