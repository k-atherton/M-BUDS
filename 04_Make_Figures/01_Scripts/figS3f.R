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

# read in sample names of cleaned data
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/02_Clean_Data/05_Transform_Data/16S/Aitchison_Distance_Tables/")
bg_aitch_16s <- read.csv("atherton_16S_allsampletypes_aitchisondistance_20240925.csv", 
                         row.names = 1)
leaf_aitch_16s <- read.csv("atherton_16S_leaf_aitchisondistance_20241112.csv", 
                           row.names = 1)

### FORMAT METADATA ##########################################################
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

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE S4C #################################################################
# Select the belowground ITS metadata
metadata_soil_its <- metadata_its[which(metadata_its$sample_type != 
                                          "Leaf"),]

# Select the belowground samples from the aitchison data
soil_aitch_its <- all_aitch_its[which(row.names(all_aitch_its) %in%
                                        metadata_soil_its$sample_name),
                                which(colnames(all_aitch_its) %in%
                                        metadata_soil_its$sample_name)]

# Order the aitchison data by the sample ID
soil_aitch_its <- soil_aitch_its[order(row.names(soil_aitch_its),
                                       colnames(soil_aitch_its))]

# Run PCA
pca_result_soil_its <- prcomp(soil_aitch_its, scale. = TRUE)

# Make a dataframe of the PCs
pca_soil_df_its <- data.frame(Sample = rownames(soil_aitch_its), 
                              PC1 = pca_result_soil_its$x[,1], 
                              PC2 = pca_result_soil_its$x[,2])

# Merge the PCs with the belowground metadata
pca_soil_df_its <- merge(pca_soil_df_its, metadata_soil_its, 
                         by.x = "Sample", by.y = "sample_name", all.x = T)

# Order the PC data by the sample ID
pca_soil_df_its <- pca_soil_df_its[order(pca_soil_df_its$Sample),]

# Calculate beta disperson of the PCA
beta_soil_its <- betadisper(as.dist(soil_aitch_its), 
                            group = pca_soil_df_its$tree_pit_type, 
                            type = "centroid")

### STATISTICS BETA DISPERSION ITS BELOWGROUND PCA ###
anova(beta_soil_its)

# Calculate and define the Tukey groups
tukey_its <- TukeyHSD(beta_soil_its)
tukey_df_its <- as.data.frame(tukey_its$group)
x_its <- tukey_df_its$`p adj`
names(x_its) <- rownames(tukey_df_its)
tukey_groups_its <- multcompLetters(x_its)$Letters

# Make a dataframe of the beta distances
beta_soil_dist_its <- as.data.frame(beta_soil_its$distances)

# Merge the distances with the metadata
beta_soil_df_its <- merge(beta_soil_dist_its, metadata_its, 
                          by.x = 0, by.y = "sample_name", soil = T)

# Order the dataframe by the tree pit type
beta_soil_df_its <- beta_soil_df_its[order(beta_soil_df_its$tree_pit_type),]

# Format the tree pit type for plotting
beta_soil_df_its$tree_pit_type <- gsub("Interior", "Forest Interior", 
                                       beta_soil_df_its$tree_pit_type)
beta_soil_df_its$tree_pit_type <- gsub("Edge", "Forest Edge", 
                                       beta_soil_df_its$tree_pit_type)
beta_soil_df_its$tree_pit_type <- gsub("Street", "Urban Street", 
                                       beta_soil_df_its$tree_pit_type)
beta_soil_df_its$tree_pit_type <- gsub(" ", "\n", 
                                       beta_soil_df_its$tree_pit_type)
beta_soil_df_its$tree_pit_type <- factor(beta_soil_df_its$tree_pit_type, 
                                         levels = c("Urban\nStreet\nTree", 
                                                    "Urban\nForest\nEdge",
                                                    "Urban\nForest\nInterior",
                                                    "Rural\nForest\nEdge",
                                                    "Rural\nForest\nInterior"))

# Define where the Tukey group labels should go
text_positions <- beta_soil_df_its %>%
  group_by(tree_pit_type) %>%
  summarise(max_value = max(`beta_soil_its$distances`), .groups = 'drop') %>%
  mutate(text_position = max_value + 0.5) # Add buffer for text
text_positions$label <- c("c", "a", "ab", "b", "c")

figs3f <- ggplot(beta_soil_df_its,
                 aes(x = tree_pit_type, y = `beta_soil_its$distances`, 
                     fill = tree_pit_type)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(aes(color = tree_pit_type, shape = sample_type), alpha = 0.75) +
  theme_classic() + xlab("") + 
  ylab("Beta dispersion") + 
  geom_text(data = text_positions, aes(x = tree_pit_type, y = text_position, 
                                       label = label,
                                       family = "Arial"), 
            inherit.aes = FALSE, vjust = -0.5) +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 0), 
        plot.title = element_text(size = 10),
        legend.position = "bottom") + 
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = categorical_pal) +
  scale_fill_manual(values = categorical_pal) +
  guides(color = "none", fill = "none") +
  ggtitle("Fungal community belowground") 

figs3f
ggsave(paste0("atherton_figures3f_ITSbetadisperbelowground_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs3f, dpi = 300, width = 6, height = 3.5, units = "in")
