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

### FORMAT METADATA ###########################################################
# Combine 16S soil horizon with sample type
metadata_16s$sample_type<- paste0(metadata_16s$soil_horizon, " ",
                                  metadata_16s$sample_type)
metadata_16s$sample_type <- gsub("NA ", "", 
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

figs3d <- ggplot(beta_leaf_df_16s,
                 aes(x = tree_pit_type, y = `beta_leaf_16s$distances`, 
                     fill = tree_pit_type)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(aes(color = tree_pit_type), alpha = 0.75, 
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

figs3d
ggsave(paste0("atherton_figures3d_16Sbetadisperleaves_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs3d, dpi = 300, width = 6, height = 3.5, units = "in")
