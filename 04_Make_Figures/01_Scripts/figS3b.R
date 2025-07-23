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
metadata_its <- read.csv("atherton_sample_metadata_ITS_20250125.csv")
metadata_its <- distinct(metadata_its)

# read in sample names of cleaned data
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/02_Clean_Data/05_Transform_Data/ITS/Aitchison_Distance_Tables/")
all_aitch_its <- read.csv("atherton_ITS_allsampletypes_aitchisondistance_20240919.csv", 
                          row.names = 1)

### FORMAT METADATA ###########################################################
# Combine ITS soil horizon with sample type
metadata_its$sample_type<- paste0(metadata_its$soil_horizon, " ", 
                                  metadata_its$sample_type)
metadata_its$sample_type <- gsub("NA ", "", 
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

figS3b <- ggplot(beta_leaf_df_its, 
                 aes(x = tree_pit_type, y = `beta_leaf_its$distances`, 
                     fill = tree_pit_type)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(aes(color = tree_pit_type), alpha = 0.75, 
              stroke = NA) + theme_classic() + xlab("") + 
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
figS3b

ggsave(paste0("atherton_figureS3b_ITSbetadisperleaves_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figS3b, dpi = 300, width = 6, height = 3.5, units = "in")
