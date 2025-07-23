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

### FIGURE S4D ################################################################
# Select the belowground 16s metadata
metadata_msoil_16s <- metadata_16s[which(metadata_16s$sample_type == 
                                           "Soil, 15-30 cm"),]

# Select the belowground samples from the aitchison data
msoil_aitch_16s <- bg_aitch_16s[which(row.names(bg_aitch_16s) %in%
                                         metadata_msoil_16s$sample_name),
                                 which(colnames(bg_aitch_16s) %in%
                                         metadata_msoil_16s$sample_name)]

# Order the aitchison data by the sample ID
msoil_aitch_16s <- msoil_aitch_16s[order(row.names(msoil_aitch_16s),
                                         colnames(msoil_aitch_16s))]

# Run PCA
pca_result_msoil_16s <- prcomp(msoil_aitch_16s, scale. = TRUE)

# Make a dataframe of the PCs
pca_msoil_df_16s <- data.frame(Sample = rownames(msoil_aitch_16s), 
                               PC1 = pca_result_msoil_16s$x[,1], 
                               PC2 = pca_result_msoil_16s$x[,2])

# Merge the PCs with the belowground metadata
pca_msoil_df_16s <- merge(pca_msoil_df_16s, metadata_msoil_16s, 
                          by.x = "Sample", by.y = "sample_name", all.x = T)

# Order the PC data by the sample ID
pca_msoil_df_16s <- pca_msoil_df_16s[order(pca_msoil_df_16s$Sample),]

# Calculate beta disperson of the PCA
beta_msoil_16s <- betadisper(as.dist(msoil_aitch_16s), 
                             group = pca_msoil_df_16s$tree_pit_type, 
                             type = "centroid")

### STATISTICS BETA DISPERSION 16s BELOWGROUND PCA ###
anova(beta_msoil_16s)

# Calculate and define the Tukey groups
tukey_16s_msoil <- TukeyHSD(beta_msoil_16s)
tukey_df_16s_msoil <- as.data.frame(tukey_16s_msoil$group)
x_16s_msoil <- tukey_df_16s_msoil$`p adj`
names(x_16s_msoil) <- rownames(tukey_df_16s_msoil)
tukey_groups_16s_msoil <- multcompLetters(x_16s_msoil)$Letters

# Make a dataframe of the beta distances
beta_msoil_dist_16s <- as.data.frame(beta_msoil_16s$distances)

# Merge the distances with the metadata
beta_msoil_df_16s <- merge(beta_msoil_dist_16s, metadata_16s, 
                           by.x = 0, by.y = "sample_name", soil = T)

# Order the dataframe by the tree pit type
beta_msoil_df_16s <- beta_msoil_df_16s[order(beta_msoil_df_16s$tree_pit_type),]

# Format the tree pit type for plotting
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

# Define where the Tukey group labels should go
text_positions_msoil <- beta_msoil_df_16s %>%
  group_by(tree_pit_type) %>%
  summarise(max_value = max(`beta_msoil_16s$distances`), .groups = 'drop') %>%
  mutate(text_position = max_value + 0.5) # Add buffer for text
text_positions_msoil$label <- c("b", "a", "ab", "ab", "b")

figs4h_msoil <- ggplot(beta_msoil_df_16s,
                       aes(x = tree_pit_type, y = `beta_msoil_16s$distances`, 
                           fill = tree_pit_type)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(aes(color = tree_pit_type), alpha = 0.75, 
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
  ggtitle("Bacterial community in soil, 15-30 cm") + ylim(c(30,150))

figs4h_msoil
ggsave(paste0("atherton_figures4h_16sbetadispermsoil_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs4h_msoil, dpi = 300, width = 6, height = 3.5, units = "in")
