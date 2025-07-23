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

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE s4g #################################################################
# Select the belowground ITS metadata
metadata_msoil_its <- metadata_its[which(metadata_its$sample_type == 
                                           "Soil, 15-30 cm"),]

# Select the belowground samples from the aitchison data
msoil_aitch_its <- all_aitch_its[which(row.names(all_aitch_its) %in%
                                         metadata_msoil_its$sample_name),
                                 which(colnames(all_aitch_its) %in%
                                         metadata_msoil_its$sample_name)]

# Order the aitchison data by the sample ID
msoil_aitch_its <- msoil_aitch_its[order(row.names(msoil_aitch_its),
                                         colnames(msoil_aitch_its))]

# Run PCA
pca_result_msoil_its <- prcomp(msoil_aitch_its, scale. = TRUE)

# Make a dataframe of the PCs
pca_msoil_df_its <- data.frame(Sample = rownames(msoil_aitch_its), 
                              PC1 = pca_result_msoil_its$x[,1], 
                              PC2 = pca_result_msoil_its$x[,2])

# Merge the PCs with the belowground metadata
pca_msoil_df_its <- merge(pca_msoil_df_its, metadata_msoil_its, 
                          by.x = "Sample", by.y = "sample_name", all.x = T)

# Order the PC data by the sample ID
pca_msoil_df_its <- pca_msoil_df_its[order(pca_msoil_df_its$Sample),]

# Calculate beta disperson of the PCA
beta_msoil_its <- betadisper(as.dist(msoil_aitch_its), 
                            group = pca_msoil_df_its$tree_pit_type, 
                            type = "centroid")

### STATISTICS BETA DISPERSION ITS BELOWGROUND PCA ###
anova(beta_msoil_its)

# Calculate and define the Tukey groups
tukey_its_msoil <- TukeyHSD(beta_msoil_its)
tukey_df_its_msoil <- as.data.frame(tukey_its_msoil$group)
x_its_msoil <- tukey_df_its_msoil$`p adj`
names(x_its_msoil) <- rownames(tukey_df_its_msoil)
tukey_groups_its_msoil <- multcompLetters(x_its_msoil)$Letters

# Make a dataframe of the beta distances
beta_msoil_dist_its <- as.data.frame(beta_msoil_its$distances)

# Merge the distances with the metadata
beta_msoil_df_its <- merge(beta_msoil_dist_its, metadata_its, 
                          by.x = 0, by.y = "sample_name", soil = T)

# Order the dataframe by the tree pit type
beta_msoil_df_its <- beta_msoil_df_its[order(beta_msoil_df_its$tree_pit_type),]

# Format the tree pit type for plotting
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

# Define where the Tukey group labels should go
text_positions_msoil <- beta_msoil_df_its %>%
  group_by(tree_pit_type) %>%
  summarise(max_value = max(`beta_msoil_its$distances`), .groups = 'drop') %>%
  mutate(text_position = max_value + 0.5) # Add buffer for text
text_positions_msoil$label <- c("c", "a", "bc", "ab", "c")


# Plot Fig S4g
figs4g_msoil <- ggplot(beta_msoil_df_its,
                       aes(x = tree_pit_type, y = `beta_msoil_its$distances`, 
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
  scale_fill_manual(values = categorical_pal)+
  scale_color_manual(values = categorical_pal)+
  ggtitle("Fungal community in soil, 15-30 cm") + ylim(c(20,75)) 

figs4g_msoil
ggsave(paste0("atherton_figures4g_ITSbetadispermsoil_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs4g_msoil, dpi = 300, width = 6, height = 3.5, units = "in")
