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

### FIGURE S4C #################################################################
# Select the belowground ITS metadata
metadata_osoil_its <- metadata_its[which(metadata_its$sample_type == 
                                          "Soil, 0-15 cm"),]

# Select the belowground samples from the aitchison data
osoil_aitch_its <- all_aitch_its[which(row.names(all_aitch_its) %in%
                                        metadata_osoil_its$sample_name),
                                which(colnames(all_aitch_its) %in%
                                        metadata_osoil_its$sample_name)]

# Order the aitchison data by the sample ID
osoil_aitch_its <- osoil_aitch_its[order(row.names(osoil_aitch_its),
                                       colnames(osoil_aitch_its))]

# Run PCA
pca_result_osoil_its <- prcomp(osoil_aitch_its, scale. = TRUE)

# Make a dataframe of the PCs
pca_osoil_df_its <- data.frame(Sample = rownames(osoil_aitch_its), 
                              PC1 = pca_result_osoil_its$x[,1], 
                              PC2 = pca_result_osoil_its$x[,2])

# Merge the PCs with the belowground metadata
pca_osoil_df_its <- merge(pca_osoil_df_its, metadata_osoil_its, 
                         by.x = "Sample", by.y = "sample_name", all.x = T)

# Order the PC data by the sample ID
pca_osoil_df_its <- pca_osoil_df_its[order(pca_osoil_df_its$Sample),]

# Calculate beta disperson of the PCA
beta_osoil_its <- betadisper(as.dist(osoil_aitch_its), 
                            group = pca_osoil_df_its$tree_pit_type, 
                            type = "centroid")

### STATISTICS BETA DISPERSION ITS BELOWGROUND PCA ###
anova(beta_osoil_its)

# Calculate and define the Tukey groups
tukey_its_osoil <- TukeyHSD(beta_osoil_its)
tukey_df_its_osoil <- as.data.frame(tukey_its_osoil$group)
x_its_osoil <- tukey_df_its_osoil$`p adj`
names(x_its_osoil) <- rownames(tukey_df_its_osoil)
tukey_groups_its_osoil <- multcompLetters(x_its_osoil)$Letters

# Make a dataframe of the beta distances
beta_osoil_dist_its <- as.data.frame(beta_osoil_its$distances)

# Merge the distances with the metadata
beta_osoil_df_its <- merge(beta_osoil_dist_its, metadata_its, 
                          by.x = 0, by.y = "sample_name", soil = T)

# Order the dataframe by the tree pit type
beta_osoil_df_its <- beta_osoil_df_its[order(beta_osoil_df_its$tree_pit_type),]

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

# Define where the Tukey group labels should go
text_positions_osoil <- beta_osoil_df_its %>%
  group_by(tree_pit_type) %>%
  summarise(max_value = max(`beta_osoil_its$distances`), .groups = 'drop') %>%
  mutate(text_position = max_value + 0.1) # Add buffer for text
text_positions_osoil$label <- c("d", "a", "ab", "bc", "cd")

figs4c_osoil <- ggplot(beta_osoil_df_its,
                 aes(x = tree_pit_type, y = `beta_osoil_its$distances`, 
                     fill = tree_pit_type)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(aes(color = tree_pit_type), alpha = 0.75, 
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
  ggtitle("Fungal community in soil, 0-15 cm") + ylim(c(20,70))

figs4c_osoil
ggsave(paste0("atherton_figures4c_ITSbetadisperosoil_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs4c_osoil, dpi = 300, width = 6, height = 3.5, units = "in")
