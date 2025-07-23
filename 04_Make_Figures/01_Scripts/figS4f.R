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

### FIGURE s4f ################################################################
# Select the belowground 16s metadata
metadata_msoil_16s <- metadata_16s[which(metadata_16s$sample_type == 
                                           "Soil, 15-30 cm"),]

length(unique(metadata_msoil_16s$tree_id))

# Select the belowground samples from the aitchison data
msoil_aitch_16s <- bg_aitch_16s[which(row.names(bg_aitch_16s) %in%
                                         metadata_msoil_16s$sample_name),
                                 which(colnames(bg_aitch_16s) %in%
                                         metadata_msoil_16s$sample_name)]
# Run PCA
pca_result_msoil_16s <- prcomp(msoil_aitch_16s, scale. = TRUE)

# Make a dataframe of the PCs
pca_msoil_df_16s <- data.frame(Sample = rownames(msoil_aitch_16s), 
                               PC1 = pca_result_msoil_16s$x[,1], 
                               PC2 = pca_result_msoil_16s$x[,2])

# Merge the PCs with the belowground metadata
pca_msoil_df_16s <- merge(pca_msoil_df_16s, metadata_msoil_16s, 
                          by.x = "Sample", by.y = "sample_name", all.x = T)

# Summarize the belowground PC data
pca_summary_msoil_16s <- summary(pca_result_msoil_16s)
pca_variance_msoil_16s <- pca_summary_msoil_16s$importance[2, ]  # Proportion of variance
pca_variance_cumulative_msoil_16s <- pca_summary_msoil_16s$importance[3, ]  # Cumulative proportion of variance

### STATISTICS BELOWGROUND PCA ###
# define the spatial matrix
dec.degrees.16s_msoil <- as.matrix(cbind(pca_msoil_df_16s$tree_longitude, 
                                         pca_msoil_df_16s$tree_latitude))
rownames(dec.degrees.16s_msoil) <- pca_msoil_df_16s$Sample 

# find the spatial distance between samples
space.matrix_16s_msoil <-rdist.earth(dec.degrees.16s_msoil[,1:2], miles=FALSE)

# define the design matrix
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

# Find the design distance between samples
design.dist_16s_msoil <- dist(design_16s_msoil)

# Test the effect of space and design on the PCA groupings
MRM_msoil <- lm(as.numeric(as.dist(msoil_aitch_16s)) ~  as.numeric(design.dist_16s_msoil) + as.numeric(as.dist(space.matrix_16s_msoil)))
anova(MRM_msoil)

pca_msoil_df_16s$tree_age <- "O"
pca_msoil_df_16s$tree_age[which(pca_msoil_df_16s$tree_dbh_2021 < 30)] <- "M"
pca_msoil_df_16s$tree_age[which(pca_msoil_df_16s$tree_dbh_2021 < 10)] <- "Y"

print(adonis2(msoil_aitch_16s~pca_msoil_df_16s$tree_pit_type))
print(adonis2(msoil_aitch_16s~pca_msoil_df_16s$tree_age))
print(adonis2(msoil_aitch_16s~pca_msoil_df_16s$tree_species))

# Plot figure s4f
figs4f_msoil <- ggplot(pca_msoil_df_16s, 
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

figs4f_msoil
ggsave(paste0("atherton_figures4f_16spcamsoil_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs4f_msoil, dpi = 300, width = 6, height = 3.5, units = "in")