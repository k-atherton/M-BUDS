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

### FIGURE s4e #################################################################
# Select the belowground ITS metadata
metadata_msoil_its <- metadata_its[which(metadata_its$sample_type == 
                                           "Soil, 15-30 cm"),]

length(unique(metadata_msoil_its$tree_id))

# Select the belowground samples from the aitchison data
msoil_aitch_its <- all_aitch_its[which(row.names(all_aitch_its) %in%
                                         metadata_msoil_its$sample_name),
                                 which(colnames(all_aitch_its) %in%
                                         metadata_msoil_its$sample_name)]

# Run PCA
pca_result_msoil_its <- prcomp(msoil_aitch_its, scale. = TRUE)

# Make a dataframe of the PCs
pca_msoil_df_its <- data.frame(Sample = rownames(msoil_aitch_its), 
                               PC1 = pca_result_msoil_its$x[,1], 
                               PC2 = pca_result_msoil_its$x[,2])

# Merge the PCs with the belowground metadata
pca_msoil_df_its <- merge(pca_msoil_df_its, metadata_msoil_its, 
                          by.x = "Sample", by.y = "sample_name", all.x = T)

# Summarize the belowground PC data
pca_summary_msoil_its <- summary(pca_result_msoil_its)
pca_variance_msoil_its <- pca_summary_msoil_its$importance[2, ]  # Proportion of variance
pca_variance_cumulative_msoil_its <- pca_summary_msoil_its$importance[3, ]  # Cumulative proportion of variance

### STATISTICS BELOWGROUND PCA ###
# define the spatial matrix
dec.degrees.its_msoil <- as.matrix(cbind(pca_msoil_df_its$tree_longitude, 
                                         pca_msoil_df_its$tree_latitude))
rownames(dec.degrees.its_msoil) <- pca_msoil_df_its$Sample 

# find the spatial distance between samples
space.matrix_its_msoil <-rdist.earth(dec.degrees.its_msoil[,1:2], miles=FALSE)

# define the design matrix
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

# Find the design distance between samples
design.dist_its_msoil <- dist(design_its_msoil)

# Test the effect of space and design on the PCA groupings
MRM_msoil <- lm(as.numeric(as.dist(msoil_aitch_its)) ~  as.numeric(design.dist_its_msoil) + as.numeric(as.dist(space.matrix_its_msoil)))
anova(MRM_msoil)

pca_msoil_df_its$tree_age <- "O"
pca_msoil_df_its$tree_age[which(pca_msoil_df_its$tree_dbh_2021 < 30)] <- "M"
pca_msoil_df_its$tree_age[which(pca_msoil_df_its$tree_dbh_2021 < 10)] <- "Y"

print(adonis2(msoil_aitch_its~pca_msoil_df_its$tree_pit_type))
print(adonis2(msoil_aitch_its~pca_msoil_df_its$tree_age))
print(adonis2(msoil_aitch_its~pca_msoil_df_its$tree_species))

# Plot figure s4e
figs4e_msoil <- ggplot(pca_msoil_df_its, 
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

figs4e_msoil
ggsave(paste0("atherton_figures4e_ITSpcamsoil_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs4e_msoil, dpi = 300, width = 6, height = 3.5, units = "in")

