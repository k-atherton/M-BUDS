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

pca_leaf_df_its$tree_age <- "O"
pca_leaf_df_its$tree_age[which(pca_leaf_df_its$tree_dbh_2021 < 30)] <- "M"
pca_leaf_df_its$tree_age[which(pca_leaf_df_its$tree_dbh_2021 < 10)] <- "Y"

print(adonis2(leaf_aitch_its~pca_leaf_df_its$tree_pit_type))
print(adonis2(leaf_aitch_its~pca_leaf_df_its$tree_age))
print(adonis2(leaf_aitch_its~pca_leaf_df_its$tree_species))

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
        plot.title = element_text(size = 10), legend.position = "none") +
  ggtitle("Fungal community in leaves") + 
  scale_color_manual(values = categorical_pal)

figs3a
ggsave(paste0("atherton_figures3a_ITSpcaleaves_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs3a, dpi = 300, width = 6, height = 3.5, units = "in")
