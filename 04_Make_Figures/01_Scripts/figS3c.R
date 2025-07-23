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

### FIGURE S3B #################################################################
# Select the leaf ITS metadata
metadata_leaf_16s <- metadata_16s[which(metadata_16s$sample_type == 
                                          "Leaf"),]

# Select the leaf samples from the aitchison data
leaf_aitch_16s <- leaf_aitch_16s[which(row.names(leaf_aitch_16s) %in%
                                         metadata_leaf_16s$sample_name),
                                 which(colnames(leaf_aitch_16s) %in%
                                         metadata_leaf_16s$sample_name)]

# Run PCA
pca_result_leaf_16s <- prcomp(leaf_aitch_16s, scale. = TRUE)

# Make a dataframe of the PCs
pca_leaf_df_16s <- data.frame(Sample = rownames(leaf_aitch_16s), 
                              PC1 = pca_result_leaf_16s$x[,1], 
                              PC2 = pca_result_leaf_16s$x[,2])

# Merge the PCs with the leaf metadata
pca_leaf_df_16s <- merge(pca_leaf_df_16s, metadata_leaf_16s, 
                         by.x = "Sample", by.y = "sample_name", all = T)

# Summarize the leaf PC data
pca_summary_leaf_16s <- summary(pca_result_leaf_16s)
pca_variance_leaf_16s <- pca_summary_leaf_16s$importance[2, ]  # Proportion of variance
pca_variance_cumulative_leaf_16s <- pca_summary_leaf_16s$importance[3, ]  # Cumulative proportion of variance

### STATISTICS LEAF PCA ###
# define the spatial matrix
dec.degrees.16s<-as.matrix(cbind(pca_leaf_df_16s$tree_longitude, 
                                 pca_leaf_df_16s$tree_latitude))
rownames(dec.degrees.16s) <-pca_leaf_df_16s$Sample 

# find the spatial distance between samples
space.matrix_16s<-rdist.earth(dec.degrees.16s[,1:2],miles=FALSE)

# define the design matrix
design_16s <- as.data.frame(matrix(nrow = nrow(pca_leaf_df_16s), ncol = 5, data = 0))
rownames(design_16s) <- pca_leaf_df_16s$Sample
colnames(design_16s) <- c("street_tree", "urban_edge", "urban_interior", 
                          "rural_edge", "rural_interior")
for(i in 1:nrow(design_16s)){
  tree_type <- pca_leaf_df_16s$tree_pit_type[which(pca_leaf_df_16s$Sample == rownames(design_16s)[i])]
  if(tree_type == "Street Tree"){
    design_16s$street_tree[i] <- 1
  } else if(tree_type == "Urban Edge"){
    design_16s$urban_edge[i] <- 1
  } else if(tree_type == "Urban Interior"){
    design_16s$urban_interior[i] <- 1
  } else if(tree_type == "Rural Edge"){
    design_16s$rural_edge[i] <- 1
  } else if(tree_type == "Rural Interior"){
    design_16s$rural_interior[i] <- 1
  }
}

# Find the design distance between samples
design.dist_16s <- dist(design_16s)

# Test the effect of space and design on the PCA groupings
MRM <- lm(as.numeric(as.dist(leaf_aitch_16s)) ~  as.numeric(design.dist_16s) + as.numeric(as.dist(space.matrix_16s)))
anova(MRM)

pca_leaf_df_16s$tree_age <- "O"
pca_leaf_df_16s$tree_age[which(pca_leaf_df_16s$tree_dbh_2021 < 30)] <- "M"
pca_leaf_df_16s$tree_age[which(pca_leaf_df_16s$tree_dbh_2021 < 10)] <- "Y"

print(adonis2(leaf_aitch_16s~pca_leaf_df_16s$tree_pit_type))
print(adonis2(leaf_aitch_16s~pca_leaf_df_16s$tree_age))
print(adonis2(leaf_aitch_16s~pca_leaf_df_16s$tree_species))

# Plot figure s3b
figs3b <- ggplot(pca_leaf_df_16s, 
                 aes(x = PC1, y = PC2, color = tree_pit_type)) +
  geom_point(shape = 18, alpha = 0.5, size = 2.5) +
  xlab(paste("Principal Component 1\n(", 
             round(pca_variance_leaf_16s[1] * 100, 2), "% variance)", sep = "")) +
  ylab(paste("Principal Component 2\n(", 
             round(pca_variance_leaf_16s[2] * 100, 2), "% variance)", sep = "")) +
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
  ggtitle("Bacterial community in leaves") + 
  scale_color_manual(values = categorical_pal)

figs3b
ggsave(paste0("atherton_figures3c_16Spcaleaves_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs3b, dpi = 300, width = 6, height = 3.5, units = "in")