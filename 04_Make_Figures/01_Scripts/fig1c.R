### LOAD IN LIBRARIES #########################################################
library(dplyr)
library(ggplot2)
library(nlme)
library(emmeans)
library(ggtext)
library(MuMIn)
library(reshape2)

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
metadata_its <- read.csv("atherton_sample_metadata_ITS_20250125.csv")
metadata_its <- distinct(metadata_its)

# read in sample names of cleaned data
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/02_Clean_Data/05_Transform_Data/ITS/Aitchison_Distance_Tables/")
all_aitch_its <- read.csv("atherton_ITS_allsampletypes_aitchisondistance_20240919.csv", 
                          row.names = 1)
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/02_Clean_Data/05_Transform_Data/16S/Aitchison_Distance_Tables/")
bg_aitch_16s <- read.csv("atherton_16S_allsampletypes_aitchisondistance_20240925.csv", 
                         row.names = 1)
leaf_aitch_16s <- read.csv("atherton_16S_leaf_aitchisondistance_20241112.csv", 
                           row.names = 1)

### FORMAT METADATA ###########################################################
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

# Normalize the data to be used in this figure
metadata_its$perc_sooty_mold <- log(metadata_its$perc_sooty_mold+1)
metadata_16s$perc_plant_pathogen <- log(metadata_16s$perc_plant_pathogen+1)

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE 1C #################################################################
# Select leaf data from ITS
metadata_leaf_its <- metadata_its[which(metadata_its$sample_type == "Leaf"),]

# For figure caption: number of trees' data used (n = 83)
print("Number of trees used for this analysis:")
length(unique(metadata_leaf_its$tree_id))

# Select leaf data from 16S
metadata_leaf_16s <- metadata_16s[which(metadata_16s$sample_type == "Leaf"),]

# For figure caption: number of trees' data used (n = 62)
print("Number of trees used for this analysis:")
length(unique(metadata_leaf_16s$tree_id))

# Take the appropriate guild abundances part of the metadata from ITS and 16S
sooty_mold <- metadata_leaf_its[,which(colnames(metadata_leaf_its) %in% 
                                         c("tree_dist_to_boston_km", 
                                           "perc_sooty_mold"))]
sooty_mold$Guild <- "Sooty Molds"
colnames(sooty_mold)[2] <- "abundance"
plant_path <- metadata_leaf_16s[,which(colnames(metadata_leaf_16s) %in% 
                                         c("tree_dist_to_boston_km", 
                                           "perc_plant_pathogen"))]
plant_path$Guild <- "Bacterial Plant Pathogens"
colnames(plant_path)[2] <- "abundance"

# Combine the data for just the guilds analyzed in this figure 
plot_data <- rbind(sooty_mold, plant_path)
plot_data$Guild <- factor(plot_data$Guild, 
                          levels = c("Bacterial Plant Pathogens", 
                                     "Sooty Molds"))

### STATISTICS SOOTY MOLDS ###
model_sootymolds <- lme(formula(perc_sooty_mold ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
                        random = ~1|plot_name/tree_id, na.action = na.exclude, 
                        method = "ML", data = metadata_leaf_its)
anova(model_sootymolds)
r.squaredGLMM(model_sootymolds)

### STATISTICS PLANT PATHOGENS ###
model_plantpath <- lme(formula(perc_plant_pathogen ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
                       random = ~1|plot_name/tree_id, na.action = na.exclude, 
                       method = "ML", data = metadata_leaf_16s)
anova(model_plantpath)
r.squaredGLMM(model_plantpath)

# Plot figure 1c
fig1c <- ggplot(plot_data, aes(x = log(tree_dist_to_boston_km+1), 
                               y = abundance, col = Guild, fill = Guild, )) + 
  xlab("log(Distance from Boston (km))") + 
  ylab("log(Functional group\nrelative abundance (%))") + 
  # this adds the individual data points
  geom_point(aes(x = log(tree_dist_to_boston_km+1),
                 y = abundance,
                 col = Guild), 
             alpha = 0.5, stroke = NA) +
  # this creates the standard error curve
  geom_smooth(method = "lm", alpha = 0.12, 
              aes(ymax = after_stat(y + se * sqrt(length(y))),
                  ymin = after_stat(y - se * sqrt(length(y))))) +
  # this creates the linear regression
  geom_smooth(method = "lm", se = FALSE, aes(x = log(tree_dist_to_boston_km+1),
                                             y = abundance,
                                             col = Guild, fill = Guild)) +
  theme_classic() + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), 
        legend.position = "none") +
  scale_fill_manual(values = rev(linear_pal)) + 
  scale_color_manual(values = rev(linear_pal)) + 
  labs(fill = "Microbial Guild", color = "Microbial Guild") +
  ggtitle("Plant pathogens in leaves")

fig1c
ggsave(paste0("atherton_figure1c_plantpathogensleaves_vs_dfb_", 
              format(Sys.Date(), "%Y%m%d"),".pdf"), device = cairo_pdf, 
       plot = fig1c, dpi = 300, width = 5.7, height = 3.5, units = "in")
