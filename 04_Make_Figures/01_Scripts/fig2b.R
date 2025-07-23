### LOAD IN LIBRARIES #########################################################
library(dplyr)
library(ggplot2)
library(nlme)
library(emmeans)
library(ggtext)
library(MuMIn)
library(reshape2)
library(ggrepel)

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

### FIGURE 2B ##################################################################
# Select the leaf data from 16S
metadata_leaf_16s <- metadata_16s[which(metadata_16s$sample_type == "Leaf"),]

# For figure caption: number of trees' data used (n = 62)
print("Number of trees used in the soil analysis:")
length(unique(metadata_leaf_16s$tree_id))

### STATISTICS ZOONOTIC PATHOGENS ###
model_zoonotic <- lme(formula(perc_zoonotic_pathogen ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
                      random = ~1|plot_name/tree_id, na.action = na.exclude, 
                      method = "ML", data = metadata_leaf_16s)
anova(model_zoonotic)
r.squaredGLMM(model_zoonotic)

fig2b <- ggplot(metadata_leaf_16s) + 
  # this makes the individual points
  geom_point(aes(x = log(tree_dist_to_boston_km+1),
                 y = perc_zoonotic_pathogen,
                 col = linear_pal[2]), 
             alpha = 0.5, stroke = NA) +
  # this makes the standard error area around the regression
  geom_smooth(method = "lm", alpha = 0.12, 
              aes(x = log(tree_dist_to_boston_km+1),
                  y = perc_zoonotic_pathogen,
                  col = linear_pal[2], 
                  fill = linear_pal[2],
                  ymax = after_stat(y + se * sqrt(length(y))),
                  ymin = after_stat(y - se * sqrt(length(y))))) +
  # this makes the linear regression line
  geom_smooth(method = "lm", se = FALSE, aes(x = log(tree_dist_to_boston_km+1),
                                             y = perc_zoonotic_pathogen,
                                             col = linear_pal[2], 
                                             fill = linear_pal[2])) +
  xlab("log(Distance from Boston (km))") + 
  ylab("Bacterial zoonotic pathogen\nrelative abundance (%)") + 
  theme_classic() + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 9.5, color="black"), 
        axis.title.y = element_text(size = 9.5), 
        axis.text.x = element_text(size = 9.5, color = "black"), 
        axis.text.y = element_text(size = 9.5, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 9.5), legend.position = "none") +
  scale_fill_manual(values = linear_pal[2]) + 
  scale_color_manual(values = linear_pal[2]) + 
  labs(fill = "Microbial Guild", color = "Microbial Guild") +
  ggtitle("Bacterial zoonotic pathogens in leaves")

fig2b
ggsave(paste0("atherton_figure2b_humanpathogensleaves_vs_dfb_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig2b, dpi = 300, width = 6, height = 3.5, units = "in")
