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

### FIGURE S2B #################################################################
# Select the soil, top layer data from 16S
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 52)
print("Number of trees used in the root analysis:")
length(unique(metadata_soil_16s$tree_id))

# Format the x-axis tick marks
metadata_soil_16s$site_edge <- paste0(metadata_soil_16s$tree_urban_type, " ", 
                                      metadata_soil_16s$tree_edge_type)
metadata_soil_16s$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                                    "Urban Street Tree",
                                    metadata_soil_16s$site_edge)
metadata_soil_16s$site_edge <- factor(metadata_soil_16s$site_edge, 
                                      levels = c("Urban Street Tree", 
                                                 "Urban Forest Edge", 
                                                 "Urban Forest Interior", 
                                                 "Rural Forest Edge", 
                                                 "Rural Forest Interior"))

# Keep the inorganic N deposition data
data <- cbind(metadata_soil_16s$rindy_2019_nh4n_deposition.kgn_ha_yr, 
              metadata_soil_16s$rindy_2019_no3n_deposition.kgn_ha_yr)
colnames(data) <- c("NH<sub>4</sub><sup>+</sup> annual throughfall", 
                    "NO<sub>3</sub><sup>-</sup> annual throughfall")

# Take the mean and standard error of the data
mean <- aggregate(data~site_edge, data = metadata_soil_16s, FUN= "mean",
                  na.action = na.omit)
se <- aggregate(data~site_edge, data = metadata_soil_16s, FUN= standard_error,
                na.action = na.omit)

# Reformat the data for plotting
melt.mean <- melt(mean)
melt.se <- melt(se)

# Combine the mean and standard error into one data frame
vdata <- cbind(melt.mean, melt.se[,ncol(melt.se)])
colnames(vdata) <- c("site_edge","Nutrient","Abundance","SE")

# Reformat the x-axis tick marks for plotting
vdata$site_edge <- gsub(" ", "\n", vdata$site_edge)
vdata$site_edge <- factor(vdata$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

# Reorder the nutrients for plotting
vdata$Nutrient <- factor(vdata$Nutrient, 
                         levels = c("NH<sub>4</sub><sup>+</sup> annual throughfall", 
                                    "NO<sub>3</sub><sup>-</sup> annual throughfall"))

# Plot figure s2b
figS2b <- ggplot() +
  # This adds the dotted line tracking average abundance across groups
  geom_line(data = vdata, 
            aes(x=site_edge, y=Abundance, shape = Nutrient, group = Nutrient), 
            linetype = 2) +
  # This adds the average abundance point
  geom_point(data = metadata_soil_16s, 
             aes(x=site_edge, y=Abundance, color=site_edge, shape = Nutrient, 
                 group = Nutrient), size=5) + 
  # This adds the error bars around the average abundance point
  geom_errorbar(data = vdata, 
                aes(x=site_edge, ymin=Abundance-SE, ymax=Abundance+SE, 
                    color=site_edge), width = 0.2) +
  theme_classic() + labs (y="Atmospheric inorganic N deposition (kg/ha)", x="") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_markdown(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_markdown(size = 10), 
        legend.position = c(0.85, 1), 
        legend.background = element_rect(fill = NA)) +
  scale_color_manual(values = categorical_pal) + 
  scale_fill_manual(values = categorical_pal) +
  guides(color = "none", shape = guide_legend(title = ""))+
  ggtitle("Annual atmospheric N throughfall")

figS2b
ggsave(paste0("atherton_figureS2b_nthroughfall_vs_treetype_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figsS2b, dpi = 300, width = 6, height = 3.5, units = "in")
