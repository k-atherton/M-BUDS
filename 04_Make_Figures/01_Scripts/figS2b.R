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
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == 
                                          "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 88)
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

nh4 <- metadata_soil_16s[,which(colnames(metadata_soil_16s) %in% 
                                  c("site_edge", 
                                    "rindy_2019_nh4n_deposition.kgn_ha_yr"))]
no3 <- metadata_soil_16s[,which(colnames(metadata_soil_16s) %in% 
                                  c("site_edge", 
                                    "rindy_2019_no3n_deposition.kgn_ha_yr"))]

colnames(nh4)[1] <- "deposition"
colnames(no3)[1] <- "deposition"

nh4$Nutrient <- "NH<sub>4</sub><sup>+</sup> annual<br>throughfall"
no3$Nutrient <- "NO<sub>3</sub><sup>-</sup> annual<br>throughfall"

vdata <- rbind(nh4, no3)

# Format the x-axis tickmarks
vdata$site_edge <- gsub(" ", "\n", vdata$site_edge)
vdata$site_edge <- factor(vdata$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

# Order the functional guilds as a factor for plotting
vdata$Nutrient <- factor(vdata$Nutrient, 
                              levels = c("NH<sub>4</sub><sup>+</sup> annual<br>throughfall",
                                         "NO<sub>3</sub><sup>-</sup> annual<br>throughfall"))

# Summarize the data for adding the error bars, mean, and standard error
summary_data <- vdata %>%
  group_by(site_edge, Nutrient) %>%
  summarise(
    mean = mean(deposition, na.rm = TRUE),
    se = sd(deposition, na.rm = TRUE) / sqrt(n()),
    sd = sd(deposition, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot Figure S2b
figS2b <- ggplot(vdata, aes(x = site_edge, color = Nutrient)) +
  # this adds the jittered points
  geom_jitter(aes(y = deposition, group = Nutrient, stroke = NA), 
              size = 0.75,  alpha = 0.2, 
              position = position_jitterdodge(jitter.width = 0.5, 
                                              dodge.width = 0.75)) +
  # this adds the lines that follow the mean across tree types
  stat_summary(data = vdata,
               aes(y = deposition, group = Nutrient, color = Nutrient),
               fun = mean, geom = "line",
               position = position_dodge(width = 0.75), linetype = "dashed",
               linewidth = 0.8, alpha = 0.5) +
  # this adds the error bars
  geom_errorbar(data = summary_data,
                aes(x = site_edge, ymin = mean - sd, ymax = mean + sd,
                    group = Nutrient, color = Nutrient),
                position = position_dodge(width = 0.75), width = 0.1,
                alpha = 0.5) +
  # this adds the box with the mean and standard error
  geom_crossbar(data = summary_data,
                aes(x = site_edge, y = mean, ymin = mean - se, ymax = mean + se,
                    fill = Nutrient),
                position = position_dodge(width = 0.75), width = 0.5,
                alpha = 0.5, fatten = 1) +
  theme_classic() +
  labs(y="Atmospheric inorganic N\nin throughfall (kg/ha)", x = "", 
       color = "", fill = "") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color = "black"), 
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.text = element_markdown(size = 10),
        legend.title = element_text(size = 0),
        legend.key.size = unit(0.15, "in"),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(size = 10),
        legend.position = c(0.85, 0.9)) +
  scale_color_manual(values = linear_pal) +
  scale_fill_manual(values = linear_pal) +
  ggtitle("Annual atmospheric N throughfall")

figS2b
ggsave(paste0("atherton_figureS2b_nthroughfall_vs_treetype_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figS2b, dpi = 300, width = 6, height = 3.5, units = "in")
