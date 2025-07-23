### LOAD IN LIBRARIES #########################################################
library(dplyr)
library(ggplot2)
library(nlme)
library(emmeans)
library(ggtext)
library(MuMIn)
library(reshape2)
library(ggsignif)
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

### FIGURE 3A #################################################################
# Make new x-axis variable: combination of urban and edge types
metadata_16s$site_edge <- paste0(metadata_16s$tree_urban_type, " ", 
                                 metadata_16s$tree_edge_type)
metadata_16s$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                               "Urban Street Tree", metadata_16s$site_edge)

# Make variable a factor and order it for the x-axis
metadata_16s$site_edge <- factor(metadata_16s$site_edge, 
                                 levels = c("Urban Street Tree", 
                                            "Urban Forest Edge", 
                                            "Urban Forest Interior", 
                                            "Rural Forest Edge", 
                                            "Rural Forest Interior"))

# Select the data for only the soil top layer
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == 
                                          "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 88)
print("Number of trees used in the soil analysis:")
length(unique(metadata_soil_16s$tree_id))

# Create statistical data frames
Results_anova <- NULL
Results_multcomp <- NULL

# Keep the data that will be plotted in this figure
data <- as.data.frame(metadata_soil_16s$soil_no3)

# Rename the column as it will appear in the figure
colnames(data) <- c("Available soil NO<sub>3</sub><sup>-</sup>")

# Take the mean and standard error of the data for each tree type
mean <- aggregate(data$`Available soil NO<sub>3</sub><sup>-</sup>` ~metadata_soil_16s$site_edge, 
                  FUN= "mean", na.action = na.omit)
se <- aggregate(data$`Available soil NO<sub>3</sub><sup>-</sup>` ~metadata_soil_16s$site_edge, 
                FUN= standard_error, na.action = na.omit)

# Reformat the data for the figure plotting
melt.mean <- melt(mean)
melt.se <- melt(se)

# Combine the mean and standard error for plotting
vdata <- cbind(melt.mean, melt.se[,ncol(melt.se)])

# Rename the columns
colnames(vdata) <- c("site_edge","Nutirent","Abundance","SE")

### STATISTICS SOIL NO3- ###
model_no3 <- lme(formula(soil_no3 ~ site_edge), 
                 random = ~1|plot_name/tree_id, na.action = na.exclude, 
                 method = "ML", data = metadata_soil_16s)
summary(model_no3)
summary(emmeans(model_no3, pairwise ~ site_edge))

# Create the Tukey groups
vdata$tukey <- c(
  "a", 
  "ab", 
  "b", 
  "ab", 
  "b"
)

# Reformat the x-axis tickmark labels
metadata_soil_16s$site_edge <- gsub(" ", "\n", metadata_soil_16s$site_edge)
metadata_soil_16s$site_edge <- factor(metadata_soil_16s$site_edge, 
                                      levels = c("Urban\nStreet\nTree", 
                                                 "Urban\nForest\nEdge",
                                                 "Urban\nForest\nInterior",
                                                 "Rural\nForest\nEdge",
                                                 "Rural\nForest\nInterior"))

vdata$site_edge <- gsub(" ", "\n", vdata$site_edge)
vdata$site_edge <- factor(vdata$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

# Visualize

fig3a <- ggplot() +
  # this creates the dotted line between points to track the average across tree types
  geom_line(data = vdata, 
            aes(x=site_edge, y=Abundance, group=Nutirent), 
            linetype = 2, color = "black") +
  # this creates the points that show the average abundance for each guild and tree type
  geom_point(data = vdata, 
             aes(x=site_edge, y=Abundance, color=site_edge, shape=Nutirent), 
             size=3) + 
  # this creates the standard error bars around the average abundance points
  geom_errorbar(data = vdata, 
                aes(x=site_edge, ymin=Abundance-SE, ymax=Abundance+SE, 
                    color=site_edge), width = 0.2, size = 0.3) +
  # this adds the Tukey groups labels defined above
  geom_text(data = vdata, direction = "y",
            aes(x=site_edge, y=Abundance + SE + 0.05, label = tukey), 
            size = 3, parse = TRUE) +
  # this adds the individual data points to the figure
  geom_point(data = metadata_soil_16s,
             aes(x = site_edge, y = soil_no3, color = site_edge),
             position = position_jitter(width = 0.2), size = 2, alpha = 0.2, 
             shape = 16, stroke = NA) +
  theme_classic() + labs (y="Soil inorganic NO<sub>3</sub><sup>-</sup> content (ug/g)", 
                          x="") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_markdown(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_markdown(size = 10, color = "black"), 
        legend.text = element_markdown(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_markdown(size = 10), legend.position = "none", 
        legend.background = element_rect(fill = NA)) +
  scale_color_manual(values = categorical_pal) + 
  guides(color = "none", shape = guide_legend(title = ""))+
  ggtitle("Available soil NO<sub>3</sub><sup>-</sup>")

fig3a
ggsave(paste0("atherton_figure3a_no3soil_vs_treetype_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig3a, dpi = 300, width = 6, height = 3.5, units = "in")
