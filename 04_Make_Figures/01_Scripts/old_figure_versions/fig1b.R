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
# Combine ITS soil horizon with sample type
metadata_its$sample_type<- paste0(metadata_its$soil_horizon, " ", 
                                  metadata_its$sample_type)
metadata_its$sample_type <- gsub("NA ", "", 
                                 metadata_its$sample_type)

# Replace ITS soil horizon category with depth measurement
metadata_its$sample_type <- gsub("O Soil", "Soil, 0-15 cm", 
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

# Normalize the data to be used in this figure
metadata_its$perc_wood_saprotroph <- log(metadata_its$perc_wood_saprotroph+1)
metadata_16s$perc_cellulolytic <- log(metadata_16s$perc_cellulolytic+1)
metadata_16s$perc_lignolytic <- log(metadata_16s$perc_lignolytic+1)

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE 1B ##################################################################
# Select soil, top layer data from ITS
metadata_soil_its <- metadata_its[which(metadata_its$sample_type == 
                                          "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 87)
print("Number of trees used for this analysis:")
length(unique(metadata_soil_its$tree_id))

# Select soil, top layer data from 16S
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == 
                                          "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 87)
print("Number of trees used for this analysis:")
length(unique(metadata_soil_its$tree_id))

# Take the fungal guild abundances part of the metadata from ITS and 16S
fungal_trait <- metadata_soil_its[, colnames(metadata_soil_its) == 
                                    "perc_wood_saprotroph"]
bacterial_trait <- metadata_soil_16s[, colnames(metadata_soil_16s) %in% 
                                       c("perc_cellulolytic", 
                                         "perc_lignolytic")]

# Calculate the mean of the guild abundances for each tree type
mean_its <- aggregate(fungal_trait,
                      by=list(metadata_soil_its$tree_urban_type,
                              metadata_soil_its$tree_edge_type), FUN= "mean")
colnames(mean_its)[3] <- "woodsap_mean"
mean_16s <- aggregate(bacterial_trait,
                      by=list(metadata_soil_16s$tree_urban_type,
                              metadata_soil_16s$tree_edge_type), FUN= "mean")
colnames(mean_16s)[c(3,4)] <- c("cell_mean", "lig_mean")

# Calculate the standard error of the guild abundances for each tree type
se_its <- aggregate(fungal_trait,
                    by=list(metadata_soil_its$tree_urban_type,
                            metadata_soil_its$tree_edge_type), 
                    FUN= standard_error)
colnames(se_its)[3] <- "woodsap_se"
se_16s <- aggregate(bacterial_trait,
                    by=list(metadata_soil_16s$tree_urban_type,
                            metadata_soil_16s$tree_edge_type), 
                    FUN= standard_error)
colnames(se_16s)[c(3,4)] <- c("cell_se", "lig_se")

# Combine the data for just the guilds analyzed in this figure 
# (wood saprotrophs, cellulolytic bacteria, and lignolytic bacteria)
vdata <- merge(mean_its, mean_16s, by = intersect(names(mean_its), 
                                                  names(mean_16s)))
vdata <- merge(vdata, se_its, by = intersect(names(vdata), 
                                                  names(se_its)))
vdata <- merge(vdata, se_16s, by = intersect(names(vdata), 
                                             names(se_16s)))

# Name the columns 
colnames(vdata)[1:2] <- c("Site", "Edge")

# Make a long version of the data for visualization
vdata_long <- as.data.frame(matrix(data = NA, nrow = 15, ncol = 5))

# Name the columns
colnames(vdata_long) <- c("Site", "Edge", "Guild", "Abundance", "SE")

# Since we have 3 guilds in the figure, repeat the site and edge information three times
vdata_long$Site <- c(rep(vdata$Site, 3))
vdata_long$Edge <- rep(vdata$Edge, 3)

# Since there are 5 tree types, repeat the guilds five times
vdata_long$Guild <- c(rep("Fungal wood decomposers", 5), 
                      rep("Cellulolytic bacteria", 5), 
                      rep("Lignolytic bacteria", 5))

# Input the abundance and standard error data for each guild
vdata_long$Abundance <- c(vdata$woodsap_mean, vdata$cell_mean, vdata$lig_mean)
vdata_long$SE <- c(vdata$woodsap_se, vdata$cell_se, vdata$lig_se)

# Create the x-axis categories
vdata_long$site_edge <- paste0(vdata_long$Site, " ", vdata_long$Edge)
vdata_long$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                             "Urban Street Tree", vdata_long$site_edge)

# Make the tree type a factor with the order for the x-axis
vdata_long$site_edge <- factor(vdata_long$site_edge, 
                               levels = c("Urban Street Tree", 
                                          "Urban Forest Edge",
                                          "Urban Forest Interior", 
                                          "Rural Forest Edge", 
                                          "Rural Forest Interior"))

# Make guild a factor for plotting the legend in order
vdata_long$Guild <- factor(vdata_long$Guild, 
                           levels = c("Cellulolytic bacteria", 
                                      "Lignolytic bacteria", 
                                      "Fungal wood decomposers"))

### STATISTICS WOOD SAPROTROPHS ###
model_woodsap <- lme(formula(perc_wood_saprotroph ~ tree_pit_type), 
                     random = ~1|plot_name/tree_id, na.action = na.exclude, 
                     method = "ML", data = metadata_soil_its)
summary(model_woodsap)
summary(emmeans(model_woodsap, pairwise ~ tree_pit_type))

### STATISTICS CELLULOLYTIC BACTERIA ###
model_cell <- lme(formula(perc_cellulolytic ~ tree_pit_type), 
                  random = ~1|plot_name/tree_id, na.action = na.exclude, 
                  method = "ML", data = metadata_soil_16s)
summary(model_cell)
summary(emmeans(model_cell, pairwise ~ tree_pit_type))


### STATISTICS LIGNOLYGIC BACTERIA ###
model_lig <- lme(formula(perc_lignolytic ~ tree_pit_type), 
                 random = ~1|plot_name/tree_id, na.action = na.exclude, 
                 method = "ML", data = metadata_soil_16s)
summary(model_lig)
summary(emmeans(model_lig, pairwise ~ tree_pit_type))

# Define Tukey groups for the figure
vdata_long$tukey <- c("b[\"\\u25A0\"]", # rural edge fungal wood saps
                      "b[\"\\u25A0\"]", # rural interior fungal wood saps
                      "b[\"\\u25A0\"]", # urban edge fungal wood saps
                      "b[\"\\u25A0\"]", # urban interior fungal wood saps
                      "a[\"\\u25A0\"]", # street tree fungal wood saps 
                      "b[\"\\u25CF\"]", # rural edge cellulolytic
                      "b[\"\\u25CF\"]", # rural interior cellulolytic
                      "b[\"\\u25CF\"]", # urban edge cellulolytic
                      "b[\"\\u25CF\"]", # urban interior cellulolytic
                      "a[\"\\u25CF\"]", # street tree cellulolytic
                      "b[\"\\u25B2\"]", # rural edge lignolytic
                      "b[\"\\u25B2\"]", # rural interior lignolytic
                      "b[\"\\u25B2\"]", # urban edge lignolytic
                      "b[\"\\u25B2\"]", # urban interior lignolytic
                      "a[\"\\u25B2\"]" # street tree lignolytic
)

# Format labels for plotting
vdata_long$site_edge <- gsub(" ", "\n", as.character(vdata_long$site_edge))
vdata_long$site_edge <- factor(vdata_long$site_edge, 
                               levels = c("Urban\nStreet\nTree", 
                                          "Urban\nForest\nEdge",
                                          "Urban\nForest\nInterior",
                                          "Rural\nForest\nEdge",
                                          "Rural\nForest\nInterior"))

# Plot Figure 1b
fig1b <- ggplot(vdata_long) +
  # this creates the dotted line between points to track the average across tree types
  geom_line(aes(x=site_edge, y=Abundance, group = Guild), 
            linetype = 2) +
  # this creates the points that show the average abundance for each guild and tree type
  geom_point(aes(x=site_edge, y=Abundance, color=site_edge, shape = Guild, 
                 group = Guild), size=3) + 
  # this creates the standard error bars around the average abundance points
  geom_errorbar(aes(x=site_edge, ymin=Abundance-SE, ymax=Abundance+SE, 
                    color=site_edge), width = 0.2, size = 0.3) +
  # this adds the Tukey groups labels defined above
  # geom_text(family = "Arial",
  #           aes(x = site_edge, y = Abundance, label = tukey, hjust = 2.2, 
  #               vjust = 0.5, # Adjust text to the left of the points
  #               nudge_x = -0.5), size = 3, parse = TRUE) +
  geom_text_repel(family = "Arial", direction = "y", segment.color = NA,
            aes(x = site_edge, y = Abundance, label = tukey, hjust = 2.1, 
                vjust = 1), size = 3, parse = TRUE) +
  theme_classic() + 
  labs (y="log(Functional group\nrelative abundance (%))", x="", shape = "") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), legend.position = c(0.5, 0.8)) +
  scale_color_manual(values = categorical_pal, guide = "none") + 
  #scale_fill_manual(values = categorical_pal, guide = "none") + 
  ggtitle("Wood decomposers in soil")

fig1b
ggsave(paste0("atherton_figure1b_wooddecomposersoil_vs_treetype_", 
              format(Sys.Date(), "%Y%m%d"),".pdf"), device = cairo_pdf, 
       plot = fig1b, dpi = 300, width = 5.75, height = 3.5, units = "in")
