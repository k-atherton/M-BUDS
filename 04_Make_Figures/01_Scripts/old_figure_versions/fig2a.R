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

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE 2A #################################################################
# Select leaf data from ITS
metadata_leaf_its <- metadata_its[which(metadata_its$sample_type == "Leaf"),]

# For figure caption: number of trees' data used (n = 83)
print("Number of trees used in the soil analysis:")
length(unique(metadata_leaf_its$tree_id))

# Select soil, top layer data from ITS 
metadata_soil_its <- metadata_its[which(metadata_its$sample_type == 
                                          "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 87)
print("Number of trees used in the soil analysis:")
length(unique(metadata_soil_its$tree_id))

# Select soil, top layer data from 16S
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == 
                                          "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 88)
print("Number of trees used in the soil analysis:")
length(unique(metadata_soil_16s$tree_id))

# Take the fungal guild abundances part of the metadata from soil and leaves
fungal_trait_leaf <- as.data.frame(metadata_leaf_its[,colnames(metadata_leaf_its) == 
                                         "perc_animal_pathogen"])
fungal_trait_soil <- as.data.frame(metadata_soil_its[,colnames(metadata_soil_its) == 
                                         "perc_animal_pathogen"])
bacterial_trait_soil <- as.data.frame(metadata_soil_16s[,colnames(metadata_soil_16s) == 
                                            "perc_zoonotic_pathogen"])

# Calculate the mean of the guild abundances for each tree type
mean_its_leaf <- aggregate(fungal_trait_leaf,
                           by=list(metadata_leaf_its$tree_urban_type,
                                   metadata_leaf_its$tree_edge_type), 
                           FUN= "mean")
colnames(mean_its_leaf)[3] <- "animal_leaf_mean"
mean_its_soil <- aggregate(fungal_trait_soil,
                           by=list(metadata_soil_its$tree_urban_type,
                                   metadata_soil_its$tree_edge_type), 
                           FUN= "mean")
colnames(mean_its_soil)[3] <- "animal_soil_mean"
mean_16s_soil <- aggregate(bacterial_trait_soil,
                           by=list(metadata_soil_16s$tree_urban_type,
                                   metadata_soil_16s$tree_edge_type), 
                           FUN= "mean")
colnames(mean_16s_soil)[3] <- "zoonotic_soil_mean"

# Calculate the standard error of the guild abundances for each tree type
se_its_leaf <- aggregate(fungal_trait_leaf,
                         by=list(metadata_leaf_its$tree_urban_type,
                                 metadata_leaf_its$tree_edge_type), 
                         FUN= standard_error)
colnames(se_its_leaf)[3] <- "animal_leaf_se"
se_its_soil <- aggregate(fungal_trait_soil,
                         by=list(metadata_soil_its$tree_urban_type,
                                 metadata_soil_its$tree_edge_type), 
                         FUN= standard_error)
colnames(se_its_soil)[3] <- "animal_soil_se"
se_16s_soil <- aggregate(bacterial_trait_soil,
                         by=list(metadata_soil_16s$tree_urban_type,
                                 metadata_soil_16s$tree_edge_type), 
                         FUN= standard_error)
colnames(se_16s_soil)[3] <- "zoonotic_soil_se"

# Combine the data for just the guilds analyzed in this figure (animal pathogens and human pathogens)
vdata <- merge(mean_its_leaf, mean_its_soil, 
               by = intersect(names(mean_its_leaf), names(mean_its_soil)))
vdata <- merge(vdata, mean_16s_soil, 
               by = intersect(names(vdata), names(mean_16s_soil)))
vdata <- merge(vdata, se_its_leaf, 
               by = intersect(names(vdata), names(se_its_leaf)))
vdata <- merge(vdata, se_its_soil, 
               by = intersect(names(vdata), names(se_its_soil)))
vdata <- merge(vdata, se_16s_soil, 
               by = intersect(names(vdata), names(se_16s_soil)))

# Name the columns 
colnames(vdata)[c(1,2)] <- c("Site", "Edge")

# Make a long version of the data for visualization
vdata_long <- as.data.frame(matrix(data = NA, nrow = 15, ncol = 5))

# Name the columns
colnames(vdata_long) <- c("Site", "Edge", "Guild", "Abundance", "SE")

# Since we have 3 guilds in the figure, repeat the site and edge information three times
vdata_long$Site <- rep(vdata$Site, 3)
vdata_long$Edge <- rep(vdata$Edge, 3)

# Since there are 5 tree types, repeat the guilds five times
vdata_long$Guild <- c(rep("Fungal animal pathogens in leaves", 5),
                      rep("Fungal animal pathogens in soil", 5),
                      rep("Bacterial zoonotic pathogens in soil", 5))

# Order the guilds for the figure
vdata_long$Guild <- factor(vdata_long$Guild, 
                           levels = c("Fungal animal pathogens in leaves", 
                                      "Bacterial zoonotic pathogens in soil",
                                      "Fungal animal pathogens in soil"))

# Input the abundance and standard error data for each guild
vdata_long$Abundance <- c(vdata$animal_leaf_mean, vdata$animal_soil_mean, 
                          vdata$zoonotic_soil_mean)
vdata_long$SE <- c(vdata$animal_leaf_se, vdata$animal_soil_se, 
                   vdata$zoonotic_soil_se)

# Create the x-axis categories
vdata_long$site_edge <- paste0(vdata_long$Site, " ", vdata_long$Edge)
vdata_long$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                             "Urban Street Tree",
                             vdata_long$site_edge)

# Make the tree type a factor with the order for the x-axis
vdata_long$site_edge <- factor(vdata_long$site_edge, 
                               levels = c("Urban Street Tree",
                                          "Urban Forest Edge", 
                                          "Urban Forest Interior", 
                                          "Rural Forest Edge", 
                                          "Rural Forest Interior"))

### STATISTICS LEAF ANIMAL PATHOGENS ###
model_leaf_animal <- lme(formula(perc_animal_pathogen ~ tree_pit_type), 
                         random = ~1|plot_name/tree_id, na.action = na.exclude, 
                         method = "ML", data = metadata_leaf_its)
summary(model_leaf_animal)
summary(emmeans(model_leaf_animal, pairwise ~ tree_pit_type))

### STATISTICS SOIL ANIMAL PATHOGENS ###
model_soil_animal <- lme(formula(perc_animal_pathogen ~ tree_pit_type), 
                         random = ~1|plot_name/tree_id, na.action = na.exclude, 
                         method = "ML", data = metadata_soil_its)
summary(model_soil_animal)
summary(emmeans(model_soil_animal, pairwise ~ tree_pit_type))

### STATISTICS SOIL HUMAN PATHOGENS ###
model_soil_zoonotic <- lme(formula(perc_zoonotic_pathogen ~ tree_pit_type), 
                           random = ~1|plot_name/tree_id, na.action = na.exclude, 
                           method = "ML", data = metadata_soil_16s)
summary(model_soil_zoonotic)
summary(emmeans(model_soil_zoonotic, pairwise ~ tree_pit_type))

# Make the text for adding Tukey groups
vdata_long$tukey <- c(
  "ab[\"\\u25CF\"]", # rural edge animal pathogens in leaves
  "b[\"\\u25CF\"]", # rural interior animal pathogens in leaves
  "b[\"\\u25CF\"]", # urban edge animal pathogens in leaves
  "ab[\"\\u25CF\"]", #urban interior animal pathogens in leaves
  "a[\"\\u25CF\"]", # street trees animal pathogens in leaves
  "ab[\"\\u25A0\"]", # rural edge animal pathogens in soil
  "ab[\"\\u25A0\"]", # rural interior animal pathogens in soil
  "b[\"\\u25A0\"]", # urban edge animal pathogens in soil
  "ab[\"\\u25A0\"]", # urban interior animal pathogens in soil
  "a[\"\\u25A0\"]", # street trees animal pathogens in soil
  "b[\"\\u25B2\"]", # rural edge zoonotic pathogens in soil
  "b[\"\\u25B2\"]", # rural interior zoonotic pathogens in soil
  "b[\"\\u25B2\"]", # urban edge zoonotic pathogens in soil
  "b[\"\\u25B2\"]",  #urban interior zoonotic pathogens in soil
  "a[\"\\u25B2\"]" # street trees zoonotic pathogens in soil
)

# Format labels for plotting
vdata_long$site_edge <- gsub(" ", "\n", vdata_long$site_edge)
vdata_long$site_edge <- factor(vdata_long$site_edge, 
                               levels = c("Urban\nStreet\nTree", 
                                          "Urban\nForest\nEdge",
                                          "Urban\nForest\nInterior",
                                          "Rural\nForest\nEdge",
                                          "Rural\nForest\nInterior"))


fig2a <- ggplot() +
  # this creates the dotted line between points to track the average across tree types
  geom_line(data = vdata_long, 
            aes(x=site_edge, y=Abundance, group = Guild), 
            linetype = 2) +
  # this creates the points that show the average abundance for each guild and tree type
  geom_point(data = vdata_long, 
             aes(x=site_edge, y=Abundance, color=site_edge, shape = Guild, 
                 group = Guild), size=3) + 
  # this creates the standard error bars around the average abundance points
  geom_errorbar(data = vdata_long, 
                aes(x=site_edge, ymin=Abundance-SE, ymax=Abundance+SE, 
                    color=site_edge), width = 0.2, size = 0.3) +
  # this adds the Tukey groups labels defined above
  geom_text_repel(data = vdata_long, family = "Arial", direction = "y",
                  aes(x=site_edge, y=Abundance, label = tukey), 
                  hjust = 1, nudge_x = -0.1, nudge_y = 0.5,   # Align left of points
                  size = 3, parse = TRUE, segment.color = NA) +
  theme_classic() + labs (y="Functional group relative abundance (%)", x="", shape = "") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 9.5, color="black"), 
        axis.title.y = element_text(size = 9.5), 
        axis.text.x = element_text(size = 9.5, color = "black"), 
        axis.text.y = element_text(size = 9.5, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 9.5), legend.position = c(0.5, 0.8)) +
  scale_color_manual(values = categorical_pal, guide = "none") + 
  # scale_fill_manual(values = pal, guide = "none") + 
  ggtitle("Fungal animal pathogens in soil and leaves and bacterial zoonotic pathogens in soil")

fig2a
ggsave(paste0("athterton_figure2a_pathogens_fragmentation_", 
              format(Sys.Date(), "%Y%m%d"),".pdf"), device = cairo_pdf, 
       plot = fig2a, dpi = 300, width = 6, height = 3.5, units = "in")
