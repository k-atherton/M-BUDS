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
metadata_soil_its <- metadata_its[which(metadata_its$sample_type == "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 87)
print("Number of trees used in the soil analysis:")
length(unique(metadata_soil_its$tree_id))

# Select soil, top layer data from 16S
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 88)
print("Number of trees used in the soil analysis:")
length(unique(metadata_soil_16s$tree_id))

# Take the fungal guild abundances part of the metadata from soil and leaves
metadata_leaf_its <- metadata_leaf_its[which(colnames(metadata_leaf_its) %in% 
                                               c("sample_name", "tree_id", 
                                                 "plot_name", "tree_edge_type", 
                                                 "tree_urban_type", 
                                                 "perc_animal_pathogen"))]
metadata_soil_its <- metadata_soil_its[which(colnames(metadata_soil_its) %in%
                                               c("sample_name", "tree_id", 
                                                 "plot_name", "tree_edge_type", 
                                                 "tree_urban_type", 
                                                 "perc_animal_pathogen"))]
metadata_soil_16s <- metadata_soil_16s[which(colnames(metadata_soil_16s) %in%
                                               c("sample_name", "tree_id", 
                                                 "plot_name", "tree_edge_type", 
                                                 "tree_urban_type", 
                                                 "perc_zoonotic_pathogen"))]

# Create the x-axis categories
metadata_leaf_its$site_edge <- paste0(metadata_leaf_its$tree_urban_type, " ", 
                                      metadata_leaf_its$tree_edge_type)
metadata_leaf_its$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                                    "Urban Street Tree", 
                                    metadata_leaf_its$site_edge)

metadata_soil_its$site_edge <- paste0(metadata_soil_its$tree_urban_type, " ", 
                                      metadata_soil_its$tree_edge_type)
metadata_soil_its$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                                    "Urban Street Tree", 
                                    metadata_soil_its$site_edge)

metadata_soil_16s$site_edge <- paste0(metadata_soil_16s$tree_urban_type, " ", 
                                      metadata_soil_16s$tree_edge_type)
metadata_soil_16s$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                                    "Urban Street Tree", 
                                    metadata_soil_16s$site_edge)

# Rename the columns
colnames(metadata_leaf_its) <- c("sample_name", "tree_id", "plot_name", 
                                 "tree_edge_type", "tree_urban_type", 
                                 "rel_abund", "site_edge")
colnames(metadata_soil_its) <- c("sample_name", "tree_id", "plot_name",
                                 "tree_edge_type", "tree_urban_type", 
                                 "rel_abund", "site_edge")
colnames(metadata_soil_16s) <- c("sample_name", "tree_id", "plot_name", 
                                 "tree_edge_type", "tree_urban_type", 
                                 "rel_abund", "site_edge")

# Add the name of the functional guild
metadata_leaf_its$guild <- "Fungal animal pathogens in leaves"
metadata_soil_its$guild <- "Fungal animal pathogens in soil"
metadata_soil_16s$guild <- "Bacterial zoonotic pathogens in soil"

# Combine the soil and leaf datasets for plotting
vdata <- rbind(metadata_leaf_its, metadata_soil_its, metadata_soil_16s)

# Make the tree type a factor with the order for the x-axis
vdata$site_edge <- factor(vdata$site_edge, levels = c("Urban Street Tree", 
                                                      "Urban Forest Edge",
                                                      "Urban Forest Interior",
                                                      "Rural Forest Edge", 
                                                      "Rural Forest Interior"))

# Make guild a factor for plotting the legend in order
vdata$guild <- factor(vdata$guild, levels = c("Fungal animal pathogens in leaves", 
                                              "Bacterial zoonotic pathogens in soil",
                                              "Fungal animal pathogens in soil"))

### STATISTICS FUNGAL ANIMAL PATHOGENS IN LEAVES ###
model_animal_leaves <- lme(formula(rel_abund ~ site_edge), 
                           random = ~1|plot_name/tree_id, na.action = na.exclude, 
                           method = "ML", data = metadata_leaf_its)
summary(model_animal_leaves)
summary(emmeans(model_animal_leaves, pairwise ~ site_edge))

### STATISTICS FUNGAL ANIMAL PATHOGENS IN SOIL ###
model_animal_soil <- lme(formula(rel_abund ~ site_edge), 
                         random = ~1|plot_name/tree_id, na.action = na.exclude, 
                         method = "ML", data = metadata_soil_its)
summary(model_animal_soil)
summary(emmeans(model_animal_soil, pairwise ~ site_edge))

### STATISTICS BACTERIAL ZOONOTIC PATHOGENS ###
model_zoonotic <- lme(formula(rel_abund ~ site_edge), 
                      random = ~1|plot_name/tree_id, na.action = na.exclude, 
                      method = "ML", data = metadata_soil_16s)
summary(model_zoonotic)
summary(emmeans(model_zoonotic, pairwise ~ site_edge))

# Make the text for adding Tukey groups and record the location of the group id on the y-axis
vdata$tukey <- NA
vdata$location <- NA

for(i in 1:nrow(vdata)){
  guild <- as.character(vdata$guild[i])
  site_edge <- as.character(vdata$site_edge[i])
  
  if(guild == "Fungal animal pathogens in leaves"){
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "ab"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "ab"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    }
  } else if(guild == "Bacterial zoonotic pathogens in soil"){
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    }
  } else if(guild == "Fungal animal pathogens in soil"){
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "ab"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "ab"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "ab"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 2
    }
  }
}

# Format labels for plotting
vdata$site_edge <- gsub(" ", "\n", vdata$site_edge)
vdata$site_edge <- factor(vdata$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

# Plot Figure 2b
figS3f <- ggplot() +
  # this creates the points that show the average abundance for each guild and tree type
  geom_boxplot(data = vdata, outlier.shape = NA,
               aes(x=site_edge, y=rel_abund, fill=guild), alpha = 0.5) + 
  geom_jitter(data=vdata, aes(x = site_edge, y = rel_abund, color = guild, 
                              group = guild, alpha = 0.75, stroke = NA),
              position = position_jitterdodge(jitter.width = 0.5, 
                                              dodge.width = 0.75)) +
  # this adds the Tukey groups labels defined above
  geom_text(data=vdata, family = "Arial",
            aes(x = site_edge, y = location, color = guild, group = guild, 
                label = tukey), position = position_dodge(width = 0.75),
            size = 3, parse = TRUE) +
  theme_classic() + labs (y="Functional group\nrelative abundance (%)", x="", color = "", 
                          fill = "") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 1), 
        legend.key.size = unit(0.15, "in"),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(size = 10), legend.position = c(0.8, 1)) +
  scale_color_manual(values = categorical_pal[c(1,3,4)], guide = "none") + 
  scale_fill_manual(values = categorical_pal[c(1,3,4)]) + 
  ggtitle("Fungal animal pathogens in soil and leaves and bacterial zoonotic pathogens in soil")

figS3f
ggsave(paste0("athterton_figureS3f_pathogens_fragmentation_boxplot_", 
              format(Sys.Date(), "%Y%m%d"),".pdf"), device = cairo_pdf, 
       plot = figS3f, dpi = 300, width = 5.75, height = 3.5, units = "in")
