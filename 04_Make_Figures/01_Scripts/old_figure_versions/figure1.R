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
# Combine ITS soil horizon with sample type
metadata_its$sample_type<- paste0(metadata_its$soil_horizon, " ", 
                                  metadata_its$sample_type)
metadata_its$sample_type <- gsub("NA ", "", 
                                 metadata_its$sample_type)

# Replace ITS soil horizon category with depth measurement
metadata_its$sample_type <- gsub("O Soil", "Soil, 0-15 cm", 
                                 metadata_its$sample_type)
metadata_its$sample_type <- gsub("M Soil", "Soil, 15-30 cm", 
                                 metadata_its$sample_type)

metadata_its$sample_type <- gsub("O Root", "Roots, 0-15 cm", 
                                 metadata_its$sample_type)
metadata_its$sample_type <- gsub("M Root", "Roots, 15-30 cm", 
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

# Format the ITS environmental data
metadata_its$soil_percent_moisture <- gsub("%", "", 
                                           metadata_its$soil_percent_moisture)
metadata_its$soil_percent_organic_matter <- gsub("%", "",
                                                 metadata_its$soil_percent_organic_matter)
metadata_its$soil_percent_moisture <- as.numeric(metadata_its$soil_percent_moisture)
metadata_its$soil_percent_organic_matter <- as.numeric(metadata_its$soil_percent_organic_matter)
metadata_its$tree_distance_from_edge <- as.numeric(metadata_its$tree_distance_from_edge)

# Combine 16S soil horizon with sample type
metadata_16s$sample_type<- paste0(metadata_16s$soil_horizon, " ",
                                      metadata_16s$sample_type)
metadata_16s$sample_type <- gsub("NA ", "", 
                                     metadata_16s$sample_type)

# Replace 16S soil horizon category with depth measurement
metadata_16s$sample_type <- gsub("O Soil", "Soil, 0-15 cm",
                                     metadata_16s$sample_type)
metadata_16s$sample_type <- gsub("M Soil", "Soil, 15-30 cm",
                                     metadata_16s$sample_type)

metadata_16s$sample_type <- gsub("O Root", "Roots, 0-15 cm",
                                     metadata_16s$sample_type)
metadata_16s$sample_type <- gsub("M Root", "Roots, 15-30 cm",
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

# Format the 16S environmental data
metadata_16s$soil_percent_moisture <- gsub("%", "",
                                           metadata_16s$soil_percent_moisture)
metadata_16s$soil_percent_organic_matter <- gsub("%", "",
                                                 metadata_16s$soil_percent_organic_matter)
metadata_16s$soil_percent_moisture <- as.numeric(metadata_16s$soil_percent_moisture)
metadata_16s$soil_percent_organic_matter <- as.numeric(metadata_16s$soil_percent_organic_matter)
metadata_16s$tree_distance_from_edge <- as.numeric(metadata_16s$tree_distance_from_edge)

# Normalize the data to be used in this figure
metadata_its$perc_epiphyte <- log(metadata_its$perc_epiphyte+1)
metadata_its$perc_plant_pathogen <- log(metadata_its$perc_plant_pathogen+1)
metadata_its$perc_sooty_mold <- log(metadata_its$perc_sooty_mold+1)
metadata_its$perc_wood_saprotroph <- log(metadata_its$perc_wood_saprotroph+1)

metadata_16s$perc_cellulolytic <- log(metadata_16s$perc_cellulolytic+1)
metadata_16s$perc_lignolytic <- log(metadata_16s$perc_lignolytic+1)
metadata_16s$perc_plant_pathogen <- log(metadata_16s$perc_plant_pathogen+1)

# Make a list of fungal guild names
fungal_trait_names <- c("perc_algal_parasite", "perc_animal_endosymbiont", 
                        "perc_animal_parasite", 
                        "perc_amf", "perc_arthropod_associated", 
                        "perc_dung_saprotroph", "perc_ecm", "perc_epiphyte",
                        "perc_foliar_endophyte", "perc_lichen_parasite",
                        "perc_lichenized", "perc_litter_saprotroph",
                        "perc_moss_symbiont", "perc_mycoparasite",
                        "perc_nectar_tap_saprotroph", 
                        "perc_plant_pathogen", 
                        "perc_pollen_saprotroph", "perc_root_endophyte",
                        "perc_soil_saprotroph", "perc_sooty_mold",
                        "perc_wood_saprotroph", 
                        "perc_animal_pathogen", 
                        "perc_human_pathogen", "perc_opportunistic_pathogen",
                        "perc_root_associated_pathogen", 
                        "perc_leaf_fruit_seed_pathogen", "perc_wood_pathogen",
                        "shannon", "shannon_nopath", "perc_all_saprotroph", 
                        "perc_all_pathotroph", "perc_all_symbiont")

# Make a list of bacterial guild names
bacterial_trait_names <- c("perc_c_fixation", "perc_cellulolytic", 
                           "perc_chitinolytic", "perc_lignolytic",
                           "perc_carbon_monoxide_oxidation", "perc_n_fixation",
                           "perc_copiotroph", "perc_denitrification",
                           "perc_dissim_nitrate_reduction", 
                           "perc_hydrocarbon_degradation",
                           "perc_sulfonate_desulfurization",
                           "perc_other_p_cycling",
                           "perc_oxidize_reduced_sulfur", "perc_methanotroph",
                           "perc_oligotroph", "perc_partial_nitrification", 
                           "perc_archaea", "perc_plant_pathogen", 
                           "perc_animal_pathogen", "perc_zoonotic_pathogen", 
                           "shannon", "shannon.x", "shannon_nopath")

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE 1A #################################################################
# Select soil, top layer data from ITS
metadata_soil_its <- metadata_its[which(metadata_its$sample_type %in% c("Soil, 0-15 cm")),]

# For figure caption: number of trees' data used (n = 87)
print("Number of trees used in the soil analysis:")
length(unique(metadata_soil_its$tree_id))

# Select leaf data from ITS
metadata_leaf_its <- metadata_its[which(metadata_its$sample_type %in% c("Leaf")),]

# For figure caption: number of trees' data used (n = 83)
print("Number of trees used in the leaf analysis:")
length(unique(metadata_leaf_its$tree_id))

# Take the fungal guild abundances part of the metadata from soil and leaves
fungal_trait_soil <- metadata_soil_its[, colnames(metadata_soil_its) %in% 
                                         fungal_trait_names]
fungal_trait_leaf <- metadata_leaf_its[, colnames(metadata_leaf_its) %in% 
                                         fungal_trait_names]

# Calculate the mean of the guild abundances for each tree type
mean_soil <- aggregate(fungal_trait_soil,
                       by=list(metadata_soil_its$tree_urban_type,
                               metadata_soil_its$tree_edge_type), FUN= "mean")
mean_leaf <- aggregate(fungal_trait_leaf,
                       by=list(metadata_leaf_its$tree_urban_type,
                               metadata_leaf_its$tree_edge_type), FUN= "mean")

# Calculate the standard error of the guild abundances for each tree type
se_soil <- aggregate(fungal_trait_soil,
                     by=list(metadata_soil_its$tree_urban_type,
                             metadata_soil_its$tree_edge_type), 
                     FUN= standard_error)
se_leaf <- aggregate(fungal_trait_leaf,
                     by=list(metadata_leaf_its$tree_urban_type,
                             metadata_leaf_its$tree_edge_type), 
                     FUN= standard_error)

# Combine the data for just the guilds analyzed in this figure (ECM and epiphytes)
vdata <- cbind(mean_soil[,1:2], mean_soil$perc_ecm, se_soil$perc_ecm, 
               mean_leaf$perc_epiphyte, se_leaf$perc_epiphyte, 
               mean_leaf$perc_foliar_endophyte, se_leaf$perc_foliar_endophyte)

# Name the columns 
colnames(vdata) <- c("Site", "Edge", "ecm", "ecm_se", "epiphyte", "epiphtye_se", 
                     "endophyte", "endophyte_se")

# Make a long version of the data for visualization
vdata_long <- as.data.frame(matrix(data = NA, nrow = 10, ncol = 5))

# Name the columns
colnames(vdata_long) <- c("Site", "Edge", "Guild", "Abundance", "SE")

# Since we have 2 guilds in the figure, repeat the site and edge information twice
vdata_long$Site <- rep(vdata$Site, 2)
vdata_long$Edge <- rep(vdata$Edge, 2)

# Since there are 5 tree types, repeat the guilds five times
vdata_long$Guild <- c(rep("ECM in soil", 5), 
                      rep("Epiphytes in leaves", 5))

# Input the abundance data for each guild, multiply the epiphyte abundance by 20
# for visualization -- this will be adjusted for in the secondary y-axis
vdata_long$Abundance <- c(vdata$ecm, 20*vdata$epiphyte)

# Input the standard error data for each guild
vdata_long$SE <- c(vdata$ecm_se, vdata$epiphyte_se)

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

# Match the x-axis labels for the remaining metadata files so that the individual points can be plotted
metadata_soil_its$site_edge <- paste0(metadata_soil_its$tree_urban_type, " ", metadata_soil_its$tree_edge_type)
metadata_soil_its$site_edge[which(metadata_soil_its$site_edge == "Urban Street Tree Urban Street Tree")] <- "Urban Street Tree"
metadata_leaf_its$site_edge <- paste0(metadata_leaf_its$tree_urban_type, " ", metadata_leaf_its$tree_edge_type)
metadata_leaf_its$site_edge[which(metadata_leaf_its$site_edge == "Urban Street Tree Urban Street Tree")] <- "Urban Street Tree"

# Make guild a factor for plotting the legend in order
vdata_long$Guild <- factor(vdata_long$Guild, levels = c("ECM in soil", 
                                                        "Epiphytes in leaves"))

### STATISTICS ECM ###
model_ecm <- lme(formula(perc_ecm ~ site_edge), 
                 random = ~1|plot_name/tree_id, na.action = na.exclude, 
                 method = "ML", data = metadata_soil_its)
summary(model_ecm)
summary(emmeans(model_ecm, pairwise ~ site_edge))

### STATISTICS EPIPHYTES ###
model_epiphyte <- lme(formula(perc_epiphyte ~ site_edge), 
                      random = ~1|plot_name/tree_id, na.action = na.exclude, 
                      method = "ML", data = metadata_leaf_its)
summary(model_epiphyte)
summary(emmeans(model_epiphyte, pairwise ~ site_edge))


# Make the text for adding Tukey groups
vdata_long$tukey <- c(
  "ab[\"\\u25CF\"]", # ecm rural edge
  "a[\"\\u25CF\"]", # ecm rural interior
  "a[\"\\u25CF\"]", # ecm urban edge
  "a[\"\\u25CF\"]", # ecm urban interior
  "b[\"\\u25CF\"]", # ecm street tree
  "a[\"\\u25B2\"]", # epiphyte rural edge
  "a[\"\\u25B2\"]", # epiphyte rural interior
  "a[\"\\u25B2\"]", # epiphyte urban edge
  "a[\"\\u25B2\"]", # epiphyte urban interior
  "b[\"\\u25B2\"]" # epiphyte street tree
)

# Format labels for plotting
metadata_soil_its$site_edge <- gsub(" ", "\n", metadata_soil_its$site_edge)
metadata_soil_its$site_edge <- factor(metadata_soil_its$site_edge, 
                                      levels = c("Urban\nStreet\nTree", 
                                                 "Urban\nForest\nEdge",
                                                 "Urban\nForest\nInterior",
                                                 "Rural\nForest\nEdge",
                                                 "Rural\nForest\nInterior"))
metadata_leaf_its$site_edge <- gsub(" ", "\n", metadata_leaf_its$site_edge)
metadata_leaf_its$site_edge <- factor(metadata_leaf_its$site_edge, 
                                      levels = c("Urban\nStreet\nTree", 
                                                 "Urban\nForest\nEdge",
                                                 "Urban\nForest\nInterior",
                                                 "Rural\nForest\nEdge",
                                                 "Rural\nForest\nInterior"))

vdata_long$site_edge <- gsub(" ", "\n", vdata_long$site_edge)
vdata_long$site_edge <- factor(vdata_long$site_edge, 
                               levels = c("Urban\nStreet\nTree", 
                                          "Urban\nForest\nEdge",
                                          "Urban\nForest\nInterior",
                                          "Rural\nForest\nEdge",
                                          "Rural\nForest\nInterior"))

# Plot Figure 1a
fig1a <- ggplot() +
  # this creates the dotted line between points to track the average across tree types
  geom_line(data = vdata_long,
            aes(x=site_edge, y=Abundance, group = Guild),
            linetype = 2) +
  # this creates the points that show the average abundance for each guild and tree type
  geom_point(data = vdata_long, 
             aes(x=site_edge, y=Abundance, color=site_edge, shape = Guild, 
                 group = Guild), size=5) + 
  # this creates the standard error bars around the average abundance points
  geom_errorbar(data = vdata_long, 
                aes(x=site_edge, ymin=Abundance-SE, ymax=Abundance+SE, 
                    color=site_edge), width = 0.2) +
  # this adds the Tukey groups labels defined above
  geom_text(data=vdata_long, family = "Arial",
            aes(x = site_edge, y = Abundance + SE + 3, label = tukey), size = 3, 
            parse = TRUE) +
  # this adds the individual data points to the figure
  geom_point(data = metadata_soil_its,
             aes(x = site_edge, y = perc_ecm, color = site_edge),
             position = position_jitter(width = 0.2), size = 3, alpha = 0.2,
             shape = 16, stroke = NA) +
  geom_point(data = metadata_leaf_its,
             aes(x = site_edge, y = perc_epiphyte * 20, color = site_edge),
             position = position_jitter(width = 0.2), size = 3, alpha = 0.2,
             shape = 2) +
  theme_classic() + labs (y="ECM relative abundance (%)", x="", shape = "") + 
  # this creates the secondary axis, dividing by 20 to adjust for the multiplication by 20 above
  scale_y_continuous(sec.axis = sec_axis(~./20, 
                                         name = "log(Epiphyte relative abundance (%))")) +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), legend.position = c(0.5, 0.1)) +
  scale_color_manual(values = categorical_pal, guide = "none") + 
  ggtitle("Tree-associated fungal symbionts")

fig1a
ggsave(paste0("atherton_figure1a_ecmsoil_epiphytesleaves_vs_treetype_",
              format(Sys.Date(), "%Y%m%d"),".pdf"), 
       plot = fig1a, dpi = 300, width = 6, height = 3.5, units = "in", 
       device = cairo_pdf)

rm(list = c("mean_leaf", "mean_soil", "metadata_leaf_its", "metadata_soil_its",
            "model_ecm", "model_epiphyte", "se_leaf", "se_soil", "vdata", 
            "vdata_long"))

### FIGURE 1B ##################################################################
# Select soil, top layer data from ITS
metadata_soil_its <- metadata_its[which(metadata_its$sample_type == "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 87)
print("Number of trees used for this analysis:")
length(unique(metadata_soil_its$tree_id))

# Select soil, top layer data from 16S
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 87)
print("Number of trees used for this analysis:")
length(unique(metadata_soil_its$tree_id))

# Take the fungal guild abundances part of the metadata from ITS and 16S
fungal_trait <- metadata_soil_its[, colnames(metadata_soil_its) %in% 
                                    fungal_trait_names]
bacterial_trait <- metadata_soil_16s[, colnames(metadata_soil_16s) %in% 
                                       bacterial_trait_names]

# Calculate the mean of the guild abundances for each tree type
mean_its <- aggregate(fungal_trait,
                     by=list(metadata_soil_its$tree_urban_type,
                             metadata_soil_its$tree_edge_type), FUN= "mean")
mean_16s <- aggregate(bacterial_trait,
                     by=list(metadata_soil_16s$tree_urban_type,
                             metadata_soil_16s$tree_edge_type), FUN= "mean")

# Calculate the standard error of the guild abundances for each tree type
se_its <- aggregate(fungal_trait,
                   by=list(metadata_soil_its$tree_urban_type,
                           metadata_soil_its$tree_edge_type), FUN= standard_error)
se_16s <- aggregate(bacterial_trait,
                   by=list(metadata_soil_16s$tree_urban_type,
                           metadata_soil_16s$tree_edge_type), FUN= standard_error)

# Combine the data for just the guilds analyzed in this figure 
# (wood saprotrophs, cellulolytic bacteria, and lignolytic bacteria)
vdata <- cbind(mean_its[,1:2], mean_its$perc_wood_saprotroph, 
               se_its$perc_wood_saprotroph, mean_16s$perc_cellulolytic, 
               se_16s$perc_cellulolytic, mean_16s$perc_lignolytic, 
               se_16s$perc_lignolytic)

# Name the columns
colnames(vdata) <- c("Site", "Edge", "woodsap", "woodsap_se", "cell", "cell_se", 
                     "lig", "lig_se")

# Make a long version of the data for visualization
vdata_long <- as.data.frame(matrix(data = NA, nrow = 15, ncol = 5))

# Name the columns
colnames(vdata_long) <- c("Site", "Edge", "Guild", "Abundance", "SE")

# Since we have 2 guilds in the figure, repeat the site and edge information three times
vdata_long$Site <- rep(vdata$Site, 3)
vdata_long$Edge <- rep(vdata$Edge, 3)

# Since there are 5 tree types, repeat the guilds five times
vdata_long$Guild <- c(rep("Fungal wood decomposers", 5), 
                      rep("Cellulolytic bacteria", 5), 
                      rep("Lignolytic bacteria", 5))

# Input the abundance and standard error data for each guild
vdata_long$Abundance <- c(vdata$woodsap, vdata$cell, vdata$lig)
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

# Match the x-axis labels for the remaining metadata files so that the individual points can be plotted
metadata_soil_its$site_edge <- paste0(metadata_soil_its$tree_urban_type, " ", metadata_soil_its$tree_edge_type)
metadata_soil_its$site_edge[which(metadata_soil_its$site_edge == "Urban Street Tree Urban Street Tree")] <- "Urban Street Tree"
metadata_soil_16s$site_edge <- paste0(metadata_soil_16s$tree_urban_type, " ", metadata_soil_16s$tree_edge_type)
metadata_soil_16s$site_edge[which(metadata_soil_16s$site_edge == "Urban Street Tree Urban Street Tree")] <- "Urban Street Tree"

# Make guild a factor for plotting the legend in order
vdata_long$Guild <- factor(vdata_long$Guild, levels = c("Cellulolytic bacteria", 
                                                        "Lignolytic bacteria", 
                                                        "Fungal wood decomposers"))

### STATISTICS WOOD SAPROTROPHS ###
model_woodsap <- lme(formula(perc_wood_saprotroph ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_its)
summary(model_woodsap)
summary(emmeans(model_woodsap, pairwise ~ site_edge))

### STATISTICS CELLULOLYTIC BACTERIA ###
model_cell <- lme(formula(perc_cellulolytic ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_16s)
summary(model_cell)
summary(emmeans(model_cell, pairwise ~ site_edge))


### STATISTICS LIGNOLYGIC BACTERIA ###
model_lig <- lme(formula(perc_lignolytic ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_16s)
summary(model_lig)
summary(emmeans(model_lig, pairwise ~ site_edge))

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
metadata_soil_its$site_edge <- gsub(" ", "\n", metadata_soil_its$site_edge)
metadata_soil_its$site_edge <- factor(metadata_soil_its$site_edge, 
                               levels = c("Urban\nStreet\nTree", 
                                          "Urban\nForest\nEdge",
                                          "Urban\nForest\nInterior",
                                          "Rural\nForest\nEdge",
                                          "Rural\nForest\nInterior"))

metadata_soil_16s$site_edge <- gsub(" ", "\n", metadata_soil_16s$site_edge)
metadata_soil_16s$site_edge <- factor(metadata_soil_16s$site_edge, 
                               levels = c("Urban\nStreet\nTree", 
                                          "Urban\nForest\nEdge",
                                          "Urban\nForest\nInterior",
                                          "Rural\nForest\nEdge",
                                          "Rural\nForest\nInterior"))

vdata_long$site_edge <- gsub(" ", "\n", vdata_long$site_edge)
vdata_long$site_edge <- factor(vdata_long$site_edge, 
                               levels = c("Urban\nStreet\nTree", 
                                          "Urban\nForest\nEdge",
                                          "Urban\nForest\nInterior",
                                          "Rural\nForest\nEdge",
                                          "Rural\nForest\nInterior"))

# Plot Figure 1b
fig1b <- ggplot() +
  # this creates the dotted line between points to track the average across tree types
  geom_line(data = vdata_long, 
            aes(x=site_edge, y=Abundance, group = Guild), 
            linetype = 2) +
  # this creates the points that show the average abundance for each guild and tree type
  geom_point(data = vdata_long, 
             aes(x=site_edge, y=Abundance, color=site_edge, shape = Guild, 
                 group = Guild), size=5) + 
  # this creates the standard error bars around the average abundance points
  geom_errorbar(data = vdata_long, 
                aes(x=site_edge, ymin=Abundance-SE, ymax=Abundance+SE, 
                    color=site_edge), width = 0.2) +
  # this adds the Tukey groups labels defined above
  geom_text(data=vdata_long, family = "Arial",
            aes(x = site_edge, y = Abundance, label = tukey, hjust = 2.2, 
                vjust = 0.5, # Adjust text to the left of the points
                nudge_x = -0.5), size = 3, parse = TRUE) +
  # geom_text_repel(data=vdata_long, family = "Arial", direction = "y", segment.color = NA,
  #           aes(x = site_edge, y = Abundance, label = tukey, hjust = 2.1, vjust = 1        # Adjust text to the left of the points
  #           ), size = 3, parse = TRUE) +
  # this adds the individual data points to the figure
  geom_point(data = metadata_soil_its,
             aes(x = site_edge, y = perc_wood_saprotroph, color = site_edge),
             position = position_jitter(width = 0.2), size = 2, alpha = 0.2, shape = 15) +
  geom_point(data = metadata_soil_16s,
             aes(x = site_edge, y = perc_cellulolytic, color = site_edge),
             position = position_jitter(width = 0.2), size = 2, alpha = 0.2, shape = 1) +
  geom_point(data = metadata_soil_16s,
             aes(x = site_edge, y = perc_lignolytic, color = site_edge),
             position = position_jitter(width = 0.2), size = 2, alpha = 0.2, shape = 2) +
  theme_classic() + labs (y="log(Functional group relative abundance (%))", x="", shape = "") +
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

rm(list = c("mean_16s", "mean_its", "metadata_soil_16s", "metadata_soil_its",
            "model_cell", "model_lig", "model_woodsap", "se_16s", "se_its", 
            "vdata", "vdata_long"))

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
sooty_mold <- metadata_leaf_its[,which(colnames(metadata_leaf_its) %in% c("tree_dist_to_boston_km", "perc_sooty_mold"))]
sooty_mold$Guild <- "Sooty Molds"
colnames(sooty_mold)[2] <- "abundance"
plant_path <- metadata_leaf_16s[,which(colnames(metadata_leaf_16s) %in% c("tree_dist_to_boston_km", "perc_plant_pathogen"))]
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
fig1c <- ggplot(plot_data) + 
  xlab("log(Distance from Boston (km))") + 
  ylab("log(Functional group relative abundance (%))") + 
  # this creates the standard error curve
  geom_smooth(method = "lm", alpha = 0.12, aes(x = log(tree_dist_to_boston_km+1),
                                               y = abundance,
                                               col = Guild, fill = Guild,
                                               ymax = after_stat(y + se * sqrt(length(y))),
                                               ymin = after_stat(y - se * sqrt(length(y))))) +
  # this creates the linear regression
  geom_smooth(method = "lm", se = FALSE, aes(x = log(tree_dist_to_boston_km+1),
                                             y = abundance,
                                             col = Guild, fill = Guild)) +
  # this adds the individual data points
  geom_point(aes(x = log(tree_dist_to_boston_km+1),
                 y = abundance,
                 col = Guild), 
             alpha = 0.5, stroke = NA, size = 2) +
  theme_classic() + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), legend.position = "none") +
  scale_fill_manual(values = rev(linear_pal)) + 
  scale_color_manual(values = rev(linear_pal)) + 
  labs(fill = "Microbial Guild", color = "Microbial Guild") +
  ggtitle("Plant pathogens in leaves")
fig1c

ggsave(paste0("atherton_figure1c_plantpathogensleaves_vs_dfb_", 
              format(Sys.Date(), "%Y%m%d"),".pdf"), device = cairo_pdf, 
       plot = fig1c, dpi = 300, width = 5.75, height = 3.5, units = "in")

rm(list = c("sooty_mold", "plant_path", "metadata_leaf_its", "metadata_leaf_16s",
            "model_sootymolds", "model_plantpath", "plot_data"))

### FIGURE 1D #################################################################
# Select soil, top layer data from ITS
metadata_soil_its <- metadata_its[which(metadata_its$sample_type == "Soil, 0-15 cm"),]

# Select soil, top layer data from 16S
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]

# Select leaf data from ITS
metadata_leaf_its <- metadata_its[which(metadata_its$sample_type == "Leaf"),]

# Select leaf data from 16S
metadata_leaf_16s <- metadata_16s[which(metadata_16s$sample_type == "Leaf"),]

# Take the functional guild abundances part of the metadata from ITS and 16S
func_soil_its <- metadata_soil_its[,colnames(metadata_soil_its) %in%
                                             c(fungal_trait_names)]
func_soil_16s <- metadata_soil_16s[,colnames(metadata_soil_16s) %in% c(bacterial_trait_names)]

func_leaf_its <- metadata_leaf_its[,colnames(metadata_leaf_its) %in%
                                             c(fungal_trait_names)]
func_leaf_16s <- metadata_leaf_16s[,colnames(metadata_leaf_16s) %in% c(bacterial_trait_names)]

# Keep the fungal guilds to be analyzed in this figure
func_soil_its <- func_soil_its[,which(colnames(func_soil_its) %in% c("perc_ecm", "perc_wood_saprotroph"))]
func_leaf_its <- func_leaf_its[,which(colnames(func_leaf_its) %in% c("perc_epiphyte", "perc_sooty_mold"))]

# Rename the column names for the figure
colnames(func_soil_its) <- c("ECM","Wood decomposers")
colnames(func_soil_its) <- paste0(colnames(func_soil_its), " in soil")
colnames(func_leaf_its) <- c("Fungal epiphytes", "Sooty molds")
colnames(func_leaf_its) <- paste0(colnames(func_leaf_its), " in leaves")

# Keep bacterial guilds to be analyzed in this figure
func_soil_16s <- func_soil_16s[,which(colnames(func_soil_16s) %in% c("perc_cellulolytic", "perc_lignolytic"))]
func_leaf_16s <- as.data.frame(func_leaf_16s[,which(colnames(func_leaf_16s) %in% c("perc_plant_pathogen"))])

# Rename the column names for the figure
colnames(func_soil_16s) <- c("Cellulolytic bacteria", "Lignolytic bacteria")
colnames(func_soil_16s) <- paste0(colnames(func_soil_16s), " in soil")
colnames(func_leaf_16s) <- c("Bacterial plant pathogens")
colnames(func_leaf_16s) <- paste0(colnames(func_leaf_16s), " in leaves")

# Keep the environmental factors for the figure
metadata_soil_heatmap_its <- metadata_soil_its[,c(17,18,20,23,30:36,38,49,53:55,59:61,63:65,67:76)]
metadata_soil_heatmap_16s <- metadata_soil_16s[,c(17,18,20,23,30:36,38,49,53:55,59:61,63:65,67:76)]
metadata_leaf_heatmap_its <- metadata_leaf_its[,c(17,18,20,23,30:36,38,49,53:55,59:61,63:65,67:76)]
metadata_leaf_heatmap_16s <- metadata_leaf_16s[,c(17,18,20,23,30:36,38,49,53:55,59:61,63:65,67:76)]

# Rename the columns for the figure
colnames(metadata_soil_heatmap_its) <- c("Tree distance from edge", 
                                         "Tree distance from Boston", 
                                         "Tree DBH", "Tree growth rate",
                                         "Soil temperature", 
                                         "Soil bulk density", 
                                         "Soil pH", "Soil moisture", 
                                         "Soil organic matter", 
                                         "Soil available NH<sub>4</sub><sup>+</sup>",
                                         "Soil available NO<sub>3</sub><sup>-</sup>", 
                                         "Normalized root biomass", 
                                         "Soil soluble salt",
                                         "NO concentration",
                                         "NO<sub>2</sub> concentration",
                                         "NO<sub>x</sub> concentration",
                                         "NH<sub>4</sub><sup>+</sup> annual throughfall",
                                         "NO<sub>3</sub><sup>-</sup> annual throughfall",
                                         "Total inorganic N annual throughfall",
                                         "O<sub>3</sub> concentration",
                                         "Foliar N %",
                                         "Foliar C:N",
                                         "Ca<sup>2+</sup> annual wet deposition",
                                         "Cl<sup>-</sup> annual wet deposition",
                                         "H<sup>+</sup> annual wet deposition",
                                         "K<sup>+</sup> annual wet deposition",
                                         "Mg<sup>2+</sup> annual wet deposition",
                                         "Na<sup>+</sup> annual wet deposition",
                                         "NH<sub>4</sub><sup>+</sup> annual wet deposition",
                                         "NO<sub>3</sub><sup>-</sup> annual wet deposition",
                                         "SO<sub>4</sub><sup>2-</sup> annual wet deposition",
                                         "Total inorganic N annual wet deposition")

# Rename the columns
colnames(metadata_soil_heatmap_16s) <- c("Tree distance from edge", 
                                         "Tree distance from Boston", 
                                         "Tree DBH", "Tree growth rate",
                                         "Soil temperature", 
                                         "Soil bulk density", 
                                         "Soil pH", "Soil moisture", 
                                         "Soil organic matter", 
                                         "Soil available NH<sub>4</sub><sup>+</sup>",
                                         "Soil available NO<sub>3</sub><sup>-</sup>", 
                                         "Normalized root biomass", 
                                         "Soil soluble salt",
                                         "NO concentration",
                                         "NO<sub>2</sub> concentration",
                                         "NO<sub>x</sub> concentration",
                                         "NH<sub>4</sub><sup>+</sup> annual throughfall",
                                         "NO<sub>3</sub><sup>-</sup> annual throughfall",
                                         "Total inorganic N annual throughfall",
                                         "O<sub>3</sub> concentration",
                                         "Foliar N %",
                                         "Foliar C:N",
                                         "Ca<sup>2+</sup> annual wet deposition",
                                         "Cl<sup>-</sup> annual wet deposition",
                                         "H<sup>+</sup> annual wet deposition",
                                         "K<sup>+</sup> annual wet deposition",
                                         "Mg<sup>2+</sup> annual wet deposition",
                                         "Na<sup>+</sup> annual wet deposition",
                                         "NH<sub>4</sub><sup>+</sup> annual wet deposition",
                                         "NO<sub>3</sub><sup>-</sup> annual wet deposition",
                                         "SO<sub>4</sub><sup>2-</sup> annual wet deposition",
                                         "Total inorganic N annual wet deposition")

# Rename the columns
colnames(metadata_leaf_heatmap_its) <- c("Tree distance from edge", 
                                         "Tree distance from Boston", 
                                         "Tree DBH", "Tree growth rate",
                                         "Soil temperature", 
                                         "Soil bulk density", 
                                         "Soil pH", "Soil moisture", 
                                         "Soil organic matter", 
                                         "Soil available NH<sub>4</sub><sup>+</sup>",
                                         "Soil available NO<sub>3</sub><sup>-</sup>", 
                                         "Normalized root biomass", 
                                         "Soil soluble salt",
                                         "NO concentration",
                                         "NO<sub>2</sub> concentration",
                                         "NO<sub>x</sub> concentration",
                                         "NH<sub>4</sub><sup>+</sup> annual throughfall",
                                         "NO<sub>3</sub><sup>-</sup> annual throughfall",
                                         "Total inorganic N annual throughfall",
                                         "O<sub>3</sub> concentration",
                                         "Foliar N %",
                                         "Foliar C:N",
                                         "Ca<sup>2+</sup> annual wet deposition",
                                         "Cl<sup>-</sup> annual wet deposition",
                                         "H<sup>+</sup> annual wet deposition",
                                         "K<sup>+</sup> annual wet deposition",
                                         "Mg<sup>2+</sup> annual wet deposition",
                                         "Na<sup>+</sup> annual wet deposition",
                                         "NH<sub>4</sub><sup>+</sup> annual wet deposition",
                                         "NO<sub>3</sub><sup>-</sup> annual wet deposition",
                                         "SO<sub>4</sub><sup>2-</sup> annual wet deposition",
                                         "Total inorganic N annual wet deposition")

# Rename the columns
colnames(metadata_leaf_heatmap_16s) <- c("Tree distance from edge", 
                                         "Tree distance from Boston", 
                                         "Tree DBH", "Tree growth rate",
                                         "Soil temperature", 
                                         "Soil bulk density", 
                                         "Soil pH", "Soil moisture", 
                                         "Soil organic matter", 
                                         "Soil available NH<sub>4</sub><sup>+</sup>",
                                         "Soil available NO<sub>3</sub><sup>-</sup>", 
                                         "Normalized root biomass", 
                                         "Soil soluble salt",
                                         "NO concentration",
                                         "NO<sub>2</sub> concentration",
                                         "NO<sub>x</sub> concentration",
                                         "NH<sub>4</sub><sup>+</sup> annual throughfall",
                                         "NO<sub>3</sub><sup>-</sup> annual throughfall",
                                         "Total inorganic N annual throughfall",
                                         "O<sub>3</sub> concentration",
                                         "Foliar N %",
                                         "Foliar C:N",
                                         "Ca<sup>2+</sup> annual wet deposition",
                                         "Cl<sup>-</sup> annual wet deposition",
                                         "H<sup>+</sup> annual wet deposition",
                                         "K<sup>+</sup> annual wet deposition",
                                         "Mg<sup>2+</sup> annual wet deposition",
                                         "Na<sup>+</sup> annual wet deposition",
                                         "NH<sub>4</sub><sup>+</sup> annual wet deposition",
                                         "NO<sub>3</sub><sup>-</sup> annual wet deposition",
                                         "SO<sub>4</sub><sup>2-</sup> annual wet deposition",
                                         "Total inorganic N annual wet deposition")

# Make coefficient matrices
# Soil ITS
Result_soil_its <- NULL
for (i in 1:ncol(func_soil_its)){
  estimate_soil_its <- NULL
  for (j in 1:ncol(metadata_soil_heatmap_its)){
    cor_soil_its <- cor.test(func_soil_its[,i],metadata_soil_heatmap_its[,j], method="p")
    estimate_soil_its <- c(estimate_soil_its, cor_soil_its$estimate)
  }
  Result_soil_its <- rbind(Result_soil_its, estimate_soil_its)
}
colnames(Result_soil_its) <- colnames(metadata_soil_heatmap_its)
rownames(Result_soil_its) <- colnames(func_soil_its) 

# Leaf ITS
Result_leaf_its <- NULL
for (i in 1:ncol(func_leaf_its)){
  estimate_leaf_its <- NULL
  for (j in 1:ncol(metadata_leaf_heatmap_its)){
    cor_leaf_its <- cor.test(func_leaf_its[,i],metadata_leaf_heatmap_its[,j], method="p")
    estimate_leaf_its <- c(estimate_leaf_its, cor_leaf_its$estimate)
  }
  Result_leaf_its <- rbind(Result_leaf_its, estimate_leaf_its)
}
colnames(Result_leaf_its) <- colnames(metadata_leaf_heatmap_its)
rownames(Result_leaf_its) <- colnames(func_leaf_its) 

# Soil 16S
Result_soil_16s <- NULL
for (i in 1:ncol(func_soil_16s)){
  estimate_soil_16s <- NULL
  for (j in 1:ncol(metadata_soil_heatmap_16s)){
    cor_soil_16s <- cor.test(func_soil_16s[,i],metadata_soil_heatmap_16s[,j], method="p")
    estimate_soil_16s <- c(estimate_soil_16s, cor_soil_16s$estimate)
  }
  Result_soil_16s <- rbind(Result_soil_16s, estimate_soil_16s)
}
colnames(Result_soil_16s) <- colnames(metadata_soil_heatmap_16s)
rownames(Result_soil_16s) <- colnames(func_soil_16s) 

# Leaf 16S
Result_leaf_16s <- NULL
estimate_leaf_16s <- NULL
for (j in 1:ncol(metadata_leaf_heatmap_16s)){
  cor_leaf_16s <- cor.test(func_leaf_16s[,1],metadata_leaf_heatmap_16s[,j], method="p")
  estimate_leaf_16s <- c(estimate_leaf_16s, cor_leaf_16s$estimate)
}
Result_leaf_16s <- rbind(Result_leaf_16s, estimate_leaf_16s)
colnames(Result_leaf_16s) <- colnames(metadata_leaf_heatmap_16s)
rownames(Result_leaf_16s) <- "Bacterial plant pathogens in leaves" 

# Make Pval matrices
# Soil ITS
Result.Pval_soil_its <- NULL
for (i in 1:ncol(func_soil_its)){
  Pval_soil_its <- NULL
  for (j in 1:ncol(metadata_soil_heatmap_its)){
    cor_soil_its <- cor.test(func_soil_its[,i],metadata_soil_heatmap_its[,j], method="p")
    Pval_soil_its <- c(Pval_soil_its, cor_soil_its$p.value)
  }
  Result.Pval_soil_its <- rbind (Result.Pval_soil_its, Pval_soil_its)
}
colnames(Result.Pval_soil_its) <- colnames(metadata_soil_heatmap_its)
rownames(Result.Pval_soil_its) <- colnames(func_soil_its) 

# Leaf ITS
Result.Pval_leaf_its <- NULL
for (i in 1:ncol(func_leaf_its)){
  Pval_leaf_its <- NULL
  for (j in 1:ncol(metadata_leaf_heatmap_its)){
    cor_leaf_its <- cor.test(func_leaf_its[,i],metadata_leaf_heatmap_its[,j], method="p")
    Pval_leaf_its <- c(Pval_leaf_its, cor_leaf_its$p.value)
  }
  Result.Pval_leaf_its <- rbind (Result.Pval_leaf_its, Pval_leaf_its)
}
colnames(Result.Pval_leaf_its) <- colnames(metadata_leaf_heatmap_its)
rownames(Result.Pval_leaf_its) <- colnames(func_leaf_its) 

# Soil 16S
Result.Pval_soil_16s <- NULL
for (i in 1:ncol(func_soil_16s)){
  Pval_soil_16s <- NULL
  for (j in 1:ncol(metadata_soil_heatmap_16s)){
    cor_soil_16s <- cor.test(func_soil_16s[,i],metadata_soil_heatmap_16s[,j], method="p")
    Pval_soil_16s <- c(Pval_soil_16s, cor_soil_16s$p.value)
  }
  Result.Pval_soil_16s <- rbind (Result.Pval_soil_16s, Pval_soil_16s)
}
colnames(Result.Pval_soil_16s) <- colnames(metadata_soil_heatmap_16s)
rownames(Result.Pval_soil_16s) <- colnames(func_soil_16s) 

# Leaf 16S
Result.Pval_leaf_16s <- NULL
Pval_leaf_16s <- NULL
for (j in 1:ncol(metadata_leaf_heatmap_16s)){
  cor_leaf_16s <- cor.test(func_leaf_16s[,1],metadata_leaf_heatmap_16s[,j], method="p")
  Pval_leaf_16s <- c(Pval_leaf_16s, cor_leaf_16s$p.value)
}
Result.Pval_leaf_16s <- rbind (Result.Pval_leaf_16s, Pval_leaf_16s)
colnames(Result.Pval_leaf_16s) <- colnames(metadata_leaf_heatmap_16s)
rownames(Result.Pval_leaf_16s) <- "Bacterial plant pathogens in leaves"

# Make Asterisk matrices
# Soil ITS
Result.Pval.combined_soil_its <- Result.Pval_soil_its
Result.asterisk_soil_its <- Result.Pval.combined_soil_its
for (k in 1:ncol(Result.Pval.combined_soil_its)){
  for (i in 1:nrow(Result.Pval.combined_soil_its)){
    if(!is.na(Result.Pval.combined_soil_its[i,k])){
      if (Result.Pval.combined_soil_its[i,k] < 0.001/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_soil_its[i,k] <- "***"
      } else if (Result.Pval.combined_soil_its[i,k] < 0.01/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_soil_its[i,k] <- "**"
      } else if (Result.Pval.combined_soil_its[i,k] < 0.05/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_soil_its[i,k] <- "*"
      } else {
        Result.asterisk_soil_its[i,k]<- " "
      }
    } 
  }
}

# Leaf ITS
Result.Pval.combined_leaf_its <- Result.Pval_leaf_its
Result.asterisk_leaf_its <- Result.Pval.combined_leaf_its
for (k in 1:ncol(Result.Pval.combined_leaf_its)){
  for (i in 1:nrow(Result.Pval.combined_leaf_its)){
    if(!is.na(Result.Pval.combined_leaf_its[i,k])){
      if (Result.Pval.combined_leaf_its[i,k] < 0.001/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_leaf_its[i,k] <- "***"
      } else if (Result.Pval.combined_leaf_its[i,k] < 0.01/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_leaf_its[i,k] <- "**"
      } else if (Result.Pval.combined_leaf_its[i,k] < 0.05/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_leaf_its[i,k] <- "*"
      } else {
        Result.asterisk_leaf_its[i,k]<- " "
      }
    } 
  }
}

# Soil 16S
Result.Pval.combined_soil_16s <- Result.Pval_soil_16s
Result.asterisk_soil_16s <- Result.Pval.combined_soil_16s
for (k in 1:ncol(Result.Pval.combined_soil_16s)){
  for (i in 1:nrow(Result.Pval.combined_soil_16s)){
    if(!is.na(Result.Pval.combined_soil_16s[i,k])){
      if (Result.Pval.combined_soil_16s[i,k] < 0.001/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_soil_16s[i,k] <- "***"
      } else if (Result.Pval.combined_soil_16s[i,k] < 0.01/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_soil_16s[i,k] <- "**"
      } else if (Result.Pval.combined_soil_16s[i,k] < 0.05/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_soil_16s[i,k] <- "*"
      } else {
        Result.asterisk_soil_16s[i,k]<- " "
      }
    } 
  }
}

# Leaf 16S
Result.Pval.combined_leaf_16s <- Result.Pval_leaf_16s
Result.asterisk_leaf_16s <- Result.Pval.combined_leaf_16s
for (k in 1:ncol(Result.Pval.combined_leaf_16s)){
  for (i in 1:nrow(Result.Pval.combined_leaf_16s)){
    if(!is.na(Result.Pval.combined_leaf_16s[i,k])){
      if (Result.Pval.combined_leaf_16s[i,k] < 0.001/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_leaf_16s[i,k] <- "***"
      } else if (Result.Pval.combined_leaf_16s[i,k] < 0.01/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_leaf_16s[i,k] <- "**"
      } else if (Result.Pval.combined_leaf_16s[i,k] < 0.05/(7*32)){ # multiple comparisons correction -- 7 functional groups * 32 environmental factors
        Result.asterisk_leaf_16s[i,k] <- "*"
      } else {
        Result.asterisk_leaf_16s[i,k]<- " "
      }
    } 
  }
}

# Combine the p values with the correlation values
Result.combined_soil_its <- Result_soil_its
melt_soil_its <- melt(Result.combined_soil_its)
melt.asterisk_soil_its <- melt(Result.asterisk_soil_its)

Result.combined_leaf_its <- Result_leaf_its
melt_leaf_its <- melt(Result.combined_leaf_its)
melt.asterisk_leaf_its <- melt(Result.asterisk_leaf_its)

Result.combined_soil_16s <- Result_soil_16s
melt_soil_16s <- melt(Result.combined_soil_16s)
melt.asterisk_soil_16s <- melt(Result.asterisk_soil_16s)

Result.combined_leaf_16s <- Result_leaf_16s
melt_leaf_16s <- melt(Result.combined_leaf_16s)
melt.asterisk_leaf_16s <- melt(Result.asterisk_leaf_16s)

melt_all <- rbind(melt_soil_its, melt_soil_16s, melt_leaf_its, melt_leaf_16s)
melt_asterisk_all <- rbind(melt.asterisk_soil_its, melt.asterisk_soil_16s, 
                           melt.asterisk_leaf_its, melt.asterisk_leaf_16s)

# Combine the correlation values with the p-values
vdata_its <- cbind(melt_all, melt_asterisk_all$value)

# Rename the columns
colnames(vdata_its) <- c("microbe", "environment", "coef", "Pval")

# Order the microbial guilds
vdata_its$microbe <- factor(vdata_its$microbe, 
                            levels = c("ECM in soil",
                                       "Fungal epiphytes in leaves",
                                       "Wood decomposers in soil", 
                                       "Cellulolytic bacteria in soil", 
                                       "Lignolytic bacteria in soil",
                                       "Sooty molds in leaves",
                                       "Bacterial plant pathogens in leaves"))

# Make an order of the environmental factors based on their correlation with the first microbial group
order_of_env_factors <- vdata_its[which(vdata_its$microbe == "ECM in soil"),]
vdata_its$environment <- factor(vdata_its$environment, 
                                levels = order_of_env_factors$environment[order(order_of_env_factors$coef)])

# plot figure 1d
figure1d <- ggplot(vdata_its) +
  # this makes the heatmap
  geom_tile(aes(x=microbe, y=environment, fill=coef)) +
  # this adds the asterisks
  geom_text(aes(x=microbe, y=environment, label = Pval), size = 3) +
  labs(x="",y="") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.3)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.3)) + 
  scale_fill_gradient2(low = "#003049", mid = "#FDF0D5", high = "#C00000", 
                       name = "Correlation\ncoefficient") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_markdown(size = 8), 
        axis.text.y = element_markdown(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8), 
        plot.title = element_text(size = 8)) + 
  ggtitle("") + coord_fixed()

figure1d
ggsave(paste0("atherton_figure1d_functionalguilds_vs_environment_heatmap_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), plot = figure1d, 
       device = cairo_pdf, width = 5, height = 10, units = "in", dpi = 300)
