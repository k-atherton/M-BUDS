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
fig2a <- ggplot() +
  stat_summary(data = vdata, aes(x = site_edge, y = rel_abund, group = guild, 
                                 color = guild),
               fun = median,
               geom = "line",
               position = position_dodge(width = 0.75),
               linetype = "dashed",
               linewidth = 0.8) +
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
  theme_classic() + labs (y="Functional group relative abundance (%)", x="", color = "", 
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
  scale_color_manual(values = categorical_pal[c(1,3,5)], guide = "none") + 
  scale_fill_manual(values = categorical_pal[c(1,3,5)]) + 
  ggtitle("Fungal animal pathogens in soil and leavers and bacterial zoonotic pathogens in soil")

fig2a
ggsave(paste0("athterton_figure2a_pathogens_fragmentation_boxplot_", 
              format(Sys.Date(), "%Y%m%d"),".pdf"), device = cairo_pdf, 
       plot = fig2a, dpi = 300, width = 6, height = 3.5, units = "in")

rm(list = c("metadata_leaf_its", "metadata_soil_its", "metadata_soil_16s", 
            "model_animal_leaves", "model_animal_soil", "model_zoonotic", 
            "vdata"))

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
  # this makes the standard error area around the regression
  geom_smooth(method = "lm", alpha = 0.12, aes(x = log(tree_dist_to_boston_km+1),
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
  # this makes the individual points
  geom_point(aes(x = log(tree_dist_to_boston_km+1),
                 y = perc_zoonotic_pathogen,
                 col = linear_pal[2]), 
             alpha = 0.5, stroke = NA, size = 2) +
  xlab("log(Distance from Boston (km))") + 
  ylab("Bacterial zoonotic pathogen relative abundance (%)") + 
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

rm(list = c("mean_its_soil", "metadata_leaf_16s", "model_zoonotic"))

### FIGURE 2C #################################################################
# Select the soil, top layer data from ITS
metadata_soil_its <- metadata_its[which(metadata_its$sample_type == "Soil, 0-15 cm"),]

# Select the soil, top layer data from 16S
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]

# Select the leaf data from ITS
metadata_leaf_its <- metadata_its[which(metadata_its$sample_type == "Leaf"),]

# Select the leaf data from 16S
metadata_leaf_16s <- metadata_16s[which(metadata_16s$sample_type == "Leaf"),]

# Take the functional guild abundances part of the metadata from ITS and 16S
func_soil_its <- metadata_soil_its[,colnames(metadata_soil_its) %in%
                                             c(fungal_trait_names)]
func_soil_16s <- metadata_soil_16s[,colnames(metadata_soil_16s) %in% 
                                     c(bacterial_trait_names)]
func_leaf_its <- metadata_leaf_its[,colnames(metadata_leaf_its) %in%
                                             c(fungal_trait_names)]
func_leaf_16s <- metadata_leaf_16s[,colnames(metadata_leaf_16s) %in% 
                                     c(bacterial_trait_names)]

# Keep the fungal guilds to be analyzed in this figure
func_soil_its <- as.data.frame(func_soil_its[,which(colnames(func_soil_its) %in% 
                                                      c("perc_animal_pathogen"))])
func_leaf_its <- as.data.frame(func_leaf_its[,which(colnames(func_leaf_its) %in% 
                                                      c("perc_animal_pathogen"))])

# Rename the column names for the figure
colnames(func_soil_its) <- c("Animal fungal pathogens")
colnames(func_soil_its) <- paste0(colnames(func_soil_its), " in soil")
colnames(func_leaf_its) <- c("Animal fungal pathogens")
colnames(func_leaf_its) <- paste0(colnames(func_leaf_its), "in leaves")

# Keep bacterial guilds to be analyzed in this figure
func_soil_16s <- as.data.frame(func_soil_16s[,which(colnames(func_soil_16s) %in% 
                                                      c("perc_zoonotic_pathogen"))])
func_leaf_16s <- as.data.frame(func_leaf_16s[,which(colnames(func_leaf_16s) %in% 
                                                      c("perc_zoonotic_pathogen"))])

# Rename the column names for the figure
colnames(func_soil_16s) <- c("Bacterial zoonotic pathogens")
colnames(func_soil_16s) <- paste0(colnames(func_soil_16s), " in soil")
colnames(func_leaf_16s) <- c("Bacterial hoonotic pathogens")
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
estimate_soil_its <- NULL
for (j in 1:ncol(metadata_soil_heatmap_its)){
  cor_soil_its <- cor.test(func_soil_its[,1],metadata_soil_heatmap_its[,j], method="p")
  estimate_soil_its <- c(estimate_soil_its, cor_soil_its$estimate)
}
Result_soil_its <- rbind(Result_soil_its, estimate_soil_its)
colnames(Result_soil_its) <- colnames(metadata_soil_heatmap_its)
rownames(Result_soil_its) <- "Fungal animal pathogens in soil"

# Leaf ITS
Result_leaf_its <- NULL
estimate_leaf_its <- NULL
for (j in 1:ncol(metadata_leaf_heatmap_its)){
  cor_leaf_its <- cor.test(func_leaf_its[,1],metadata_leaf_heatmap_its[,j], method="p")
  estimate_leaf_its <- c(estimate_leaf_its, cor_leaf_its$estimate)
}
Result_leaf_its <- rbind(Result_leaf_its, estimate_leaf_its)
colnames(Result_leaf_its) <- colnames(metadata_leaf_heatmap_its)
rownames(Result_leaf_its) <- "Fungal animal pathogens in leaves"

# Soil 16S
Result_soil_16s <- NULL
estimate_soil_16s <- NULL
for (j in 1:ncol(metadata_soil_heatmap_16s)){
  cor_soil_16s <- cor.test(func_soil_16s[,1],metadata_soil_heatmap_16s[,j], method="p")
  estimate_soil_16s <- c(estimate_soil_16s, cor_soil_16s$estimate)
}
Result_soil_16s <- rbind(Result_soil_16s, estimate_soil_16s)
colnames(Result_soil_16s) <- colnames(metadata_soil_heatmap_16s)
rownames(Result_soil_16s) <- "Bacterial zoonotic pathogens in soil"

# Leaf 16S
Result_leaf_16s <- NULL
estimate_leaf_16s <- NULL
for (j in 1:ncol(metadata_leaf_heatmap_16s)){
  cor_leaf_16s <- cor.test(func_leaf_16s[,1],metadata_leaf_heatmap_16s[,j], method="p")
  estimate_leaf_16s <- c(estimate_leaf_16s, cor_leaf_16s$estimate)
}
Result_leaf_16s <- rbind(Result_leaf_16s, estimate_leaf_16s)
colnames(Result_leaf_16s) <- colnames(metadata_leaf_heatmap_16s)
rownames(Result_leaf_16s) <- "Bacterial zoonotic pathogens in leaves" 

# Make Pval matrices
# Soil ITS
Result.Pval_soil_its <- NULL
Pval_soil_its <- NULL
for (j in 1:ncol(metadata_soil_heatmap_its)){
  cor_soil_its <- cor.test(func_soil_its[,1],metadata_soil_heatmap_its[,j], method="p")
  Pval_soil_its <- c(Pval_soil_its, cor_soil_its$p.value)
}
Result.Pval_soil_its <- rbind (Result.Pval_soil_its, Pval_soil_its)
colnames(Result.Pval_soil_its) <- colnames(metadata_soil_heatmap_its)
rownames(Result.Pval_soil_its) <- "Fungal animal pathogens in soil"

# Leaf ITS
Result.Pval_leaf_its <- NULL
Pval_leaf_its <- NULL
for (j in 1:ncol(metadata_leaf_heatmap_its)){
  cor_leaf_its <- cor.test(func_leaf_its[,1],metadata_leaf_heatmap_its[,j], method="p")
  Pval_leaf_its <- c(Pval_leaf_its, cor_leaf_its$p.value)
}
Result.Pval_leaf_its <- rbind (Result.Pval_leaf_its, Pval_leaf_its)
colnames(Result.Pval_leaf_its) <- colnames(metadata_leaf_heatmap_its)
rownames(Result.Pval_leaf_its) <- "Fungal animal pathogens in leaves"

# Soil 16S
Result.Pval_soil_16s <- NULL
Pval_soil_16s <- NULL
for (j in 1:ncol(metadata_soil_heatmap_16s)){
  cor_soil_16s <- cor.test(func_soil_16s[,1],metadata_soil_heatmap_16s[,j], method="p")
  Pval_soil_16s <- c(Pval_soil_16s, cor_soil_16s$p.value)
}
Result.Pval_soil_16s <- rbind (Result.Pval_soil_16s, Pval_soil_16s)
colnames(Result.Pval_soil_16s) <- colnames(metadata_soil_heatmap_16s)
rownames(Result.Pval_soil_16s) <- "Bacterial zoonotic pathogens in soil"

# Leaf 16S
Result.Pval_leaf_16s <- NULL
Pval_leaf_16s <- NULL
for (j in 1:ncol(metadata_leaf_heatmap_16s)){
  cor_leaf_16s <- cor.test(func_leaf_16s[,1],metadata_leaf_heatmap_16s[,j], method="p")
  Pval_leaf_16s <- c(Pval_leaf_16s, cor_leaf_16s$p.value)
}
Result.Pval_leaf_16s <- rbind (Result.Pval_leaf_16s, Pval_leaf_16s)
colnames(Result.Pval_leaf_16s) <- colnames(metadata_leaf_heatmap_16s)
rownames(Result.Pval_leaf_16s) <- "Bacterial zoonotic pathogens in leaves"

# Make Asterisk matrices
# Soil ITS
Result.Pval.combined_soil_its <- Result.Pval_soil_its
Result.asterisk_soil_its <- Result.Pval.combined_soil_its
for (k in 1:ncol(Result.Pval.combined_soil_its)){
  for (i in 1:nrow(Result.Pval.combined_soil_its)){
    if(!is.na(Result.Pval.combined_soil_its[i,k])){
      if (Result.Pval.combined_soil_its[i,k] < 0.001/(32*4)){ # multiple comparisons correction -- 4 functional groups * 32 environmental factors
        Result.asterisk_soil_its[i,k] <- "***"
      } else if (Result.Pval.combined_soil_its[i,k] < 0.01/(32*4)){ # multiple comparisons correction -- 4 functional groups * 32 environmental factors
        Result.asterisk_soil_its[i,k] <- "**"
      } else if (Result.Pval.combined_soil_its[i,k] < 0.05/(32*4)){ # multiple comparisons correction -- 4 functional groups * 32 environmental factors
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
      if (Result.Pval.combined_leaf_its[i,k] < 0.001/(32*4)){ # multiple comparisons correction -- 4 functional groups * 32 environmental factors
        Result.asterisk_leaf_its[i,k] <- "***"
      } else if (Result.Pval.combined_leaf_its[i,k] < 0.01/(32*4)){ # multiple comparisons correction -- 4 functional groups * 32 environmental factors
        Result.asterisk_leaf_its[i,k] <- "**"
      } else if (Result.Pval.combined_leaf_its[i,k] < 0.05/(32*4)){ # multiple comparisons correction -- 4 functional groups * 32 environmental factors
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
      if (Result.Pval.combined_soil_16s[i,k] < 0.001/(32*4)){ # multiple comparisons correction -- 4 functional groups * 32 environmental factors
        Result.asterisk_soil_16s[i,k] <- "***"
      } else if (Result.Pval.combined_soil_16s[i,k] < 0.01/(32*4)){ # multiple comparisons correction -- 4 functional groups * 32 environmental factors
        Result.asterisk_soil_16s[i,k] <- "**"
      } else if (Result.Pval.combined_soil_16s[i,k] < 0.05/(32*4)){ # multiple comparisons correction -- 4 functional groups * 32 environmental factors
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
      if (Result.Pval.combined_leaf_16s[i,k] < 0.001/(32*4)){ # multiple comparisons correction -- 4 functional groups * 32 environmental factors
        Result.asterisk_leaf_16s[i,k] <- "***"
      } else if (Result.Pval.combined_leaf_16s[i,k] < 0.01/(32*4)){ # multiple comparisons correction -- 4 functional groups * 32 environmental factors
        Result.asterisk_leaf_16s[i,k] <- "**"
      } else if (Result.Pval.combined_leaf_16s[i,k] < 0.05/(32*4)){ # multiple comparisons correction -- 4 functional groups * 32 environmental factors
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

melt_all <- rbind(melt_leaf_its, melt_soil_its, melt_leaf_16s, melt_soil_16s)
melt_asterisk_all <- rbind(melt.asterisk_leaf_its, melt.asterisk_soil_its, 
                           melt.asterisk_leaf_16s, melt.asterisk_soil_16s)

# Combine the correlation values with the p-values
vdata_its <- cbind(melt_all, melt_asterisk_all$value)

# Rename the columns
colnames(vdata_its) <- c("microbe", "environment", "coef", "Pval")

# Make an order of the environmental factors based on their correlation with the first microbial group
order_of_env_factors <- vdata_its[which(vdata_its$microbe == "Fungal animal pathogens in leaves"),]
vdata_its$environment <- factor(vdata_its$environment, 
                                levels = order_of_env_factors$environment[order(order_of_env_factors$coef)])

# plot figure 2c
fig2c <- ggplot(vdata_its) +
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

fig2c
ggsave(paste0("atherton_figure2c_functionalguilds_vs_environment_heatmap", 
              format(Sys.Date(), "%Y%m%d"),".pdf"), plot = fig2c, 
       device = cairo_pdf, width = 5, height = 10, units = "in", dpi = 300)