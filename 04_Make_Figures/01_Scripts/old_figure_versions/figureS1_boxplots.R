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

# read in sequence data
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/02_Clean_Data/05_Transform_Data/ITS/Rarefied_ASV_Tables/")
o_rare_its <- read.csv("atherton_ITS_osoil_ASV_table_rare_taxonomy_guild_ASVID_20250130.csv",
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

### FIGURE S1A #################################################################
# Select the root data from ITS
metadata_root_its <- metadata_its[which(metadata_its$sample_type %in% 
                                          c("Roots, 0-15 cm", 
                                            "Roots, 15-30 cm")),]

# For figure caption: number of trees' data used (n = 52)
print("Number of trees used in the root analysis:")
length(unique(metadata_root_its$tree_id))

# Select the root data from 16S
metadata_root_16s <- metadata_16s[which(metadata_16s$sample_type %in% 
                                          c("Roots, 0-15 cm", 
                                            "Roots, 15-30 cm")),]

# For figure caption: number of trees' data used (n = 56)
print("Number of trees used in the root analysis:")
length(unique(metadata_root_16s$tree_id))

# Create the x-axis categories
metadata_root_its$site_edge <- paste0(metadata_root_its$tree_urban_type, " ", 
                                      metadata_root_its$tree_edge_type)
metadata_root_its$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                                    "Urban Street Tree",
                                    metadata_root_its$site_edge)
metadata_root_its$site_edge <- factor(metadata_root_its$site_edge, 
                                      levels = c("Urban Street Tree", 
                                                 "Urban Forest Edge",
                                                 "Urban Forest Interior",
                                                 "Rural Forest Edge", 
                                                 "Rural Forest Interior"))

metadata_root_16s$site_edge <- paste0(metadata_root_16s$tree_urban_type, " ", 
                                      metadata_root_16s$tree_edge_type)
metadata_root_16s$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                                    "Urban Street Tree",
                                    metadata_root_16s$site_edge)
metadata_root_16s$site_edge <- factor(metadata_root_16s$site_edge, 
                                      levels = c("Urban Street Tree", 
                                                 "Urban Forest Edge",
                                                 "Urban Forest Interior",
                                                 "Rural Forest Edge", 
                                                 "Rural Forest Interior"))

# Take the functional guild abundances part of the metadata from roots
metadata_sap <- metadata_root_its[which(colnames(metadata_root_its) %in% 
                                               c("sample_name", "tree_id", 
                                                 "plot_name", "tree_edge_type", 
                                                 "tree_urban_type", 
                                                 "perc_wood_saprotroph", 
                                                 "site_edge"))]
metadata_root_cell <- metadata_root_16s[which(colnames(metadata_root_16s) %in% 
                                                c("sample_name", "tree_id", 
                                                  "plot_name", "tree_edge_type", 
                                                  "tree_urban_type", 
                                                  "perc_cellulolytic", 
                                                  "site_edge"))]
metadata_root_lig <- metadata_root_16s[which(colnames(metadata_root_16s) %in% 
                                               c("sample_name", "tree_id", 
                                                 "plot_name", "tree_edge_type", 
                                                 "tree_urban_type", 
                                                 "perc_lignolytic", 
                                                 "site_edge"))]

# Rename the columns
colnames(metadata_sap) <- c("sample_name", "tree_id", "plot_name", 
                                 "tree_edge_type", "tree_urban_type", 
                                 "rel_abund", "site_edge")
colnames(metadata_root_cell) <- c("sample_name", "tree_id", "plot_name", 
                                 "tree_edge_type", "tree_urban_type", 
                                 "rel_abund", "site_edge")
colnames(metadata_root_lig) <- c("sample_name", "tree_id", "plot_name", 
                                  "tree_edge_type", "tree_urban_type", 
                                  "rel_abund", "site_edge")

# Add the name of the functional guild
metadata_sap$guild <- "Fungal wood decomposers"
metadata_root_cell$guild <- "Cellulolytic bacteria"
metadata_root_lig$guild <- "Lignolytic bacteria"

# Combine the soil and leaf datasets for plotting
vdata <- rbind(metadata_sap, metadata_root_cell, metadata_root_lig)

# Order the guild as a factor for the figure
vdata$guild <- factor(vdata$guild, levels = c("Cellulolytic bacteria", 
                                              "Lignolytic bacteria",
                                              "Fungal wood decomposers"))

### STATISTICS WOOD DECOMPOSERS ###
model_wood <- lme(formula(rel_abund ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_sap)
summary(model_wood)
summary(emmeans(model_wood, pairwise ~ site_edge))

### STATISTICS CELLULOLYTIC BACTERIA ###
model_cell <- lme(formula(rel_abund ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_root_cell)
summary(model_cell)
summary(emmeans(model_cell, pairwise ~ site_edge))

### STATISTICS LIGNOLYTIC BACTERIA ###
model_lig <- lme(formula(rel_abund ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_root_lig)
summary(model_lig)
summary(emmeans(model_lig, pairwise ~ site_edge))

# Make the text for adding Tukey groups and record the location of the group id on the y-axis
vdata$tukey <- NA
vdata$location <- NA

for(i in 1:nrow(vdata)){
  guild <- as.character(vdata$guild[i])
  site_edge <- as.character(vdata$site_edge[i])
  
  if(guild == "Fungal wood decomposers"){
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "ab"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "ab"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "ab"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
    }
  } else if(guild == "Cellulolytic bacteria"){
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
    }
  } else if(guild == "Lignolytic bacteria"){
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.25
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

# Plot Figure 1a
figs1a <- ggplot() +
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
                              group = guild), alpha = 0.75, stroke = NA,
              position = position_jitterdodge(jitter.width = 0.5, 
                                              dodge.width = 0.75)) +
  # this adds the Tukey groups labels defined above
  geom_text(data=vdata, family = "Arial",
            aes(x = site_edge, y = location, color = guild, group = guild, 
                label = tukey), position = position_dodge(width = 0.75),
            size = 3, parse = TRUE) +
  theme_classic() + labs (y="log(Functional group relative abundance (%)", 
                          x="", color = "", fill = "") + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), 
        legend.position = c(0.6, 1), 
        legend.background = element_rect(fill = NA),
        legend.direction = "horizontal") +
  scale_color_manual(values = categorical_pal[c(1,3,5)]) + 
  scale_fill_manual(values = categorical_pal[c(1,3,5)]) +
  guides(color = "none") +
  ggtitle("Wood decomposers on tree roots")

figs1a
ggsave(paste0("atherton_figures1a_wooddecomposersroots_vs_treetype_boxplot_",
              format(Sys.Date(), "%Y%m%d"),".pdf"), 
       plot = figs1a, dpi = 300, width = 6, height = 3.5, units = "in", 
       device = cairo_pdf)

rm(list = c("metadata_root_16s", "metadata_root_cell", "metadata_root_its", 
            "metadata_root_lig", "metadata_sap", "model_cell", "model_lig", 
            "model_wood", "vdata"))

### FIGURE S1B ##################################################################
# Select belowground data from ITS
metadata_soil_its <- metadata_its[which(metadata_its$sample_type != "Leaf"),]

### STATISTICS ECM VS WOOD DECOMPOSERS ###
model_ecm_wood <- lme(formula(perc_ecm ~ perc_wood_saprotroph), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_its)
anova(model_ecm_wood)
r.squaredGLMM(model_ecm_wood)

figs1b <- ggplot(metadata_soil_its, aes(x = perc_ecm, 
                                        y = perc_wood_saprotroph)) + 
  geom_point(alpha = 0.5, col = "#003049", aes(shape = sample_type)) +
  xlab("ECM relative abundance (%)))") + 
  ylab("log(Wood decomposer relative abundance (%))") +
  labs(shape = "") +
  geom_smooth(method = "lm", alpha = 0.12, col = "#003049", fill = "#003049", 
              aes(x = perc_ecm,
                  y = perc_wood_saprotroph,
                  col = Guild, fill = Guild,
                  ymax = after_stat(y + se * sqrt(length(y))),
                  ymin = after_stat(y - se * sqrt(length(y))))) +
  theme_classic() + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10),
        legend.position = "bottom") + 
  scale_shape_manual(values = shapes[2:5]) +
  ggtitle("Fungal wood decomposer and ECM relative abundance belowground")

figs1b
ggsave(paste0("atherton_figures1b_wooddecomposersroots_vs_ecm_",
              format(Sys.Date(), "%Y%m%d"),".pdf"), 
       plot = figs1b, dpi = 300, width = 6, height = 3.5, units = "in", 
       device = cairo_pdf)

rm(list = c("metadata_soil_its", "model_ecm_wood"))

### FIGURE S1C #################################################################
# Subset the wood saprotrophs from the reads in the rarefied O Soil ITS data
o_sap_paths <- subset(o_rare_its, primary_lifestyle %in% c("wood_saprotroph"))

# Subset the Plant pathogens from the wood saprotrophs
o_sap_wood_paths <- subset(o_sap_paths, Plant_pathogenic_capacity_template %in% 
                             c("leaf/fruit/seed_pathogen", "root_pathogen", 
                               "wood_pathogen"))

# Subset the human pathogens from the wood saprotrophs
o_sap_human_paths <- subset(o_sap_paths, Animal_biotrophic_capacity_template 
                            %in% c("opportunistic_human_parasite"))

# Find the total sum of the plant pathogen / wood saprotrophs
plant_path_abund <- colSums(o_sap_wood_paths[12:253])
plant_path_abund <- data.frame(plant_path_abund)

# Find the total sum of the human pathogen / wood saprotrophs
human_path_abund <- colSums(o_sap_human_paths[12:253])
human_path_abund <- data.frame(human_path_abund)

# Make the human and plant pathogen data one dataframe
microbe.data<-merge(human_path_abund, plant_path_abund, by = 0)

# Format the metadata
metadata_its$sample_name <- gsub(" ","", metadata_its$sample_name)

# Pull the relevant metadata for the samples for plotting
sample_type<-metadata_its$sample_type[which(metadata_its$sample_name %in% 
                                              microbe.data$Row.names)]
soil_horizon<-metadata_its$soil_horizon[which(metadata_its$sample_name %in% 
                                                microbe.data$Row.names)]
DFB<-metadata_its$tree_dist_to_boston_km[which(metadata_its$sample_name %in% 
                                                 microbe.data$Row.names)]
DFE<-metadata_its$tree_distance_from_edge[which(metadata_its$sample_name %in% 
                                                  microbe.data$Row.names)]
plot_name<-metadata_its$plot_name[which(metadata_its$sample_name %in% 
                                          microbe.data$Row.names)]
site_name<-metadata_its$site_name[which(metadata_its$sample_name %in% 
                                          microbe.data$Row.names)]
tree_id<-metadata_its$tree_id[which(metadata_its$sample_name %in% 
                                      microbe.data$Row.names)]
tree_edge_type<-metadata_its$tree_edge_type[which(metadata_its$sample_name %in%
                                                    microbe.data$Row.names)]
tree_urban_type<-metadata_its$tree_urban_type[which(metadata_its$sample_name %in% 
                                                      microbe.data$Row.names)]
new_tree_type<-paste(tree_urban_type, tree_edge_type, sep = " ")
new_tree_type<-gsub("Urban Street Tree Urban Street Tree", "Urban Street Tree", 
                    new_tree_type)
wood_paths<-microbe.data$plant_path_abund[which(microbe.data$Row.names %in% 
                                                  metadata_its$sample_name)]
human_paths<-microbe.data$human_path_abund[which(microbe.data$Row.names %in% 
                                                   metadata_its$sample_name)]

wood_paths_df <- as.data.frame(cbind(wood_paths, new_tree_type))
human_paths_df <- as.data.frame(cbind(human_paths, new_tree_type))
wood_paths_df$wood_paths <- as.numeric(wood_paths)
human_paths_df$human_paths <- as.numeric(human_paths)
colnames(wood_paths_df) <- c("rel_abund", "tree_type")
colnames(human_paths_df) <- c("rel_abund", "tree_type")
wood_paths_df$pathogen_type <- "Potential plant pathogen"
human_paths_df$pathogen_type <- "Potential human pathogen"
pathogens_df <- rbind(wood_paths_df, human_paths_df)

# Find the average abundance and standard error of the human and wood pathogens
mean_root_hp <- aggregate(log(human_paths+1)~new_tree_type, FUN= "mean", 
                          na.action = na.omit)
se_root_hp <- aggregate(log(human_paths+1)~new_tree_type, FUN= standard_error,
                        na.action = na.omit)
mean_root_pp <- aggregate(log(wood_paths+1)~new_tree_type, FUN= "mean",
                          na.action = na.omit)
se_root_pp <- aggregate(log(wood_paths+1)~new_tree_type, FUN= standard_error,
                        na.action = na.omit)

# Merge the means and standard errors together
root_hp <- cbind(mean_root_hp, se_root_hp[,ncol(se_root_hp)],
                 "Potential human pathogen")
root_pp <- cbind(mean_root_pp, se_root_pp[,ncol(se_root_pp)],
                 "Potential plant pathogen")

# Rename the columns
colnames(root_hp) <- c("site_edge","Abundance","SE","Microbe")
colnames(root_pp) <- c("site_edge","Abundance","SE","Microbe")

# Make one dataframe for plotting
vdata <- rbind(root_hp, root_pp)
vdata$max <- vdata$Abundance + vdata$SE

### STATISTICS HUMAN PATHOGENS ###
model_human <- lme(formula(log(human_paths+1) ~ new_tree_type), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML")
summary(model_human)
summary(emmeans(model_human, pairwise ~ new_tree_type))

### STATISTICS WOOD PATHOGENS ###
model_plant <- lme(formula(log(wood_paths+1) ~ new_tree_type), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML")
summary(model_plant)
summary(emmeans(model_plant, pairwise ~ new_tree_type))

# Define the Tukey groups
vdata$tukey <- c("", "", "", "", "", "", "", "", "", "")

# Format the x-axis tickmarks
vdata$site_edge <- gsub(" ", "\n", vdata$site_edge)
vdata$site_edge <- factor(vdata$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

pathogens_df$tree_type <- gsub(" ", "\n", pathogens_df$tree_type)
pathogens_df$tree_type <- factor(pathogens_df$tree_type, 
                                     levels = c("Urban\nStreet\nTree",
                                                "Urban\nForest\nEdge",
                                                "Urban\nForest\nInterior",
                                                "Rural\nForest\nEdge",
                                                "Rural\nForest\nInterior"))

# Order the functional guilds as a factor for plotting
vdata$Microbe <- factor(vdata$Microbe, levels = c("Potential plant pathogen", 
                                                  "Potential human pathogen"))

figs1c <- ggplot() +
  stat_summary(data = pathogens_df, aes(x = tree_type, y = log(rel_abund+1), 
                                        group = pathogen_type, 
                                        color = pathogen_type),
               fun = median,
               geom = "line",
               position = position_dodge(width = 0.75),
               linetype = "dashed",
               linewidth = 0.8) +
  # this creates the points that show the average abundance for each guild and tree type
  geom_boxplot(data = pathogens_df, outlier.shape = NA,
               aes(x=tree_type, y=log(rel_abund+1), fill=pathogen_type), 
               alpha = 0.5) + 
  geom_jitter(data=pathogens_df, aes(x = tree_type, y = log(rel_abund+1), 
                                     color = pathogen_type, 
                              group = pathogen_type), alpha = 0.75, stroke = NA,
              position = position_jitterdodge(jitter.width = 0.5, 
                                              dodge.width = 0.75)) +
  theme_classic() + labs (y="log(Potential pathogens relative abundance (%))", 
                          x="", fill = "") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_markdown(size = 8), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10),
        legend.position = c(0.75, 1),
        legend.direction = "horizontal") +
  scale_color_manual(values = linear_pal[c(2,1)]) + 
  scale_fill_manual(values = linear_pal[c(2,1)]) +
  guides(color = "none", shape = guide_legend(title = ""))+
  ggtitle("Wood decomposer potential human and plant pathogen relative abundance in soils")

figs1c
ggsave(paste0("atherton_figures1c_pathogenicwooddecomposersoils_vs_treetype_boxplot_",
              format(Sys.Date(), "%Y%m%d"),".pdf"), 
       plot = figs1c, dpi = 300, width = 6, height = 3.5, units = "in", 
       device = cairo_pdf)
