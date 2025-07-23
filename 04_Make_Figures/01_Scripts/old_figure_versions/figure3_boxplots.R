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
metadata_its$perc_dung_saprotroph <- log(metadata_its$perc_dung_saprotroph+1)
metadata_its$perc_wood_saprotroph <- log(metadata_its$perc_wood_saprotroph+1)
metadata_its$perc_litter_saprotroph <- log(metadata_its$perc_litter_saprotroph+1)
metadata_its$perc_soil_saprotroph <- log(metadata_its$perc_soil_saprotroph+1)

metadata_16s$perc_chitinolytic <- log(metadata_16s$perc_chitinolytic+1)
metadata_16s$perc_partial_nitrification <- log(metadata_16s$perc_partial_nitrification+1)
metadata_16s$perc_dissim_nitrate_reduction <- log(metadata_16s$perc_dissim_nitrate_reduction+1)
metadata_16s$perc_denitrification <- log(metadata_16s$perc_denitrification+1)
metadata_16s$perc_methanotroph <- log(metadata_16s$perc_methanotroph+1)

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
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 88)
print("Number of trees used in the soil analysis:")
length(unique(metadata_soil_16s$tree_id))

### STATISTICS SOIL NO3- ###
model_no3 <- lme(formula(soil_no3 ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_16s)
summary(model_no3)
summary(emmeans(model_no3, pairwise ~ site_edge))

# Create the Tukey groups
for(i in 1:nrow(metadata_soil_16s)){
  site_edge <- as.character(metadata_soil_16s$site_edge[i])
    if(site_edge == "Urban Street Tree"){
      metadata_soil_16s$tukey[i] <- "a"
      metadata_soil_16s$location[i] <- max(metadata_soil_16s$soil_no3[which(metadata_soil_16s$site_edge == site_edge)]) + 0.05
    } else if(site_edge == "Urban Forest Edge"){
      metadata_soil_16s$tukey[i] <- "ab"
      metadata_soil_16s$location[i] <- max(metadata_soil_16s$soil_no3[which(metadata_soil_16s$site_edge == site_edge)]) + 0.05
    } else if(site_edge == "Urban Forest Interior"){
      metadata_soil_16s$tukey[i] <- "b"
      metadata_soil_16s$location[i] <- max(metadata_soil_16s$soil_no3[which(metadata_soil_16s$site_edge == site_edge)]) + 0.05
    } else if(site_edge == "Rural Forest Edge"){
      metadata_soil_16s$tukey[i] <- "ab"
      metadata_soil_16s$location[i] <- max(metadata_soil_16s$soil_no3[which(metadata_soil_16s$site_edge == site_edge)]) + 0.05
    } else if(site_edge == "Rural Forest Interior"){
      metadata_soil_16s$tukey[i] <- "b"
      metadata_soil_16s$location[i] <- max(metadata_soil_16s$soil_no3[which(metadata_soil_16s$site_edge == site_edge)]) + 0.05
    }
}

# Reformat the x-axis tickmark labels
metadata_soil_16s$site_edge <- gsub(" ", "\n", metadata_soil_16s$site_edge)
metadata_soil_16s$site_edge <- factor(metadata_soil_16s$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

# Visualize
fig3a <- ggplot() +
  stat_summary(data = metadata_soil_16s, aes(x = site_edge, y = soil_no3, 
                                             group = 1),
               fun = median,
               geom = "line",
               linetype = "dashed",
               linewidth = 0.8) +
  geom_boxplot(data = metadata_soil_16s, outlier.shape = NA,
               aes(x=site_edge, y=soil_no3, fill=site_edge), alpha = 0.5) + 
  geom_jitter(data=metadata_soil_16s, aes(x = site_edge, y = soil_no3, 
                                          color = site_edge, alpha = 0.75, 
                                          stroke = NA)) +
  # this adds the Tukey groups labels defined above
  geom_text(data = metadata_soil_16s, direction = "y",
            aes(x=site_edge, y=location, label = tukey), 
            size = 3) +
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
  scale_fill_manual(values = categorical_pal) + 
  scale_color_manual(values = categorical_pal) + 
  guides(color = "none", shape = guide_legend(title = ""))+
  ggtitle("Available soil NO<sub>3</sub><sup>-</sup>")

fig3a
ggsave(paste0("atherton_figure3a_no3soil_vs_treetype_boxplot_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig3a, dpi = 300, width = 6, height = 3.5, units = "in")

rm(list = c("metadata_soil_16s", "model_no3"))

### FIGURE 1B ##################################################################
# Select the soil, top layer data from 16S
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 88)
print("Number of trees used in the soil analysis:")
length(unique(metadata_soil_16s$tree_id))

# Create the x-axis tick marks
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

# Take the functional guild abundances part of the metadata from soil
metadata_chitin <- metadata_soil_16s[which(colnames(metadata_soil_16s) %in% 
                                               c("sample_name", "tree_id", 
                                                 "plot_name", "tree_edge_type", 
                                                 "tree_urban_type", 
                                                 "perc_chitinolytic", 
                                                 "site_edge"))]
metadata_nit <- metadata_soil_16s[which(colnames(metadata_soil_16s) %in%
                                               c("sample_name", "tree_id", 
                                                 "plot_name", "tree_edge_type", 
                                                 "tree_urban_type", 
                                                 "perc_partial_nitrification", 
                                                 "site_edge"))]
metadata_dnr<- metadata_soil_16s[which(colnames(metadata_soil_16s) %in%
                                               c("sample_name", "tree_id", 
                                                 "plot_name", "tree_edge_type", 
                                                 "tree_urban_type", 
                                                 "perc_dissim_nitrate_reduction", 
                                                 "site_edge"))]
metadata_denit <- metadata_soil_16s[which(colnames(metadata_soil_16s) %in%
                                          c("sample_name", "tree_id", 
                                            "plot_name", "tree_edge_type", 
                                            "tree_urban_type", 
                                            "perc_denitrification", 
                                            "site_edge"))]

# Rename the columns
colnames(metadata_chitin) <- c("sample_name", "tree_id", "plot_name", 
                                 "tree_edge_type", "tree_urban_type", 
                                 "rel_abund", "site_edge")
colnames(metadata_nit) <- c("sample_name", "tree_id", "plot_name",
                                 "tree_edge_type", "tree_urban_type", 
                                 "rel_abund", "site_edge")
colnames(metadata_dnr) <- c("sample_name", "tree_id", "plot_name", 
                                 "tree_edge_type", "tree_urban_type", 
                                 "rel_abund", "site_edge")
colnames(metadata_denit) <- c("sample_name", "tree_id", "plot_name", 
                            "tree_edge_type", "tree_urban_type", 
                            "rel_abund", "site_edge")

# Add the name of the functional guild
metadata_chitin$guild <- "Chitinolytic bacteria"
metadata_nit$guild <- "Nitrifying bacteria"
metadata_dnr$guild <- "Dissimilatory nitrate reducing bacteria"
metadata_denit$guild <- "Denitrifying bacteria"

# Combine the soil and leaf datasets for plotting
vdata <- rbind(metadata_chitin, metadata_nit, metadata_denit, metadata_denit)

# Order the microbial groups for the plot
vdata$guild <- factor(vdata$guild, levels = c("Chitinolytic bacteria", 
                                              "Nitrifying bacteria",
                                              "Dissimilatory nitrate reducing bacteria",
                                              "Denitrifying bacteria"))

### STATISTICS CHITINOLYTIC BACTERIA ###
model_chitin <- lme(formula(perc_chitinolytic ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_16s)
summary(model_chitin)
summary(emmeans(model_chitin, pairwise ~ site_edge))

### STATISTICS DISSIMILATORY NITRATE REDUCING BACTERIA ###
model_dnr <- lme(formula(perc_dissim_nitrate_reduction ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_16s)
summary(model_dnr)
summary(emmeans(model_dnr, pairwise ~ site_edge))

### STATISTICS NITRIFYING BACTERIA ###
model_nit <- lme(formula(perc_partial_nitrification ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_16s)
summary(model_nit)
summary(emmeans(model_nit, pairwise ~ site_edge))

### STATISTICS DENITRIFYING BACTERIA ###
model_denit <- lme(formula(perc_denitrification ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_16s)
summary(model_denit)
summary(emmeans(model_denit, pairwise ~ site_edge))

# Make the text for adding Tukey groups and record the location of the group id on the y-axis
vdata$tukey <- NA
vdata$location <- NA

for(i in 1:nrow(vdata)){
  guild <- as.character(vdata$guild[i])
  site_edge <- as.character(vdata$site_edge[i])
  
  if(guild == "Chitinolytic bacteria"){
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    }
  } else if(guild == "Nitrifying bacteria"){
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    }
  } else if(guild == "Dissimilatory nitrate reducing bacteria"){
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    }
  } else if(guild == "Denitrifying bacteria"){
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    }
  }
}

# Reformat the x-axis tick marks for the plot
vdata$site_edge <- gsub(" ", "\n", vdata$site_edge)
vdata$site_edge <- factor(vdata$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

# Summarize the data for adding the error bars, mean, and standard error
summary_data <- vdata %>%
  group_by(site_edge, guild) %>%
  summarise(
    mean = mean(rel_abund, na.rm = TRUE),
    se = sd(rel_abund, na.rm = TRUE) / sqrt(n()),
    sd = sd(rel_abund, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot Figure 1a
fig3b <- ggplot(vdata, aes(x = site_edge, color = guild)) +
  # this adds the jittered points
  geom_jitter(aes(y = rel_abund, group = guild, stroke = NA), 
              size = 0.75,  alpha = 0.2, 
              position = position_jitterdodge(jitter.width = 0.5, 
                                              dodge.width = 0.75)) +
  # this adds the lines that follow the mean across tree types
  stat_summary(data = vdata,
               aes(y = rel_abund, group = guild, color = guild),
               fun = mean, geom = "line",
               position = position_dodge(width = 0.75), linetype = "dashed",
               linewidth = 0.8, alpha = 0.5) +
  # this adds the error bars
  geom_errorbar(data = summary_data,
                aes(x = site_edge, ymin = mean - sd, ymax = mean + sd,
                    group = guild, color = guild),
                position = position_dodge(width = 0.75), width = 0.1,
                alpha = 0.5) +
  # this adds the box with the mean and standard error
  geom_crossbar(data = summary_data,
                aes(x = site_edge, y = mean, ymin = mean - se, ymax = mean + se,
                    fill = guild),
                position = position_dodge(width = 0.75), width = 0.5,
                alpha = 0.5, fatten = 1) +
  # this adds the tukey groups
  geom_text(data = vdata, family = "Arial",
            aes(x = site_edge, y = location, color = guild, group = guild,
                label = tukey),
            position = position_dodge(width = 0.75), size = 3, parse = TRUE,
            show.legend = FALSE) +
  theme_classic() +
  labs(y = "ECM relative abundance (%)", x = "", color = "", fill = "") +
  scale_y_continuous(
    sec.axis = sec_axis(~./20, name = "log(Epiphyte relative abundance (%))")) +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color = "black"), 
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 0),
        legend.key.size = unit(0.15, "in"),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(size = 10),
        legend.position = c(0.8, 1.05)) +
  scale_color_manual(values = categorical_pal[c(1,2,3,4)]) +
  scale_fill_manual(values = categorical_pal[c(1,2,3,4)]) +
  ggtitle("Nitrogen cycling bacteria in soil")

fig3b
ggsave(paste0("atherton_figure3b_ncyclingsoil_vs_treetype_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig3b, dpi = 300, width = 5.75, height = 3.5, units = "in")

rm(list = c("metadata_chitin", "metadata_denit", "metadata_dnr", "metadata_nit",
            "metadata_soil_16s", "model_chitin", "model_denit", "model_dnr",
            "model_nit", "vdata"))

### FIGURE 3C #################################################################
# Select the soil, top layer from ITS
metadata_soil_its <- metadata_its[which(metadata_its$sample_type == "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 87)
print("Number of trees used in the soil analysis:")
length(unique(metadata_soil_its$tree_id))

# Create the x-axis tick marks
metadata_soil_its$site_edge <- paste0(metadata_soil_its$tree_urban_type, " ", 
                                      metadata_soil_its$tree_edge_type)
metadata_soil_its$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                                    "Urban Street Tree",
                                    metadata_soil_its$site_edge)
metadata_soil_its$site_edge <- factor(metadata_soil_its$site_edge, 
                                      levels = c("Urban Street Tree", 
                                                 "Urban Forest Edge", 
                                                 "Urban Forest Interior", 
                                                 "Rural Forest Edge", 
                                                 "Rural Forest Interior"))

# Select the abundances of the saprotrophs
data <- cbind(metadata_soil_its$perc_dung_saprotroph, 
              metadata_soil_its$perc_wood_saprotroph, 
              metadata_soil_its$perc_litter_saprotroph, 
              metadata_soil_its$perc_soil_saprotroph)

# Rename the columns
colnames(data) <- c("Dung","Wood","Litter","Soil")

# Calculate the mean abundance
mean <- aggregate(data~site_edge, data = metadata_soil_its, FUN= "mean", 
                  na.action = na.omit)

# Reformat the dataframe
vdata <- melt(mean)

# Define the color pallette
barchart_pal <- c("#B3CDDE", "#669BBC", "#003049", "#001825")

# Order the decomposers in decreasing abundance
vdata <- vdata[order(vdata$variable, decreasing = TRUE),]

### STATISTICS TOTAL DECOMPOSER ABUNDANCE ###
metadata_soil_its$sap_tot <- metadata_soil_its$perc_dung_saprotroph + 
  metadata_soil_its$perc_wood_saprotroph + 
  metadata_soil_its$perc_litter_saprotroph + 
  metadata_soil_its$perc_soil_saprotroph

model_total <- lme(formula(sap_tot ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_its)
summary(model_total)
summary(emmeans(model_total, pairwise ~ site_edge))

### STATISTICS DUNG DECOMPOSERS ###
model_dung <- lme(formula(perc_dung_saprotroph ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_its)
summary(model_dung)
summary(emmeans(model_dung, pairwise ~ site_edge))

### STATISTICS WOOD DECOMPOSERS ###
model_wood <- lme(formula(perc_wood_saprotroph ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_its)
summary(model_wood)
summary(emmeans(model_wood, pairwise ~ site_edge))

### STATISTICS LITTER DECOMPOSERS ###
model_litter <- lme(formula(perc_litter_saprotroph ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_its)
summary(model_litter)
summary(emmeans(model_litter, pairwise ~ site_edge))

### STATISTICS SOIL DECOMPOSERS ###
model_soil <- lme(formula(perc_soil_saprotroph ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_its)
summary(model_soil)
summary(emmeans(model_soil, pairwise ~ site_edge))

# Format the data and find the point where the Tukey groups should be placed
vdata <- vdata %>%
  group_by(site_edge) %>%
  mutate(cumulative = cumsum(value),            # Cumulative sum
         midpoint = cumulative - value / 2)     # Midpoint of each segment
vdata$midpoint[17:20] <- vdata$cumulative[17:20] + 0.4
vdata_cumsum <- vdata %>%
  group_by(site_edge) %>%
  summarise(cumulative_sum = sum(value))

# Make the Tukey groups text labels
vdata$tukey <- c("b", "a", "a", "a", "a", rep(c("a", "b", "b", "b", "b"), 3))
vdata$color <- c(rep("white",10), rep("black",10))

# Format the x-axis labels
vdata$site_edge <- gsub(" ", "\n", vdata$site_edge)
vdata$site_edge <- factor(vdata$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

# Plot Figure 3c
fig3c <- ggplot(vdata,  mapping=aes(x=site_edge, y=value, fill=variable))+
  # this makes the stacked barchart
  geom_bar(aes(), stat="identity", position="stack",color="black")+
  # this makes the Tukey group labels
  geom_text(aes(y = midpoint, label = tukey), color = vdata$color, family = "Arial") +
  # this makes the between-group whole bar asterisks
  theme_classic()+ geom_signif(
    comparisons = list(c("Urban\nStreet\nTree", "Urban\nForest\nEdge"),
                       c("Urban\nStreet\nTree", "Urban\nForest\nInterior"),
                       c("Urban\nStreet\nTree", "Rural\nForest\nEdge"),
                       c("Urban\nStreet\nTree", "Rural\nForest\nInterior")),  # Specify pairwise comparisons
    y_position = c(max(vdata_cumsum$cumulative_sum) + 0.1,
                   max(vdata_cumsum$cumulative_sum) + 0.7,
                   max(vdata_cumsum$cumulative_sum) + 1.3,
                   max(vdata_cumsum$cumulative_sum) + 1.9),  # Position above the bars
    annotations = c("**","***", "**", "***"), textsize = 2.5) +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_markdown(size = 8), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), 
        legend.position = c(0.7,1), 
        legend.direction="horizontal") + 
  scale_fill_manual(values = barchart_pal) + 
  labs (x="",y="log(Fungal decomposer relative abundance (%))", fill="") + 
  ggtitle("Fungal decomposers in soil")

fig3c
ggsave(paste0("atherton_figure3c_decomposerssoil_vs_treetype_", 
              format(Sys.Date(), "%Y%m%d"),".pdf"), device = cairo_pdf, 
       plot = fig3c, dpi = 300, width = 5.75, height = 3.5, units = "in")

rm(list = c("data", "mean", "metadata_soil_its", "vdata", "vdata_cumsum"))

### FIGURE 3D #################################################################
# Select leaf data from 16S
metadata_leaf_16s <- metadata_16s[which(metadata_16s$sample_type == "Leaf"),]

# For figure caption: number of trees' data used (n = 62)
print("Number of trees used in the soil analysis:")
length(unique(metadata_leaf_16s$tree_id))

# Format the x-axis tick marks
metadata_leaf_16s$site_edge <- paste0(metadata_leaf_16s$tree_urban_type, " ", 
                                      metadata_leaf_16s$tree_edge_type)
metadata_leaf_16s$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                                    "Urban Street Tree",
                                    metadata_leaf_16s$site_edge)
metadata_leaf_16s$site_edge <- factor(metadata_leaf_16s$site_edge, 
                                      levels = c("Urban Street Tree", 
                                                 "Urban Forest Edge", 
                                                 "Urban Forest Interior",
                                                 "Rural Forest Edge", 
                                                 "Rural Forest Interior"))

### STATISTICS METHANOTROPHIC BACTERIA ###
model_meth <- lme(formula(perc_methanotroph ~ site_edge), 
                 random = ~1|plot_name/tree_id, na.action = na.exclude, 
                 method = "ML", data = metadata_leaf_16s)
summary(model_meth)
summary(emmeans(model_meth, pairwise ~ site_edge))

# Create the Tukey groups
for(i in 1:nrow(metadata_leaf_16s)){
  site_edge <- as.character(metadata_leaf_16s$site_edge[i])
  if(site_edge == "Urban Street Tree"){
    metadata_leaf_16s$tukey[i] <- "b"
    metadata_leaf_16s$location[i] <- max(metadata_leaf_16s$perc_methanotroph[which(metadata_leaf_16s$site_edge == site_edge)]) + 0.05
  } else if(site_edge == "Urban Forest Edge"){
    metadata_leaf_16s$tukey[i] <- "ab"
    metadata_leaf_16s$location[i] <- max(metadata_leaf_16s$perc_methanotroph[which(metadata_leaf_16s$site_edge == site_edge)]) + 0.05
  } else if(site_edge == "Urban Forest Interior"){
    metadata_leaf_16s$tukey[i] <- "a"
    metadata_leaf_16s$location[i] <- max(metadata_leaf_16s$perc_methanotroph[which(metadata_leaf_16s$site_edge == site_edge)]) + 0.05
  } else if(site_edge == "Rural Forest Edge"){
    metadata_leaf_16s$tukey[i] <- "ab"
    metadata_leaf_16s$location[i] <- max(metadata_leaf_16s$perc_methanotroph[which(metadata_leaf_16s$site_edge == site_edge)]) + 0.05
  } else if(site_edge == "Rural Forest Interior"){
    metadata_leaf_16s$tukey[i] <- "a"
    metadata_leaf_16s$location[i] <- max(metadata_leaf_16s$perc_methanotroph[which(metadata_leaf_16s$site_edge == site_edge)]) + 0.05
  }
}

# Reformat the x-axis tickmark labels
metadata_leaf_16s$site_edge <- gsub(" ", "\n", metadata_leaf_16s$site_edge)
metadata_leaf_16s$site_edge <- factor(metadata_leaf_16s$site_edge, 
                                      levels = c("Urban\nStreet\nTree", 
                                                 "Urban\nForest\nEdge",
                                                 "Urban\nForest\nInterior",
                                                 "Rural\nForest\nEdge",
                                                 "Rural\nForest\nInterior"))

# Visualize
fig3d <- ggplot() +
  stat_summary(data = metadata_leaf_16s, aes(x = site_edge, 
                                             y = perc_methanotroph, group = 1),
               fun = median,
               geom = "line",
               linetype = "dashed",
               linewidth = 0.8) +
  geom_boxplot(data = metadata_leaf_16s, outlier.shape = NA,
               aes(x=site_edge, y=perc_methanotroph, fill=site_edge), 
               alpha = 0.5) + 
  geom_jitter(data=metadata_leaf_16s, aes(x = site_edge, y = perc_methanotroph, 
                                          color = site_edge, alpha = 0.75, 
                                          stroke = NA)) +
  # this adds the Tukey groups labels defined above
  geom_text(data = metadata_leaf_16s, direction = "y",
            aes(x=site_edge, y=location, label = tukey), 
            size = 3) +
  theme_classic() + labs (y="log(Methanotrophic bacteria\nrelative abundance (%))", 
                          x="") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_markdown(size = 10, color = "black"), 
        legend.text = element_markdown(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_markdown(size = 10), legend.position = "none", 
        legend.background = element_rect(fill = NA)) +
  scale_fill_manual(values = categorical_pal) + 
  scale_color_manual(values = categorical_pal) + 
  guides(color = "none", shape = guide_legend(title = ""))+
  ggtitle("Methanotrophic bacteria in leaves")

fig3d
ggsave(paste0("atherton_figure3d_methanotrophicleaves_vs_treetype_boxplot_", 
              format(Sys.Date(), "%Y%m%d"),".pdf"), device = cairo_pdf, 
       plot = fig3d, dpi = 300, width = 5.75, height = 3.5, units = "in")
